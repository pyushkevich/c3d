#ifndef __itkSLICSuperVoxelImageFilter_hxx_
#define __itkSLICSuperVoxelImageFilter_hxx_

#include "itkSLICSuperVoxelImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include <algorithm>

namespace itk {


template <typename TInputImage, typename TLabelImage, typename TRealImage>
SLICSuperVoxelImageFilter<TInputImage, TLabelImage, TRealImage>
::SLICSuperVoxelImageFilter()
{
  m_MParameter = 20;
  m_SeedsPerDimension = 20;

  this->SetNumberOfRequiredInputs(2);

#if ITK_VERSION_MAJOR < 5
  m_Barrier = Barrier::New();
#endif
}

template <typename TInputImage, typename TLabelImage, typename TRealImage>
void
SLICSuperVoxelImageFilter<TInputImage, TLabelImage, TRealImage>
::SetGradientImage(TRealImage *image)
{
  m_GradientImage = image;
  this->SetNthInput(1, image);
}

template <typename TInputImage, typename TLabelImage, typename TRealImage>
void
SLICSuperVoxelImageFilter<TInputImage, TLabelImage, TRealImage>
#if ITK_VERSION_MAJOR >= 5
::GenerateData()
#else
::BeforeThreadedGenerateData(void)
#endif
{

#if ITK_VERSION_MAJOR >= 5
  typedef ImageRegionConstIteratorWithIndex<InputImageType> InputIter;
  typedef ImageRegionIterator<OutputImageType> OutputIter;
  typedef ImageRegionIterator<TRealImage> DistIter;
#endif

  // Input and output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename InputImageType::ConstPointer input = this->GetInput();

#if ITK_VERSION_MAJOR >= 5
  // Call AllocateOutputs is required
  this->AllocateOutputs();
#endif

  // Allocate the distance image
  m_DistanceImage = TRealImage::New();
  m_DistanceImage->SetRegions(output->GetBufferedRegion());
  m_DistanceImage->CopyInformation(output);
  m_DistanceImage->Allocate();
  m_DistanceImage->FillBuffer(1.0e100);

  // Initialize the label images
  output->FillBuffer(-1);

  // To sample the initial cluster centers, we create a dummy image of cluster
  // centers
  typename TLabelImage::Pointer imgCluster = TLabelImage::New();
  RegionType rgnOutput;
  for(int i = 0; i < ImageDimension; i++)
    {
    rgnOutput.SetSize(i, m_SeedsPerDimension);
    rgnOutput.SetIndex(i, 0);
    }

  imgCluster->SetRegions(rgnOutput);
  imgCluster->SetOrigin(input->GetOrigin());
  imgCluster->SetDirection(input->GetDirection());

  Vector<double,ImageDimension> newspc;

  m_SearchRegionMaxWidth = 0;
  for(int i = 0; i < ImageDimension; i++)
    {
    int dim = input->GetBufferedRegion().GetSize(i);
    double spc = input->GetSpacing()[i];
    double width = dim * spc;
    newspc[i] = width / m_SeedsPerDimension;

    int pixWidth = dim / m_SeedsPerDimension;
    m_SearchRegionOffset[i] = - pixWidth;
    m_SearchRegionSize[i] = 1 + 2 * pixWidth;

    if(m_SearchRegionMaxWidth < pixWidth)
      m_SearchRegionMaxWidth = pixWidth;
    }

  imgCluster->SetSpacing(newspc);
  imgCluster->SetPixelContainer(output->GetPixelContainer());

  // Create a neighborhood iterator for the gradient image search
  typedef itk::ConstNeighborhoodIterator<TRealImage> HoodIter;
  SizeType radius; radius.Fill(1);
  HoodIter itGrad(radius, m_GradientImage, m_GradientImage->GetBufferedRegion());

  // Create the cluster array
  m_Clusters.resize(imgCluster->GetBufferedRegion().GetNumberOfPixels());
  typedef ImageRegionIteratorWithIndex<TLabelImage> OutIter;
  int k = 0;
  for(OutIter it(imgCluster, imgCluster->GetBufferedRegion()); !it.IsAtEnd(); ++it, ++k)
    {
    // Map the cluster center into a voxel position
    itk::ContinuousIndex<double, ImageDimension> cidx;
    for(int i = 0; i < ImageDimension; i++)
      cidx[i] = it.GetIndex()[i] + 0.5;
    IndexType idxInput;
    PointType ptx;
    imgCluster->TransformContinuousIndexToPhysicalPoint(cidx, ptx);
    input->TransformPhysicalPointToIndex(ptx, idxInput);

    // Search the gradient magnitude image for the
    itGrad.SetLocation(idxInput);

    bool inbounds;
    int jBest = 0;
    RealPixelType gBest = 1e100;

    for(int j = 0; j < itGrad.Size(); j++)
      {
      RealPixelType gval = itGrad.GetPixel(j, inbounds);
      if(inbounds && gval < gBest)
        {
        gBest = gval; jBest = j;
        }
      }

    // Now we have finally an index
    IndexType idxCenter = itGrad.GetIndex(jBest);
    InputPixelType pxCenter = input->GetPixel(idxCenter);
    m_Clusters[k].Index = idxCenter;
    m_Clusters[k].Pixel = pxCenter;
    m_Clusters[k].Count = 1;
    m_Clusters[k].MaxPixelDist = m_MParameter;
    }

#if ITK_VERSION_MAJOR >= 5
  // A mutex for accelerating cluster info
  std::mutex cluster_mutex;

  // Run for some number of iterations
  for(int iter = 0; iter < 10; iter++)
    {
    // An value keeping track of whether the global cluster data structure is being updated for
    // the first time by a thread or not
    bool is_first_cluster_update = true;

    // Use the new ITK5 code for parallelization. We will compute clusters for each threaded
    // region, then integrate into centralized repository
    itk::MultiThreaderBase::Pointer mt = itk::MultiThreaderBase::New();
    mt->ParallelizeImageRegion<Self::ImageDimension>(
          this->GetOutput()->GetBufferedRegion(),
          [this,&input,&output,&cluster_mutex,&is_first_cluster_update](const RegionType &thread_region)
      {
      // Local copy of the clusters
      ClusterVector cv_local = m_Clusters;

      // Update voxel memberships
      for(int i = 0; i < m_Clusters.size(); i++)
        {
        // The cluster
        Cluster &cluster = m_Clusters[i];

        // Get the search region for this cluster
        RegionType region;
        region.SetSize(this->m_SearchRegionSize);
        region.SetIndex(cluster.Index + this->m_SearchRegionOffset);

        // Check if the region overlaps with the output region
        if(region.Crop(thread_region))
          {
          // Create an iterator for the cropped region
          InputIter itInput(input, region);
          OutputIter itOutput(output, region);
          DistIter itDist(m_DistanceImage, region);

          while(!itInput.IsAtEnd())
            {
            double dist = this->ComputeDistance(itInput.GetIndex(), itInput.Get(), cluster);
            if(itDist.Value() > dist)
              {
              itDist.Set(dist);
              itOutput.Set(i);
              }
            ++itInput; ++itOutput; ++itDist;
            }
          }

        // While we are in this loop, we clear the local cluster information
        cv_local[i].Reset();
        }

      // We have computed the membership of each voxel in our region. Now, we need to update the
      // cluster means.

      // Create an iterator for the cropped region
      InputIter itInput(input, thread_region);
      OutputIter itOutput(output, thread_region);

      while(!itInput.IsAtEnd())
        {
        int clid = itOutput.Value();
        InputPixelType pix = itInput.Get();
        Cluster &cl = cv_local[clid];
        for(int d = 0; d < ImageDimension; d++)
          cl.Index[d] += itInput.GetIndex()[d];

        if(pix > cl.MaxPixel)
          cl.MaxPixel = pix;
        if(pix < cl.MinPixel)
          cl.MinPixel = pix;

        cl.Pixel += pix;
        cl.Count++;

        ++itInput; ++itOutput;
        }

      // Now use the mutex lock to update the global cluster info
      std::lock_guard<std::mutex> guard(cluster_mutex);
      if(is_first_cluster_update)
        {
        m_Clusters = cv_local;
        is_first_cluster_update = false;
        }
      else
        {
        for(int i = 0; i < m_Clusters.size(); i++)
          {
          for(int d = 0; d < ImageDimension; d++)
            m_Clusters[i].Index[d] += cv_local[i].Index[d];
          m_Clusters[i].Pixel += cv_local[i].Pixel;
          m_Clusters[i].Count += cv_local[i].Count;
          m_Clusters[i].MinPixel = std::min(m_Clusters[i].MinPixel, cv_local[i].MinPixel);
          m_Clusters[i].MaxPixel = std::min(m_Clusters[i].MaxPixel, cv_local[i].MaxPixel);
          }
        }
      }, nullptr); // End of parallel block

    // Normalize the cluster centers
    for(int i = 0; i < m_Clusters.size(); i++)
      {
      for(int d = 0; d < ImageDimension; d++)
        m_Clusters[i].Index[d] = m_Clusters[i].Index[d] / m_Clusters[i].Count;
      m_Clusters[i].Pixel = m_Clusters[i].Pixel / m_Clusters[i].Count;
      m_Clusters[i].MaxPixelDist = m_Clusters[i].MaxPixel - m_Clusters[i].MinPixel;
      }
    }
}
#else
  // Initialize the per-thread cluster data
  m_PerThreadClusters.resize(this->GetNumberOfThreads(), m_Clusters);

  // Initialize the barrier
  m_Barrier->Initialize(this->GetNumberOfThreads());
}

template <typename TInputImage, typename TLabelImage, typename TRealImage>
void
SLICSuperVoxelImageFilter<TInputImage, TLabelImage, TRealImage>
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  typedef ImageRegionConstIteratorWithIndex<InputImageType> InputIter;
  typedef ImageRegionIterator<OutputImageType> OutputIter;
  typedef ImageRegionIterator<TRealImage> DistIter;

  // Input and output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename InputImageType::ConstPointer input = this->GetInput();

  // A thread-specific copy of the cluster centers
  ClusterVector &cvLocal = m_PerThreadClusters[threadId];

  // This is the iterative part of the filter

  // Standard size of the search neighborhood

  
  // Here there is a little divergence in behavior. The main thread regulates
  // the rest of the threads
  for(int iter = 0; iter < 10; iter++)
    {
    // Update voxel memberships
    for(int i = 0; i < m_Clusters.size(); i++)
      {
      // The cluster
      Cluster &cluster = m_Clusters[i];

      // Get the search region for this cluster
      RegionType region;
      region.SetSize(m_SearchRegionSize);
      region.SetIndex(cluster.Index + m_SearchRegionOffset);

      // Check if the region overlaps with the output region
      if(region.Crop(outputRegionForThread))
        {
        // Create an iterator for the cropped region
        InputIter itInput(input, region);
        OutputIter itOutput(output, region);
        DistIter itDist(m_DistanceImage, region);

        while(!itInput.IsAtEnd())
          {
          double dist = this->ComputeDistance(itInput.GetIndex(), itInput.Get(), cluster);
          if(itDist.Value() > dist)
            {
            itDist.Set(dist);
            itOutput.Set(i);
            }
          ++itInput; ++itOutput; ++itDist;
          }
        }

      // While we are in this loop, we clear the local cluster information
      cvLocal[i].Reset();
      }

    // We have computed the membership of each voxel in our region. Now, we need to update the
    // cluster means. However, this requires storing a separate copy of the cluster for each 
    // thread

    // Create an iterator for the cropped region
    InputIter itInput(input, outputRegionForThread);
    OutputIter itOutput(output, outputRegionForThread);

    while(!itInput.IsAtEnd())
      {
      int clid = itOutput.Value();
      InputPixelType pix = itInput.Get();
      Cluster &cl = cvLocal[clid];
      for(int d = 0; d < ImageDimension; d++)
        cl.Index[d] += itInput.GetIndex()[d];

      if(pix > cl.MaxPixel)
        cl.MaxPixel = pix;
      if(pix < cl.MinPixel)
        cl.MinPixel = pix;

      cl.Pixel += pix;
      cl.Count++;

      ++itInput; ++itOutput;
      }

    // Now we need to wait for all the threads to finish
    m_Barrier->Wait();

    // Now, only the first thread is going to combine the per-thread cluster means into
    // common cluster mean
    if(threadId == 0)
      {
      m_Clusters = cvLocal;
      for(int i = 1; i < m_Clusters.size(); i++)
        {
        for(int k = 0; k < m_PerThreadClusters.size(); k++)
          {
          for(int d = 0; d < ImageDimension; d++)
            m_Clusters[i].Index[d] += m_PerThreadClusters[k][i].Index[d];
          m_Clusters[i].Pixel += m_PerThreadClusters[k][i].Pixel;
          m_Clusters[i].Count += m_PerThreadClusters[k][i].Count;
          m_Clusters[i].MinPixel = 
            std::min(m_Clusters[i].MinPixel, m_PerThreadClusters[k][i].MinPixel);
          m_Clusters[i].MaxPixel = 
            std::min(m_Clusters[i].MaxPixel, m_PerThreadClusters[k][i].MaxPixel);
          }
        for(int d = 0; d < ImageDimension; d++)
          m_Clusters[i].Index[d] = m_Clusters[i].Index[d] / m_Clusters[i].Count;
        m_Clusters[i].Pixel = m_Clusters[i].Pixel / m_Clusters[i].Count;
        m_Clusters[i].MaxPixelDist = m_Clusters[i].MaxPixel - m_Clusters[i].MinPixel;
        }

      std::cout << "Cluster center 55: " << m_Clusters[55].Index << " " << m_Clusters[55].Pixel << std::endl;
      }

    // All threads must wait for thread 0 to get here
    m_Barrier->Wait();
    }
}
#endif

template <typename TInputImage, typename TLabelImage, typename TRealImage>
double
SLICSuperVoxelImageFilter<TInputImage, TLabelImage, TRealImage>
::ComputeDistance(const IndexType &index, const InputPixelType &pixel, Cluster &cluster)
{
  // Compute the squared spatial distance
  double ds2 = 0.0;
  for(int d = 0; d < ImageDimension; d++)
    {
    ds2 += (index[d] - cluster.Index[d]) * (index[d] - cluster.Index[d]);
    }

  double dc2 = (pixel - cluster.Pixel) * (pixel - cluster.Pixel);

  return (dc2 * m_SearchRegionMaxWidth * m_SearchRegionMaxWidth 
    + ds2 * m_MParameter * m_MParameter);

}

}

#endif

