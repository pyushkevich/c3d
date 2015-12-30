#ifndef RANDOMFORESTCLASSIFYIMAGEFILTER_TXX
#define RANDOMFORESTCLASSIFYIMAGEFILTER_TXX

#include "RandomForestClassifyImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "RandomForestClassifier.h"
#include "ImageCollectionToImageFilter.h"
#include <itkProgressReporter.h>

#include "Library/data.h"
#include "Library/forest.h"
#include "Library/statistics.h"
#include "Library/classifier.h"


template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::RandomForestClassifyImageFilter()
{
  // m_MixtureModel = NULL;
  m_GenerateClassProbabilities = false;
  m_Classifier = NULL;
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::~RandomForestClassifyImageFilter()
{
}


template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::AddScalarImage(InputImageType *image)
{
  this->AddInput(image);
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::AddVectorImage(InputVectorImageType *image)
{
  this->AddInput(image);
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::SetClassifier(ClassifierType *classifier)
{
  m_Classifier = classifier;
  this->Modified();
  this->UpdateOutputs();
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::SetGenerateClassProbabilities(bool flag)
{
  m_GenerateClassProbabilities = flag;
  this->Modified();
  this->UpdateOutputs();
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::UpdateOutputs()
{
  if(m_Classifier)
    {
    if(m_GenerateClassProbabilities)
      {
      this->SetNumberOfIndexedOutputs(m_Classifier->GetClassToLabelMapping().size());
      for(int i = 1; i < m_Classifier->GetClassToLabelMapping().size(); i++)
        this->SetNthOutput(i, this->MakeOutput(i));
      }
    else
      {
      this->SetNumberOfIndexedOutputs(1);
      }
    }
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::GenerateInputRequestedRegion()
{
  itk::ImageSource<TOutputImage>::GenerateInputRequestedRegion();

  for( itk::InputDataObjectIterator it( this ); !it.IsAtEnd(); it++ )
    {
    // Check whether the input is an image of the appropriate dimension
    InputImageType *input = dynamic_cast< InputImageType * >( it.GetInput() );
    InputVectorImageType *vecInput = dynamic_cast< InputVectorImageType * >( it.GetInput() );
    if (input)
      {
      InputImageRegionType inputRegion;
      this->CallCopyOutputRegionToInputRegion( inputRegion, this->GetOutput()->GetRequestedRegion() );
      inputRegion.PadByRadius(m_Classifier->GetPatchRadius());
      inputRegion.Crop(input->GetLargestPossibleRegion());
      input->SetRequestedRegion(inputRegion);
      }
    else if(vecInput)
      {
      InputImageRegionType inputRegion;
      this->CallCopyOutputRegionToInputRegion( inputRegion, this->GetOutput()->GetRequestedRegion() );
      inputRegion.PadByRadius(m_Classifier->GetPatchRadius());
      inputRegion.Crop(vecInput->GetLargestPossibleRegion());
      vecInput->SetRequestedRegion(inputRegion);
      }
    }
}


template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
  os << indent << "RandomForestClassifyImageFilter" << std::endl;
}

template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
void
RandomForestClassifyImageFilter<TInputImage, TInputVectorImage, TOutputImage, TLabel>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  assert(m_Classifier);

  // Number of outputs generated
  int nClass = m_Classifier->GetClassToLabelMapping().size();
  int nOut = m_GenerateClassProbabilities ? nClass : 1;

  // Primary output
  OutputImagePointer outputPtr = this->GetOutput(0);

  // Fill the output region with zeros in all outputs
  for(int k = 0; k < nOut; k++)
    {
    itk::ImageRegionIterator<OutputImageType> zit(this->GetOutput(k), outputRegionForThread);
    for(; !zit.IsAtEnd(); ++zit)
      zit.Set((OutputPixelType) 0);
    }

  // Adjust the output region so that we don't touch image boundaries.
  OutputImageRegionType crop_region = outputPtr->GetLargestPossibleRegion();
  crop_region.ShrinkByRadius(m_Classifier->GetPatchRadius());
  OutputImageRegionType out_region = outputRegionForThread;
  bool can_crop = out_region.Crop(crop_region);

  if(!can_crop)
    return;

  // Create output iterators
  typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutputIter;
  std::vector<OutputIter> it_out;
  for(int k = 0; k < nOut; k++)
    it_out.push_back(OutputIter(this->GetOutput(k), out_region));

  // Create a collection iterator for the inputs
  typedef ImageCollectionConstRegionIteratorWithIndex<
      TInputImage, TInputVectorImage> CollectionIter;

  // Configure the input collection iterator
  CollectionIter cit(out_region);
  for( itk::InputDataObjectIterator it( this ); !it.IsAtEnd(); it++ )
    cit.AddImage(it.GetInput());

  // TODO: This is hard-coded
  cit.SetRadius(m_Classifier->GetPatchRadius());

  // Get the number of components
  int nComp = cit.GetTotalComponents();
  int nPatch = cit.GetNeighborhoodSize();
  int nColumns = nComp * nPatch;

  // Are coordinate features used?
  if(m_Classifier->GetUseCoordinateFeatures())
    nColumns += 3;

  // Get the class weights (as they are assigned to foreground/background)
  const typename ClassifierType::WeightArray &class_weights = m_Classifier->GetClassWeights();

  // Create the MLdata representing each voxel (?)
  typedef Histogram<InputPixelType,LabelType> HistogramType;
  typedef MLData<InputPixelType,HistogramType *> TestingDataType;
  TestingDataType testData(1, nColumns);

  // Get the number of trees
  int nTrees = m_Classifier->GetForest()->trees_.size();

  // Create and allocate the test result vector
  typedef Vector<Vector<HistogramType *> > TestingResultType;
  TestingResultType testResult;
  testResult.Resize(nTrees);
  for(int i = 0; i < nTrees; i++)
    testResult[i].Resize(1);

  // Some vectors that are allocated for speed
  std::vector<size_t> vIndex(1);
  std::vector<bool> vResult(1);
  std::vector<double> vClassProb(nClass);

  // Iterate through all the voxels
  while(!it_out[0].IsAtEnd())
    {
    // Assign the data to the testData vector
    int k = 0;
    for(int i = 0; i < nComp; i++)
      for(int j = 0; j < nPatch; j++)
        testData.data[0][k++] = cit.NeighborValue(i,j);

    // Add the coordinate features
    if(m_Classifier->GetUseCoordinateFeatures())
      for(int d = 0; d < 3; d++)
        testData.data[0][k++] = it_out[0].GetIndex()[d];

    // Perform classification on this data
    m_Classifier->GetForest()->ApplyFast(testData, testResult, vIndex, vResult);

    // How we update outputs depends on the m_GenerateClassProbabilities flag
    if(m_GenerateClassProbabilities)
      {
      std::fill(vClassProb.begin(), vClassProb.end(), 0.0);
      double p_total = 0;
      for(int i = 0; i < testResult.Size(); i++)
        {
        HistogramType *hist = testResult[i][0];
        for(int j = 0; j < nClass; j++)
          {
          double p = hist->prob_[j];
          vClassProb[j] += p;
          p_total += p;
          }
        }

      for(int j = 0; j < nClass; j++)
        {
        if(p_total > 0)
          it_out[j].Set((OutputPixelType)(vClassProb[j] / p_total));
        }
      }
    else
      {
      // New code: compute output map with a bias parameter. The bias parameter q is such
      // that p_fore = q maps to 0 speed value. For the time being we just shift the linear
      // mapping from p_fore to speed and cap speed between -1 and 1

      // First we compute p_fore - for some reason not all trees in the forest have probabilities
      // summing up to one (some are zero), so we need to use division
      double p_fore_total = 0, p_total = 0;
      for(int i = 0; i < testResult.Size(); i++)
        {
        HistogramType *hist = testResult[i][0];
        for(int j = 0; j < nClass; j++)
          {
          double p = hist->prob_[j];
          if(class_weights[j] > 0.0)
            p_fore_total += p;
          p_total += p;
          }
        }

      // Set output only if the total probability is non-zero
      if(p_total > 0)
        {
        double q = m_Classifier->GetBiasParameter();
        double p_fore = p_fore_total / p_total;
        double speed = 2 * (p_fore - q);
        if(speed < -1.0)
          speed = -1.0;
        else if(speed > 1.0)
          speed = 1.0;

        it_out[0].Set((OutputPixelType)(speed * 0x7fff));
        }
      }

    ++cit;
    for(int k = 0; k < nOut; k++)
      ++it_out[k];
    }
}


#endif
