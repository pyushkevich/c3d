/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MRFVote.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include "MRFVote.h"
#include "itkImageRandomNonRepeatingIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNaryFunctorImageFilter.h"
#include "itkImageRegionIterator.h"

#include "GCoptimization.h"

/**
 * Functor for computing the mask
 */
template <class TPixel, class TMaskPixel = TPixel>
class MRFVoteMaskFunctor
{
public:
  typedef MRFVoteMaskFunctor<TPixel> Self;
  inline TPixel operator() (const std::vector<TPixel> &x)
    {
    // Test validity of the posterior
    for(int d = 0; d < x.size(); d++)
      if(!(x[d] >= 0 && x[d] <= 1))
        return 0;
    return 1;
    }

  TMaskPixel operator == (const Self &) const { return true; }
  TMaskPixel operator != (const Self &) const { return false; }
};


/**
 * The voting approach is as follows:
 *
 * Each input image is a `posterior` map for a given label. More
 * specifically, we consider each posterior to be the sum of votes
 * for that label. The votes do not have to add up to one, allowing
 * for missing data. Negative posterior values mean that the voxel
 * should be excluded from the voting (i.e. a way to supply a mask).
 */

template <class TPixel, unsigned int VDim>
void
MRFVote<TPixel, VDim>
::operator() (Mode mode, double beta)
{
  // Number of imaged voted over
  int nl = c->GetStackSize();

  // The first image
  ImagePointer iFirst  = c->PeekImage(0);

  // The mask computation filter
  typedef itk::Image<int, VDim> MaskImageType;
  typedef MRFVoteMaskFunctor<TPixel, int> MaskFunctor;
  typedef itk::NaryFunctorImageFilter<ImageType, MaskImageType, MaskFunctor> MaskFilter;
  typename MaskFilter::Pointer maskFilter = MaskFilter::New();
  for(int i = 0; i < nl; i++)
    maskFilter->SetInput(i, c->PeekImage(i));
  maskFilter->Update();
  typename MaskImageType::Pointer mask = maskFilter->GetOutput();

  // Create an image for looking up nodes
  typename MaskImageType::Pointer nodeImage = MaskImageType::New();
  nodeImage->CopyInformation(mask);
  nodeImage->SetRegions(mask->GetBufferedRegion());
  nodeImage->Allocate();

  // Fill the node index image
  typedef itk::ImageRegionIterator<MaskImageType> MaskIter;
  MaskIter itMask(mask, mask->GetBufferedRegion());
  MaskIter itNode(nodeImage, mask->GetBufferedRegion());

  int N = 0;
  while(!itNode.IsAtEnd())
    {
    if(itMask.Value() > 0)
      itNode.Set(++N);
    else
      itNode.Set(0);
    ++itNode; 
    ++itMask;
    }

  // Start reporting
  *c->verbose << "Performing MRF regularized voting using Graph Cuts" << std::endl;
  *c->verbose << "  Vertices: " << N << std::endl;
  *c->verbose << "  Labels:   " << nl << std::endl;

  // Must have some nodes!
  if(N == 0)
    throw ConvertException("No nodes for the MRF graph in -vote-mrf");

  // Iterator for enumerating all edges in the image (without duplication)
  itk::Size<VDim> radius; radius.Fill(1);
  typedef itk::ShapedNeighborhoodIterator<MaskImageType> ShapeIter;
  ShapeIter itShaped(radius, nodeImage, nodeImage->GetBufferedRegion());
  itShaped.ClearActiveList();

  // Set the allowed offsets for looking up neighbors. We only add offsets in one direction
  // to avoid duplicate edges in the graph.
  for(int d = 0; d < VDim; d++)
    {
    itk::Offset<VDim> offset;
    offset.Fill(0);
    offset[d] = -1;
    itShaped.ActivateOffset(offset);

    // Both direction edges?
    offset.Fill(0);
    offset[d] = 1;
    itShaped.ActivateOffset(offset);
    }

  // Create iterators for the lihelihood images
  std::vector<Iterator> lit;
  for(int i = 0; i < nl; i++)
    lit.push_back(Iterator(c->PeekImage(i), mask->GetBufferedRegion()));

  // Create the graph
  GCoptimizationGeneralGraph graph(N, nl);

  // Temporary storage
  double *pvec = new double[nl];

  // Iterate over the mask region, creating a graph
  int n_edges = 0;
  while(!itShaped.IsAtEnd())
    {
    // The central pixel must be non-zero in the mask
    int node_id = itShaped.GetCenterPixel();
    if(node_id > 0)
      {
      // Check nodes in the neighborhood
      for(typename ShapeIter::Iterator itn = itShaped.Begin(); !itn.IsAtEnd(); ++itn)
        {
        bool in_bounds;
        int nbr_id = itShaped.GetPixel(itn.GetNeighborhoodIndex(), in_bounds);
        if(in_bounds && nbr_id > 0)
          {
          // We use zero-based indexing for the graph class, 1-based in the image
          graph.setNeighbors(node_id - 1, nbr_id - 1);
          n_edges++;
          }
        }

      // Compute the energy term 
      double psum = 0.0;
      for(int i = 0; i < nl; i++)
        {
        pvec[i] = lit[i].Value();
        psum += pvec[i];
        }

      // Assign the energy terms using (psum - pvec[i]) formula. This assigns the
      // cost based on the number of 'votes against'
      for(int i = 0; i < nl; i++)
        {
        if(mode == VOTES_AGAINST)
          {
          graph.setDataCost(node_id - 1, i, psum - pvec[i]);
          }
        else
          {
          graph.setDataCost(node_id - 1, i, -log(pvec[i]));
          }
        }
      }

    ++itShaped;
    for(int i = 0; i < nl; i++)
      ++lit[i];
    }

  // Set the label costs
  for(int i = 0; i < nl; i++)
    for(int j = 0; j < nl; j++)
      graph.setSmoothCost(i, j, (i == j) ? 0.0 : beta);

  // Print the number of edges
  *c->verbose << "  Edges:    " << n_edges << std::endl;

  // Random order seems to make the most sense
  graph.setLabelOrder(true);

  // Run the alpha-expansion - with some verbose reporting
  double e_last = graph.compute_energy();
  c->PrintF("  Initial Energy:    Data = %9.4f   Smooth = %9.4f   Total = %9.4f\n",
    graph.giveDataEnergy(), graph.giveSmoothEnergy(), e_last);

  // List of labels in the graph
  std::vector<int> labels;
  for(int i = 0; i < nl; i++)
    labels.push_back(i);
 
  for(int iter = 0; iter < 20; iter++)
    {
    // Perform expansion in random order
    std::random_shuffle(labels.begin(), labels.end());

    // Perform alpha-expansion for each label
    for(int i = 0; i < labels.size(); i++)
      {
      graph.alpha_expansion(labels[i]);
      } 

    double e_now = graph.compute_energy();
    c->PrintF("  Iteration %4d:    Data = %9.4f   Smooth = %9.4f   Total = %9.4f\n",
      iter, graph.giveDataEnergy(), graph.giveSmoothEnergy(), e_now);

    if(e_now == e_last)
      break;

    e_last = e_now;
    }

  // Collect the final labeling
  ImagePointer result = ImageType::New();
  result->CopyInformation(mask);
  result->SetRegions(nodeImage->GetBufferedRegion());
  result->Allocate();
  result->FillBuffer(0.0);

  itk::ImageRegionConstIterator<MaskImageType> itPost(nodeImage, nodeImage->GetBufferedRegion());
  Iterator itResult(result, result->GetBufferedRegion());
  while(!itResult.IsAtEnd())
    {
    int node_id = itPost.Value();
    if(node_id > 0)
      itResult.Set(graph.whatLabel(node_id - 1) + 1);
    ++itPost; ++itResult;
    }

  // Clear histogram
  delete[] pvec;

  // Put result on stack
  c->m_ImageStack.clear();
  c->PushImage(result);
}

// Invocations
template class MRFVote<double, 2>;
template class MRFVote<double, 3>;
template class MRFVote<double, 4>;
