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
#include "itkNeighborhoodIterator.h"


template <class TPixel, unsigned int VDim>
void
MRFVote<TPixel, VDim>
::operator() (double beta, size_t niter, bool flagUseSplitLabelSet)
{
  // Create a maximum image
  ImagePointer i0 = c->m_ImageStack[0];

  // Number of imaged voted over
  size_t nl = c->m_ImageStack.size();

  // If this is in response to a split command, retrieve the label set
  typename Converter::LabelSet lset;
  if(flagUseSplitLabelSet)
    {
    if(nl != c->GetSplitLabelSet().size())
      throw ConvertException(
        "Merge failed: number of images on the stack (%i) different "
        "from the number of split labels (%i)", 
        nl, c->GetSplitLabelSet().size());
    lset = c->GetSplitLabelSet();
    }

  // Otherwise, the label mapping is identity
  else
    {
    for(size_t i = 0; i < nl; i++)
      lset.push_back((double) i);
    }

  // Create a vote image (short for now)
  ImagePointer ilabel = ImageType::New();
  ilabel->SetRegions(i0->GetBufferedRegion());
  ilabel->SetOrigin(i0->GetOrigin());
  ilabel->SetSpacing(i0->GetSpacing());
  ilabel->SetDirection(i0->GetDirection());
  ilabel->Allocate();
  ilabel->FillBuffer(lset[0]);

  // Say something
  *c->verbose << "Collapsing " << c->m_ImageStack.size() << 
    " images into a multi-label image via maximum voting" << endl;

  // For each of the images, update the vote
  for(size_t j = 1; j < nl; j++)
    {
    // Get the next image pointer
    ImagePointer ij = c->m_ImageStack[j];

    // Check the image dimensions
    if(ij->GetBufferedRegion() != ilabel->GetBufferedRegion())
      throw ConvertException("All voting images must have same dimensions");

    // Pairwise voting
    size_t n = ilabel->GetBufferedRegion().GetNumberOfPixels();
    for(size_t k = 0; k < n; k++)
      {
      double ibest = i0->GetBufferPointer()[k];
      if(ij->GetBufferPointer()[k] > ibest)
        {
        ibest = ij->GetBufferPointer()[k];
        ilabel->GetBufferPointer()[k] = j;
        }
      }
    }

  // The result of the voting gives us the initial guess. Now perform
  // ICM (Besag 1986, http://www.stat.duke.edu/~scs/Courses/Stat376/Papers/GibbsFieldEst/BesagDirtyPicsJRSSB1986.pdf)
  
  // Define an inner region (no boundary voxels)
  RegionType r_inner = ilabel->GetBufferedRegion();
  for(size_t d = 0; d < VDim; d++)
    {
    r_inner.SetIndex(d, r_inner.GetIndex(d)+1);
    r_inner.SetSize(d, r_inner.GetSize(d)-2);
    }

  // Allocate array for neighborhood histogram
  double *nhist = new double[nl];

  // Create neighborhood iterator
  typedef itk::NeighborhoodIterator<ImageType> NeighborIterType;
  typename NeighborIterType::RadiusType radius;
  radius.Fill(1);
  NeighborIterType nit(radius, ilabel, r_inner);
  nit.SetNeedToUseBoundaryCondition(false);

  // Do some iterations
  for(size_t q = 0; q < niter; q++)
    {
    // Keep track of number of updates
    size_t n_upd = 0;
    
    // Iterate over the inner region
    typedef itk::ImageRandomNonRepeatingIteratorWithIndex<ImageType> RandIterType;
    RandIterType rit(ilabel, r_inner);
    rit.SetNumberOfSamples(r_inner.GetNumberOfPixels());
    for(; !rit.IsAtEnd(); ++rit)
      {
      // Current pixel value
      TPixel x_i = rit.Value();

      // Clear the neighborhood histogram
      for(size_t j = 0; j < nl; j++)
        nhist[j] = 0.0;

      // Make up for the fact that the current voxel will be counted
      nhist[(int)x_i] = -1;
      
      // Iterate the neighborhood
      nit.SetLocation(rit.GetIndex()); // Slow , change later
      for(size_t k = 0; k < nit.Size(); k++)
        nhist[(int)nit.GetPixel(k)]++;

      // For each candidate label value, compute conditional posterior
      int j_best = x_i;
      double post_best = 1e100;
      for(int j = 0; j < nl; j++)
        {
        double p = c->m_ImageStack[j]->GetPixel(rit.GetIndex());
        if(p > 0)
          {
          double likelihood = -log(p);
          double prior = - beta * nhist[j];
          double post = likelihood + prior;
          if(post_best > post)
            { post_best = post; j_best=j; }
          }
        }
      if(x_i != j_best)
        {
        rit.Set(j_best);
        n_upd++;
        }
      }
    if(n_upd == 0)
      {
      *c->verbose << "  Early convergence after " << q << " iterations" << endl;
      break;
      }
    }

  // Relabel the image using /split labels
  if(flagUseSplitLabelSet)
    for(size_t k = 0; k < ilabel->GetBufferedRegion().GetNumberOfPixels(); k++)
      ilabel->GetBufferPointer()[k] = lset[(int) ilabel->GetBufferPointer()[k]];
  
  // Clear histogram
  delete[] nhist;

  // Put result on stack
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(ilabel);
}

// Invocations
template class MRFVote<double, 2>;
template class MRFVote<double, 3>;
template class MRFVote<double, 4>;
