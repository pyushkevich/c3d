/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ExportPatches.cxx
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

#include "ExportPatches.h"
#include "itkNeighborhoodIterator.h"
#include <vnl/vnl_random.h>

template <class TPixel, unsigned int VDim>
void
ExportPatches<TPixel, VDim>
::operator() (const char *out_file, const SizeType &radius, double sample_frequency)
{
  // The last image is the mask - used to decide which voxels to sample
  ImagePointer mask = c->PopImage();

  // Number of channels to export
  int n_chan = c->GetStackSize();

  // Create a neighborhood iterator for each channel
  typedef itk::NeighborhoodIterator<ImageType> NIterType;
  NIterType **itn = new NIterType*[n_chan];
  for(int i = 0; i < n_chan; i++)
    itn[i] = new NIterType(radius, c->PeekImage(i), mask->GetBufferedRegion());

  // Size of neighborhood
  int nb_size = itn[0]->Size();

  // Create the sample vector
  std::vector<double> sample_vec(n_chan * nb_size);

  // Open a file for writing
  FILE *f = fopen(out_file, "wb");

  // Iterate over the mask
  vnl_random rand;
  for(Iterator iter(mask, mask->GetBufferedRegion()); !iter.IsAtEnd(); ++iter)
    {
    if(iter.Value() != c->m_Background)
      {
      double drand = rand.drand32(sample_frequency);
      if(drand <= 1.0)
        {
        // Sample at this location
        int p = 0;
        for(int k = 0; k < n_chan; k++)
          {
          NIterType *it_curr = itn[k];
          for(int j = 0; j < nb_size; j++)
            sample_vec[p++] = it_curr->GetPixel(j);
          }

        // Write sample to file
        for(int m = 0; m < sample_vec.size()-1; m++)
          fprintf(f, "%f,",sample_vec[m]);
        fprintf(f, "%f\n",sample_vec.back());
        }
      }
    for(int k = 0; k < n_chan; k++)
      ++(*itn[k]);
    }

  // Restore the mask
  c->PushImage(mask);
}

// Invocations
template class ExportPatches<double, 2>;
template class ExportPatches<double, 3>;
template class ExportPatches<double, 4>;
