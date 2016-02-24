/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NormalizeLocalWindow.cxx
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
#ifndef __VectorImageTools_h_
#define __VectorImageTools_h_

#include <itkVectorImage.h>

template <class TInputImage>
typename itk::VectorImage<typename TInputImage::PixelType::ValueType, TInputImage::ImageDimension>::Pointer
WrapImageOfVectorAsVectorImage(TInputImage *image)
{
  typedef typename TInputImage::PixelType VectorType;
  typedef typename VectorType::ValueType ValueType;
  typedef itk::VectorImage<ValueType, TInputImage::ImageDimension> VectorImageType;

  // Masquerade as a vector image
  typename VectorImageType::Pointer vecImage = VectorImageType::New();
  vecImage->CopyInformation(image);
  vecImage->SetRegions(image->GetBufferedRegion());
  vecImage->SetNumberOfComponentsPerPixel(VectorType::Dimension);

  // Transfer the pixel container without copying memory
  typedef typename VectorImageType::PixelContainer OutPixelContainer;
  typename OutPixelContainer::Pointer pcon = OutPixelContainer::New();
  pcon->SetImportPointer(
        reinterpret_cast<ValueType *>(image->GetPixelContainer()->GetImportPointer()),
        image->GetPixelContainer()->Size() * VectorType::Dimension);
  pcon->SetContainerManageMemory(false);
  vecImage->SetPixelContainer(pcon);

  // Return the image
  return vecImage;
}



#endif
