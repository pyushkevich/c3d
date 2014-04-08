/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    UpdateMetadataKey.cxx
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

#include "UpdateMetadataKey.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template <class TPixel, unsigned int VDim>
void
UpdateMetadataKey<TPixel, VDim>
::operator() (const char *key, const char *value)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Report what we are doing
  *c->verbose << "Updating metadata in #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Setting key " << key << " to value " << value << endl;

  // Update the metadata values
  itk::MetaDataDictionary &mdd = img->GetMetaDataDictionary();
  typedef itk::MetaDataObject<string> StringMetaData;
  typename StringMetaData::Pointer mdval = StringMetaData::New();
  mdval->SetMetaDataObjectValue(value);
  mdd[key] = mdval;
}

// Invocations
template class UpdateMetadataKey<double, 2>;
template class UpdateMetadataKey<double, 3>;
