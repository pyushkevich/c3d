/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WriteImage.cxx
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

#include "WriteImage.h"
#include "itksys/SystemTools.hxx"
#include "itkImageFileWriter.h"
#include "itkIOCommon.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template<class TPixel, unsigned int VDim>
template<class TOutPixel>
void 
WriteImage<TPixel, VDim>
::TemplatedWriteImage(const char *file, double xRoundFactor, int pos)
{
  // Get the input image
  if(c->m_ImageStack.size() == 0)
    throw ConvertException("No data has been generated! Can't write to %s", file);

  // Get the image at the given position
  if(pos < 0) pos = c->m_ImageStack.size() - 1;
  ImagePointer input = c->m_ImageStack[pos];
  
  // Create the output image 
  typedef itk::OrientedRASImage<TOutPixel, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetDirection(input->GetDirection());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Describe what we are doing
  *c->verbose << "Writing #" << pos+1 << " to file " << file << endl;
  *c->verbose << "  Output voxel type: " << c->m_TypeId << "[" << typeid(TOutPixel).name() << "]" << endl;
  *c->verbose << "  Rounding off: " << (xRoundFactor == 0.0 ? "Disabled" : "Enabled") << endl;
  
  // Set the SPM originator header
  MakeSPMOriginFix(input, output);

  // Copy everything, rounding if the pixel type is integer
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = (TOutPixel) (input->GetBufferPointer()[i] + xRoundFactor);

  // Set the file notes for this image
  itk::EncapsulateMetaData<string>(
    output->GetMetaDataDictionary(),itk::ITK_FileNotes,
        std::string("Created by Convert3D"));

  // Write the image out
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(output);
  writer->SetFileName(file);
  writer->SetUseCompression(c->m_UseCompression);
  try { writer->Update(); }
  catch (itk::ExceptionObject &exc) {
    cerr << "Error writing image to " << file << endl;
    cerr << "ITK Exception: " << exc << endl;
    throw -1;
  }
}

template<class TPixel, unsigned int VDim>
void 
WriteImage<TPixel, VDim>
::MakeSPMOriginFix(itk::ImageBase<VDim> *input, itk::ImageBase<VDim> *output)
{
  // Set the SPM originator header
  if(c->m_FlagSPM)
    {
    size_t i;
    string originator;
    originator.resize(VDim * 2);
    
    // Compute the SPM-style origin of the image
    *c->verbose << "  Setting SPM origin field to:";
    for(i = 0; i < VDim; i++)
      {
      short ospm = (short)(0.5 - input->GetOrigin()[i] / input->GetSpacing()[i]);
      originator[2*i] = (char)(ospm & 0x00ff);
      originator[2*i+1] = (char)(ospm >> 8);
      *c->verbose << ospm << " ";
      }
    originator[2*i] = 0;
    *c->verbose << endl;

    itk::EncapsulateMetaData<string>(
      output->GetMetaDataDictionary(),itk::ITK_FileOriginator,originator);
    }
}


template<class TPixel, unsigned int VDim>
template<class TOutPixel>
void 
WriteImage<TPixel, VDim>
::TemplatedWriteMultiComponentImage(const char *file, double xRoundFactor, int pstart, int ncomp)
{
  if(ncomp <= 0)
    throw ConvertException("No data has been generated! Can't write to %s", file);

  // Get the top image on the stack (for reference information)
  ImagePointer itop = c->m_ImageStack.back();

  // Check compatibility
  for(size_t i = 0; i < ncomp-1; i++)
    if(c->m_ImageStack[pstart+i]->GetBufferedRegion().GetSize() != 
      itop->GetBufferedRegion().GetSize())
      {
      throw ConvertException("Multicomponent output error: mismatch in image dimensions");
      }

  // Define the output image type
  typedef itk::VectorImage<TOutPixel, VDim> OutputImageType;
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->CopyInformation(itop);
  output->SetRegions(itop->GetBufferedRegion());
  output->SetNumberOfComponentsPerPixel(ncomp);
  output->Allocate();

  // Describe what we are doing
  *c->verbose << "Writing Images " << pstart+1 << " to " << (pstart+ncomp) << " to multicomponent file " << file << endl;
  *c->verbose << "  Output voxel type: " << c->m_TypeId << "[" << typeid(TOutPixel).name() << "]" << endl;
  *c->verbose << "  Rounding off: " << (xRoundFactor == 0.0 ? "Disabled" : "Enabled") << endl;
    
  // Set the SPM originator header
  MakeSPMOriginFix(itop, output);

  // Copy everything, rounding if the pixel type is integer
  size_t n = itop->GetBufferedRegion().GetNumberOfPixels();
  for(size_t j = 0; j < ncomp; j++)
    {
    TPixel *buf = c->m_ImageStack[pstart+j]->GetBufferPointer();
    TOutPixel *out = output->GetBufferPointer() + j;
    for(size_t i = 0; i < n; i++, buf++, out+=ncomp)
      *out = (TOutPixel) (*buf + xRoundFactor);
    }

  // Write the image out
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(output);
  writer->SetFileName(file);
  writer->SetUseCompression(c->m_UseCompression);
  try { writer->Update(); }
  catch (itk::ExceptionObject &exc) 
    {
    throw ConvertException("Error writing image to %s\nITK Exception: %s", file, exc.GetDescription());
    }
}

template <class TPixel, unsigned int VDim>
void
WriteImage<TPixel, VDim>
::WriteMultiComponent(const char *file, int ncomp, int pos)
{
  // Get the position of the first image to include
  if(pos < 0)
    pos = c->m_ImageStack.size() - ncomp;

  if(c->m_TypeId == "char" || c->m_TypeId == "byte")
    TemplatedWriteMultiComponentImage<char>(file, c->m_RoundFactor, pos, ncomp);
  if(c->m_TypeId == "uchar" || c->m_TypeId == "ubyte")
    TemplatedWriteMultiComponentImage<unsigned char>(file, c->m_RoundFactor, pos, ncomp);
  
  if(c->m_TypeId == "short") 
    TemplatedWriteMultiComponentImage<short>(file, c->m_RoundFactor, pos, ncomp);
  if(c->m_TypeId == "ushort")
    TemplatedWriteMultiComponentImage<unsigned short>(file, c->m_RoundFactor, pos, ncomp);

  if(c->m_TypeId == "int") 
    TemplatedWriteMultiComponentImage<int>(file, c->m_RoundFactor, pos, ncomp);
  if(c->m_TypeId == "uint")
    TemplatedWriteMultiComponentImage<unsigned int>(file, c->m_RoundFactor, pos, ncomp);

  if(c->m_TypeId == "float") 
    TemplatedWriteMultiComponentImage<float>(file, 0.0, pos, ncomp);
  if(c->m_TypeId == "double")
    TemplatedWriteMultiComponentImage<double>(file, 0.0, pos, ncomp);
}


template <class TPixel, unsigned int VDim>
void
WriteImage<TPixel, VDim>
::operator() (const char *file, bool force, int pos)
{
  // Unless in 'force' mode, check if the image already exists
  if(!force && itksys::SystemTools::FileExists(file))
    {
    cerr << "File " << file << " already exists. Use -o option to override!" << endl;
    throw -1;
    }

  if(c->m_TypeId == "char" || c->m_TypeId == "byte")
    TemplatedWriteImage<char>(file, c->m_RoundFactor, pos);
  if(c->m_TypeId == "uchar" || c->m_TypeId == "ubyte")
    TemplatedWriteImage<unsigned char>(file, c->m_RoundFactor, pos);
  
  if(c->m_TypeId == "short") 
    TemplatedWriteImage<short>(file, c->m_RoundFactor, pos);
  if(c->m_TypeId == "ushort")
    TemplatedWriteImage<unsigned short>(file, c->m_RoundFactor, pos);

  if(c->m_TypeId == "int") 
    TemplatedWriteImage<int>(file, c->m_RoundFactor, pos);
  if(c->m_TypeId == "uint")
    TemplatedWriteImage<unsigned int>(file, c->m_RoundFactor, pos);

  if(c->m_TypeId == "float") 
    TemplatedWriteImage<float>(file, 0.0, pos);
  if(c->m_TypeId == "double")
    TemplatedWriteImage<double>(file, 0.0, pos);
}



// Invocations
template class WriteImage<double, 2>;
template class WriteImage<double, 3>;
template class WriteImage<double, 4>;
