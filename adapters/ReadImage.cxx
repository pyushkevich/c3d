/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ReadImage.cxx
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

#include "ReadImage.h"
#include "itkIOCommon.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"

template <class TPixel, unsigned int VDim>
void
ReadImage<TPixel, VDim>
::operator() (const char *file, const ImageInfo &info)
{
  // IO Base object
  typename itk::ImageIOBase::Pointer iobase;

  // Special handling for DICOM series
  if(info.dicom_series_id)
    {
    // Get the directory where to search for the series
    std::string series_dir = file;
    if(!itksys::SystemTools::FileIsDirectory(file))
      series_dir = itksys::SystemTools::GetParentDirectory(file);

    // NOTE: for the time being, we are relying on GDCMSeriesFileNames for
    // proper sorting of the DICOM data. This is marked as deprecated in GDCM 2
    typename itk::GDCMSeriesFileNames::Pointer gdcm_series = itk::GDCMSeriesFileNames::New();
    gdcm_series->SetUseSeriesDetails(true);
    gdcm_series->SetDirectory(series_dir);

    // Use the series provided by the user
    std::vector<std::string> dicom_files = gdcm_series->GetFileNames(info.dicom_series_id);

    // Read the information from the first filename
    if(dicom_files.size() == 0)
      throw ConvertException("Error: DICOM series not found. "
                             "Directory '%s' does not appear to contain a "
                             "series of DICOM images.", series_dir.c_str());

    // Report
    *c->verbose << "Reading #" << (1 + c->m_ImageStack.size()) << " from DICOM series "
                << info.dicom_series_id << " in " << series_dir << endl;

    iobase = itk::GDCMImageIO::New();
    iobase->SetFileName(dicom_files[0]);

    // Read the image information
    iobase->ReadImageInformation();

    // Create an image series reader
    typedef itk::ImageSeriesReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    // Set the filenames and read
    reader->SetFileNames(dicom_files);
    reader->SetImageIO(iobase);

    // Try reading this file
    try { reader->Update(); }
    catch(itk::ExceptionObject &exc)
      {
      throw ConvertException("Error reading DICOM series %s\nITK exception: %s", file, exc.GetDescription());
      }

    ImagePointer image = reader->GetOutput();
    c->m_ImageStack.push_back(image);

    // Copy the metadata from the first scan in the series
    /*
    const typename ReaderType::DictionaryArrayType *darr =
      reader->GetMetaDataDictionaryArray();
    if(darr->size() > 0)
      m_NativeImage->SetMetaDataDictionary(*((*darr)[0]));
      */
    }
  else
    {
    // Report
    *c->verbose << "Reading #" << (1 + c->m_ImageStack.size()) << " from " << file << endl;

    // Delegate to the Factory (based on filename, etc)
    iobase = itk::ImageIOFactory::CreateImageIO(file, itk::ImageIOFactory::ReadMode);

    if(!iobase)
      throw ConvertException("Unable to read image %s; IO factory can not create IO object.", file);

    iobase->SetFileName(file);

    // Read the image information
    iobase->ReadImageInformation();

    // Handle SPM origin if SPM flag is specified
    // TODO: get rid of this code, nifti is more reliable
    string ext = itksys::SystemTools::GetFilenameExtension(file);
    if((ext == ".hdr" || ext == ".img.gz" || ext == ".img") && c->m_FlagSPM)
      {
      string temp;
      if(itk::ExposeMetaData<std::string>(
          iobase->GetMetaDataDictionary(), itk::ITK_FileOriginator, temp))
        {
        // Read the SPM-style origin
        *c->verbose << "  Applying SPM origin :";
        for(size_t i=0; i < VDim; i++)
          {
          double sitk = iobase->GetSpacing(i);
          short xspm = (temp[2*i+1] << 8) + temp[2*i];
          double oitk = -sitk * xspm;
          *c->verbose << xspm << " ";
          iobase->SetOrigin(i, oitk);
          }
        *c->verbose << endl;
        }
      }

    // If the image has multiple components, we need to handle it specially
    if(iobase->GetNumberOfComponents() > 1 && c->m_MultiComponentSplit)
      {
      // Read the multi-component image
      typedef itk::VectorImage<TPixel, VDim> VectorImageType;
      typedef itk::ImageFileReader<VectorImageType> VectorReader;
      typename VectorReader::Pointer reader = VectorReader::New();
      reader->SetFileName(file);
      reader->SetImageIO(iobase);
      try { reader->Update(); }
      catch(itk::ExceptionObject &exc)
        {
        throw ConvertException("Error reading image %s\nITK exception: %s", file, exc.GetDescription());
        }

      // Report
      *c->verbose << "  Splitting " << iobase->GetNumberOfComponents() << "-component image." << endl;

      // Split the vector image into component images
      typename VectorImageType::Pointer vec = reader->GetOutput();
      size_t ncomp = vec->GetVectorLength();
      for(size_t i = 0; i < ncomp; i++)
        {
        ImagePointer icomp = ImageType::New();
        icomp->CopyInformation(vec);
        icomp->SetRegions(vec->GetBufferedRegion());
        icomp->Allocate();
        TPixel *src = vec->GetBufferPointer() + i;
        TPixel *dst = icomp->GetBufferPointer();
        TPixel *end = dst + vec->GetBufferedRegion().GetNumberOfPixels();
        for(; dst < end; dst++, src+=ncomp) *dst = *src;

        // Push the image on the stack
        c->m_ImageStack.push_back(icomp);
        }
      }
    else
      {
      // Set up the reader
      typedef itk::ImageFileReader<ImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(file);
      reader->SetImageIO(iobase);

      // Try reading this file
      try { reader->Update(); }
      catch(itk::ExceptionObject &exc)
        {
        throw ConvertException("Error reading image %s\nITK exception: %s", file, exc.GetDescription());
        }

      ImagePointer image = reader->GetOutput();
      c->m_ImageStack.push_back(image);
      }
    }
}

// Invocations
template class ReadImage<double, 2>;
template class ReadImage<double, 3>;
template class ReadImage<double, 4>;
