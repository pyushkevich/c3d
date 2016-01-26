/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    DicomSeriesList.cxx
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

#include "gdcmTag.h"
#include "gdcmFile.h"
#include "gdcmReader.h"
#include "gdcmStringFilter.h"

#include "DicomSeriesList.h"
#include "itkIOCommon.h"
#include <itksys/SystemTools.hxx>
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <set>


template <class TPixel, unsigned int VDim>
void
DicomSeriesList<TPixel, VDim>
::operator() (const char *dicom_dir)
{
  const gdcm::Tag tagRows(0x0028, 0x0010);
  const gdcm::Tag tagCols(0x0028, 0x0011);
  const gdcm::Tag tagDesc(0x0008, 0x103e);
  const gdcm::Tag tagTextDesc(0x0028, 0x0010);
  const gdcm::Tag tagSeriesInstanceUID(0x0020,0x000E);
  const gdcm::Tag tagSeriesNumber(0x0020,0x0011);
  const gdcm::Tag tagAcquisitionNumber(0x0020,0x0012);
  const gdcm::Tag tagInstanceNumber(0x0020,0x0013);

  // Get the directory where to search for the series
  std::string series_dir = dicom_dir;
  if(!itksys::SystemTools::FileIsDirectory(dicom_dir))
    series_dir = itksys::SystemTools::GetParentDirectory(dicom_dir);

  // Use the ITK stuff for parsing
  typename itk::GDCMSeriesFileNames::Pointer gdcm_series = itk::GDCMSeriesFileNames::New();
  gdcm_series->SetUseSeriesDetails(true);
  gdcm_series->SetDirectory(series_dir);

  std::cout
      << "SeriesNumber"
      << "\tDimensions"
      << "\tNumImages"
      << "\tSeriesDescription"
      << "\tSeriesID"
      << std::endl;

  // List all the unique series ids
  const itk::SerieUIDContainer uids = gdcm_series->GetSeriesUIDs();
  for(int i = 0; i < uids.size(); i++)
    {
    // Get the filenames for this serie
    const itk::FilenamesContainer &fc = gdcm_series->GetFileNames(uids[i]);
    if(!fc.size())
      continue;

    // Print out each series in order

    // Get tags for this file
    std::set<gdcm::Tag> tagset;
    tagset.insert(tagRows);
    tagset.insert(tagCols);
    tagset.insert(tagSeriesNumber);
    tagset.insert(tagDesc);

    // Read the tags
    gdcm::Reader reader;
    reader.SetFileName(fc.front().c_str());
    bool read = false;

    try { read = reader.ReadSelectedTags(tagset); }
    catch(...) { read = false; }

    if(read)
      {
      gdcm::StringFilter sf;
      sf.SetFile(reader.GetFile());

      // Read series description
      std::cout << sf.ToString(tagSeriesNumber);

      // Read the dimensions
      ostringstream oss;
      oss << sf.ToString(tagRows) << "x"
          << sf.ToString(tagCols) << "x"
          << fc.size();

      std::cout << "\t" << oss.str();
      std::cout << "\t" << fc.size();
      std::cout << "\t" << sf.ToString(tagDesc);
      std::cout << "\t" << uids[i];
      }

    std::cout << std::endl;
    }
}

// Invocations
template class DicomSeriesList<double, 2>;
template class DicomSeriesList<double, 3>;
template class DicomSeriesList<double, 4>;
