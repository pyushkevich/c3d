/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    PrintImageInfo.cxx
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

#include "PrintImageInfo.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkSpatialOrientation.h"
#include <nifti1_io.h>

template<class AnyType>
bool
try_print_metadata(ostream &sout, itk::MetaDataDictionary &mdd, string key, AnyType deflt)
  {
  AnyType value = deflt;
  if(itk::ExposeMetaData<AnyType>(mdd, key, value))
    {
    sout << "    " << key << " = " << value << endl;
    return true;
    }
  else return false;
  }

string
get_rai_code(itk::SpatialOrientation::ValidCoordinateOrientationFlags code)
  {
  std::map<itk::SpatialOrientation::ValidCoordinateOrientationFlags, string> m_CodeToString;
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP] = "RIP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP] = "LIP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP] = "RSP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP] = "LSP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA] = "RIA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA] = "LIA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA] = "RSA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA] = "LSA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP] = "IRP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP] = "ILP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP] = "SRP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP] = "SLP";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA] = "IRA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA] = "ILA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA] = "SRA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA] = "SLA";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI] = "RPI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI] = "LPI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI] = "RAI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI] = "LAI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS] = "RPS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS] = "LPS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS] = "RAS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS] = "LAS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI] = "PRI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI] = "PLI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI] = "ARI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI] = "ALI";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS] = "PRS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS] = "PLS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS] = "ARS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS] = "ALS";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR] = "IPR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR] = "SPR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR] = "IAR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR] = "SAR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL] = "IPL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL] = "SPL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL] = "IAL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL] = "SAL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR] = "PIR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR] = "PSR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR] = "AIR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR] = "ASR";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL] = "PIL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL] = "PSL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL] = "AIL";
  m_CodeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL] = "ASL";
  return m_CodeToString[code];
  }


template <class TPixel, unsigned int VDim>
std::string
PrintImageInfo<TPixel, VDim>
::GetRAICodeFromDirectionMatrix(vnl_matrix_fixed<double, VDim, VDim> dir)
{
  // Compute the RAI code from the direction matrix
  char rai[VDim+1]; rai[VDim] = 0;
  bool oblique = false;

  char codes[3][2] = { {'R', 'L'}, {'A', 'P'}, {'I', 'S'}};

  for(size_t i = 0; i < VDim; i++)
    {
    // Get the direction cosine for voxel direction i
    vnl_vector_fixed<double, VDim> dcos = dir.get_column(i);
    double dabsmax = dcos.inf_norm();

    for(size_t j = 0; j < VDim; j++)
      {
      double dabs = fabs(dcos[j]);
      size_t dsgn = dcos[j] > 0 ? 0 : 1;

      if(dabs == 1.0)
        {
        rai[i] = codes[j][dsgn];
        }
      else if(dabs == dabsmax)
        {
        oblique = true;
        rai[i] = codes[j][dsgn];
        }
      }
    }
      
  if(oblique)
    {
    std::ostringstream sout;
    sout << "Oblique, closest to " << rai;
    return sout.str();
    }
  else return string(rai);
}

template <class TPixel, unsigned int VDim>
void
PrintImageInfo<TPixel, VDim>
::operator() (bool full)
{
  // Get the input image
  ImagePointer image = c->m_ImageStack.back();

  // Print the image number
  c->sout() << "Image #" << c->m_ImageStack.size() << ":";

  // Compute the bounding box
  RealVector bb0, bb1, ospm;
  for(size_t i = 0; i < VDim; i++)
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    ospm[i] = -image->GetOrigin()[i] / image->GetSpacing()[i];
    }

  // Compute the intensity range of the image
  size_t n = image->GetBufferedRegion().GetNumberOfPixels();
  double *vox = image->GetBufferPointer();
  double iMax = vox[0], iMin = vox[0], iMean = vox[0];
  for(size_t i = 1; i < n; i++)
    {
    iMax = (iMax > vox[i]) ? iMax : vox[i];
    iMin = (iMin < vox[i]) ? iMin : vox[i];
    iMean += vox[i];
    }
  iMean /= n;

  // Short or long?
  if(!full) 
    {
    c->sout() << " dim = " << image->GetBufferedRegion().GetSize() << "; ";
    c->sout() << " bb = {[" << bb0 << "], [" << bb1 << "]}; ";
    c->sout() << " vox = " << image->GetSpacing() << "; ";
    c->sout() << " range = [" << iMin << ", " << iMax << "]; ";
    c->sout() << " orient = " << GetRAICodeFromDirectionMatrix(image->GetDirection().GetVnlMatrix()) << endl;
    c->sout() << endl;
    }
  else
    {
    c->sout() << endl;
    c->sout() << "  Image Dimensions   : " << image->GetBufferedRegion().GetSize() << endl;
    c->sout() << "  Bounding Box       : " << "{[" << bb0 << "], [" << bb1 << "]}" << endl;
    c->sout() << "  Voxel Spacing      : " << image->GetSpacing() << endl;
    c->sout() << "  Intensity Range    : [" << iMin << ", " << iMax << "]" << endl;
    c->sout() << "  Mean Intensity     : " << iMean << endl;
    c->sout() << "  Canon. Orientation : " << GetRAICodeFromDirectionMatrix(image->GetDirection().GetVnlMatrix()) << endl;
    c->sout() << "  Direction Cos Mtx. : " << endl;
    c->PrintMatrix(c->sout(), image->GetDirection().GetVnlMatrix());

    // Print NIFTI s-form matrix (check against freesurfer's MRIinfo)
    c->sout() << "  Voxel->RAS x-form  : " << endl;
    c->PrintMatrix(c->sout(), 
      image->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix(), "%12.5f ", "    ");

    //
    // Print metadata
    c->sout() << "  Image Metadata: " << endl;
    itk::MetaDataDictionary &mdd = image->GetMetaDataDictionary();
    itk::MetaDataDictionary::ConstIterator itMeta;
    for(itMeta = mdd.Begin(); itMeta != mdd.End(); ++itMeta)
      {
      // Get the metadata as a generic object
      string key = itMeta->first, v_string;
      itk::SpatialOrientation::ValidCoordinateOrientationFlags v_oflags = 
        itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID;

      if(itk::ExposeMetaData<string>(mdd, key, v_string))
        {
        // For some weird reason, some of the strings returned by this method
        // contain '\0' characters. We will replace them by spaces
        std::ostringstream sout("");
        for(unsigned int i=0;i<v_string.length();i++)
          if(v_string[i] >= ' ') sout << v_string[i];
        v_string = sout.str();

        // Make sure the value has more than blanks
        if(v_string.find_first_not_of(" ") != v_string.npos)
          {
          c->sout() << "    " << key << " = " << v_string;

          bool v_is_integer = v_string.find_first_not_of("0123456789") == v_string.npos;
          // Augment fields/lines with extra infomation
          if (key == string("datatype"))  // itkNiftiImageIO field datatype
            {
            if (v_is_integer)
              {
              int v_int = atoi(v_string.c_str());
              // Convert to string using nifti1_io.c function nifti_datatype_to_string
              c->sout() << " (" << nifti_datatype_to_string(v_int) << ")";
              }
            }
          c->sout() << endl;
          }
        }
      else if(itk::ExposeMetaData(mdd, key, v_oflags))
        {
        c->sout() << "    " << key << " = " << get_rai_code(v_oflags) << endl;
        }
      else 
        {
        bool rc = false;
        if(!rc) rc |= try_print_metadata<double>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<float>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<int>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<unsigned int>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<long>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<unsigned long>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<short>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<unsigned short>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<char>(c->sout(), mdd, key, 0);
        if(!rc) rc |= try_print_metadata<unsigned char>(c->sout(), mdd, key, 0);

        if(!rc)
          {
          c->sout() << "    " << key << " of unsupported type " 
            << itMeta->second->GetMetaDataObjectTypeName() << endl;
          }
        }
      }

    }

}

// Invocations
template class PrintImageInfo<double, 2>;
template class PrintImageInfo<double, 3>;
template class PrintImageInfo<double, 4>;
