/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPovRayDF3ImageIO.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/08 01:42:03 $
  Version:   $Revision: 1.1 $  

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkPovRayDF3ImageIO.h"
#include "itkIOCommon.h"
#include "itkExceptionObject.h"
#include "itkMetaDataObject.h"
#include "itkByteSwapper.h"
#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <math.h>

using namespace std;

namespace itk {

bool PovRayDF3ImageIO::CanWriteFile( const char * name )
{
  return CheckExtension(name);
}

void 
PovRayDF3ImageIO
::WriteImageInformation(void)
{
  // Check that the number of dimensions is correct
  if(GetNumberOfDimensions() != 3)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Unsupported number of dimensions");
    throw exception;
    }

  // Check that the component type is legal
  if(m_ComponentType != UINT && m_ComponentType != UCHAR && m_ComponentType != USHORT)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Unsupported pixel component type (uchar, ushort, uint are supported)");
    throw exception;
    }

  // Write the 3 shorts with image dimensions
  unsigned short xyz[3];
  xyz[0] = (unsigned short) m_Dimensions[0];  
  xyz[1] = (unsigned short) m_Dimensions[1];
  xyz[2] = (unsigned short) m_Dimensions[2];

  // Open a file handle
  ofstream fout(this->GetFileName(), ios_base::out | ios_base::binary);

  // Write the header
  itk::ByteSwapper<unsigned short>::SwapWriteRangeFromSystemToBigEndian(
    xyz, 3, &fout);

  // Write the data to file
  fout.close();
}

/** The write function is not implemented */
void 
PovRayDF3ImageIO
::Write( const void* buffer) 
{
  // Write the image information
  WriteImageInformation();
  
  // Open a file handle
  ofstream fout(this->GetFileName(), ios_base::app | ios_base::binary);

  // Write with a byte swap
  if(m_ComponentType == USHORT)
    itk::ByteSwapper<unsigned short>::SwapWriteRangeFromSystemToBigEndian(
      (unsigned short *) buffer, GetImageSizeInBytes() / sizeof(unsigned short), &fout);
  else if(m_ComponentType == UINT)
    itk::ByteSwapper<unsigned int>::SwapWriteRangeFromSystemToBigEndian(
      (unsigned int *) buffer, GetImageSizeInBytes() / sizeof(unsigned int), &fout);
  else if(m_ComponentType == UCHAR)
    itk::ByteSwapper<unsigned char>::SwapWriteRangeFromSystemToBigEndian(
      (unsigned char *) buffer, GetImageSizeInBytes() / sizeof(unsigned char), &fout);

  // Close the stream
  fout.close();
}

/** Print Self Method */
void PovRayDF3ImageIO::PrintSelf(ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
}


bool PovRayDF3ImageIO::CheckExtension(const char* filename)
{
  string fname = filename;
  if ( fname == "" )
  {
    itkDebugMacro(<< "No filename specified.");
    return false;
  }

  string::size_type giplPos = fname.rfind(".df3");
  return ((giplPos != string::npos) && (giplPos == fname.length() - 4));
}

} // end namespace itk
