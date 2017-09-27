/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPovRayDF3ImageIO.h,v $
  Language:  C++
  Date:      $Date: 2008/11/08 01:42:03 $
  Version:   $1.0$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPovRayDF3ImageIO_h
#define __itkPovRayDF3ImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <fstream>
#include <string>
#include <map>
#include "itkImageIOBase.h"
#include "itkSpatialOrientation.h"
#include <stdio.h>

namespace itk
{
  
/** \class PovRayDF3ImageIO
 *
 *  \brief Read PovRayDF3Image file format. 
 *
 *  \ingroup IOFilters
 *
 */
class ITK_EXPORT PovRayDF3ImageIO : public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef PovRayDF3ImageIO            Self;
  typedef ImageIOBase  Superclass;
  typedef SmartPointer<Self>  Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PovRayDF3ImageIO, Superclass);

  /*-------- This part of the interfaces deals with reading data. ----- */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char*) ITK_OVERRIDE { return false; }

  /** Set the spacing and dimension information for the set filename. */
  virtual void ReadImageInformation() ITK_OVERRIDE { }
  
  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void* buffer) ITK_OVERRIDE { }

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can write the
   * file specified. */
  virtual bool CanWriteFile(const char*) ITK_OVERRIDE;

  /** Set the spacing and dimension information for the set filename. */
  virtual void WriteImageInformation() ITK_OVERRIDE;
  
  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegions has been set properly. */
  virtual void Write(const void* buffer) ITK_OVERRIDE;

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

protected:
  PovRayDF3ImageIO() {}
  ~PovRayDF3ImageIO() {}
  
private:
  PovRayDF3ImageIO(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool CheckExtension(const char*);
};

} // end namespace itk

#endif // __itkPovRayDF3ImageIO_h
