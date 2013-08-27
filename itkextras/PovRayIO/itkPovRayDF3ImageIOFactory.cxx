/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPovRayDF3ImageIOFactory.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/08 01:42:03 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkPovRayDF3ImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkPovRayDF3ImageIO.h"
#include "itkVersion.h"

  
namespace itk
{

PovRayDF3ImageIOFactory::PovRayDF3ImageIOFactory()
{
  this->RegisterOverride("itkImageIOBase",
                         "itkPovRayDF3ImageIO",
                         "PovRay DF3 Image IO",
                         1,
                         CreateObjectFunction<PovRayDF3ImageIO>::New());
}
  
PovRayDF3ImageIOFactory::~PovRayDF3ImageIOFactory()
{
}

const char* 
PovRayDF3ImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char* 
PovRayDF3ImageIOFactory::GetDescription() const
{
  return "PovRay DF3 ImageIO Factory, allows the saving of PovRayDF3 images from Insight";
}

} // end namespace itk

