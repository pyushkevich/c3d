/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPovRayDF3ImageIOFactory.h,v $
  Language:  C++
  Date:      $Date: 2008/11/08 01:42:03 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPovRayDF3ImageIOFactory_h
#define __itkPovRayDF3ImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class PovRayDF3ImageIOFactory
 * \brief Create instances of PovRayDF3ImageIO objects using an object factory.
 */
class ITK_EXPORT PovRayDF3ImageIOFactory : public ObjectFactoryBase
{
public:  
  /** Standard class typedefs. */
  typedef PovRayDF3ImageIOFactory   Self;
  typedef ObjectFactoryBase  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Class methods used to interface with the registered factories. */
  virtual const char* GetITKSourceVersion(void) const ITK_OVERRIDE;
  virtual const char* GetDescription(void) const ITK_OVERRIDE;
  
  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PovRayDF3ImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    PovRayDF3ImageIOFactory::Pointer PovRayDF3Factory = PovRayDF3ImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(PovRayDF3Factory);
  }

protected:
  PovRayDF3ImageIOFactory();
  ~PovRayDF3ImageIOFactory();

private:
  PovRayDF3ImageIOFactory(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
  
  
} // end namespace itk

#endif
