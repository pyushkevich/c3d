/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SwapDimensions.cxx
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

#include "SwapDimensions.h"
#include "itkSpatialOrientation.h"
#include "itkOrientImageFilter.h"

using namespace itk::SpatialOrientation;

ValidCoordinateOrientationFlags GetOrientationFlagFromString(const std::string &code)
{
  // Convert the code to upper case
  string code_upper = code;
  for(int i = 0; i < code_upper.size(); i++)
    code_upper[i] = std::toupper(code_upper[i]);

  // Lazy :)
  std::map<std::string, ValidCoordinateOrientationFlags> emap;
  emap["RIP"] = ITK_COORDINATE_ORIENTATION_RIP;
  emap["LIP"] = ITK_COORDINATE_ORIENTATION_LIP;
  emap["RSP"] = ITK_COORDINATE_ORIENTATION_RSP;
  emap["LSP"] = ITK_COORDINATE_ORIENTATION_LSP;
  emap["RIA"] = ITK_COORDINATE_ORIENTATION_RIA;
  emap["LIA"] = ITK_COORDINATE_ORIENTATION_LIA;
  emap["RSA"] = ITK_COORDINATE_ORIENTATION_RSA;
  emap["LSA"] = ITK_COORDINATE_ORIENTATION_LSA;
  emap["IRP"] = ITK_COORDINATE_ORIENTATION_IRP;
  emap["ILP"] = ITK_COORDINATE_ORIENTATION_ILP;
  emap["SRP"] = ITK_COORDINATE_ORIENTATION_SRP;
  emap["SLP"] = ITK_COORDINATE_ORIENTATION_SLP;
  emap["IRA"] = ITK_COORDINATE_ORIENTATION_IRA;
  emap["ILA"] = ITK_COORDINATE_ORIENTATION_ILA;
  emap["SRA"] = ITK_COORDINATE_ORIENTATION_SRA;
  emap["SLA"] = ITK_COORDINATE_ORIENTATION_SLA;
  emap["RPI"] = ITK_COORDINATE_ORIENTATION_RPI;
  emap["LPI"] = ITK_COORDINATE_ORIENTATION_LPI;
  emap["RAI"] = ITK_COORDINATE_ORIENTATION_RAI;
  emap["LAI"] = ITK_COORDINATE_ORIENTATION_LAI;
  emap["RPS"] = ITK_COORDINATE_ORIENTATION_RPS;
  emap["LPS"] = ITK_COORDINATE_ORIENTATION_LPS;
  emap["RAS"] = ITK_COORDINATE_ORIENTATION_RAS;
  emap["LAS"] = ITK_COORDINATE_ORIENTATION_LAS;
  emap["PRI"] = ITK_COORDINATE_ORIENTATION_PRI;
  emap["PLI"] = ITK_COORDINATE_ORIENTATION_PLI;
  emap["ARI"] = ITK_COORDINATE_ORIENTATION_ARI;
  emap["ALI"] = ITK_COORDINATE_ORIENTATION_ALI;
  emap["PRS"] = ITK_COORDINATE_ORIENTATION_PRS;
  emap["PLS"] = ITK_COORDINATE_ORIENTATION_PLS;
  emap["ARS"] = ITK_COORDINATE_ORIENTATION_ARS;
  emap["ALS"] = ITK_COORDINATE_ORIENTATION_ALS;
  emap["IPR"] = ITK_COORDINATE_ORIENTATION_IPR;
  emap["SPR"] = ITK_COORDINATE_ORIENTATION_SPR;
  emap["IAR"] = ITK_COORDINATE_ORIENTATION_IAR;
  emap["SAR"] = ITK_COORDINATE_ORIENTATION_SAR;
  emap["IPL"] = ITK_COORDINATE_ORIENTATION_IPL;
  emap["SPL"] = ITK_COORDINATE_ORIENTATION_SPL;
  emap["IAL"] = ITK_COORDINATE_ORIENTATION_IAL;
  emap["SAL"] = ITK_COORDINATE_ORIENTATION_SAL;
  emap["PIR"] = ITK_COORDINATE_ORIENTATION_PIR;
  emap["PSR"] = ITK_COORDINATE_ORIENTATION_PSR;
  emap["AIR"] = ITK_COORDINATE_ORIENTATION_AIR;
  emap["ASR"] = ITK_COORDINATE_ORIENTATION_ASR;
  emap["PIL"] = ITK_COORDINATE_ORIENTATION_PIL;
  emap["PSL"] = ITK_COORDINATE_ORIENTATION_PSL;
  emap["AIL"] = ITK_COORDINATE_ORIENTATION_AIL;
  emap["ASL"] = ITK_COORDINATE_ORIENTATION_ASL;

  if(emap.find(code) != emap.end())
    return emap[code];

  else
    return ITK_COORDINATE_ORIENTATION_INVALID;
}

template <class TPixel, unsigned int VDim>
class SwapDimensions_OrientWorker
{
public:
  typedef typename ConvertAdapter<TPixel, VDim>::Converter ConverterType;
  void operator() (ConverterType *c, const std::string &code)
    {
    throw ConvertException("Swapping Dimensions in non-3D images is not yet implemented");
    }
};

template <class TPixel>
class SwapDimensions_OrientWorker<TPixel, 3>
{
public:
  typedef typename ConvertAdapter<TPixel, 3>::Converter ConverterType;
  typedef typename ConvertAdapter<TPixel, 3>::ImageType ImageType;
  void operator() (ConverterType *c, const std::string &code)
    {
    // Convert the orientation code to an ITK 
    ValidCoordinateOrientationFlags oflag = GetOrientationFlagFromString(code);
    if(oflag == itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID)
      throw ConvertException("Orientation flag %s is not valid", code.c_str());

    // Create the filter
    typedef itk::OrientImageFilter<ImageType, ImageType> FilterType;
    typename FilterType::Pointer fltOrient = FilterType::New();
    fltOrient->SetInput(c->PopImage());
    fltOrient->UseImageDirectionOn();
    fltOrient->SetDesiredCoordinateOrientation(oflag);

    // Report something
    c->PrintF("Swapping dimensions of #%d to achieve orientation %s\n", 
              c->GetStackSize(), code.c_str()); 

    fltOrient->Update();
    c->PushImage(fltOrient->GetOutput());
    }
};


template <class TPixel, unsigned int VDim>
void
SwapDimensions<TPixel, VDim>
::operator() (const std::vector<std::string> &code)
{
  // If the code is a single character string and VDim==3, we use OrientImageFilter
  if(code.size() == 1)
    {
    // Delegate to the specialized worker - this is because the underlying ITK filter 
    // only supports 3D images
    SwapDimensions_OrientWorker<TPixel, VDim> worker;
    worker(c, code.front());
    }
  else
    {
    throw ConvertException("Swapping Dimensions with arbitrary code is not yet implemented");
    }
}

// Invocations
template class SwapDimensions<double, 2>;
template class SwapDimensions<double, 3>;
template class SwapDimensions<double, 4>;
