/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ScalarToRGB.cxx
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

#include "ScalarToRGB.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include <map>

template <class TPixel, unsigned int VDim>
void
ScalarToRGB<TPixel, VDim>
::operator() (const std::string &colormap, double int_min, double int_max)
{
  // Check the input parameters
  typedef itk::RGBPixel<unsigned char> RGBPixel;
  typedef itk::Image<RGBPixel, VDim> RGBImageType;
  typedef itk::ScalarToRGBColormapImageFilter<ImageType, RGBImageType> RGBFilterType; 
#if (ITK_VERSION_MAJOR == 5 && ITK_VERSION_MINOR >= 1) || ITK_VERSION_MAJOR > 5
  typedef typename RGBFilterType::RGBColormapFilterEnum ColormapEnumType;
#else
  typedef typename RGBFilterType::ColormapEnumType ColormapEnumType;
#endif
  typedef std::map<std::string, ColormapEnumType> ColormapNameMap;
  ColormapNameMap clmap;

  // Set up the map
  clmap["red"]=ColormapEnumType::Red;
  clmap["green"]=ColormapEnumType::Green;
  clmap["blue"]=ColormapEnumType::Blue;
  clmap["grey"]=ColormapEnumType::Grey;
  clmap["hot"]=ColormapEnumType::Hot;
  clmap["cool"]=ColormapEnumType::Cool;
  clmap["spring"]=ColormapEnumType::Spring;
  clmap["summer"]=ColormapEnumType::Summer;
  clmap["autumn"]=ColormapEnumType::Autumn;
  clmap["winter"]=ColormapEnumType::Winter;
  clmap["copper"]=ColormapEnumType::Copper;
  clmap["jet"]=ColormapEnumType::Jet;
  clmap["hsv"]=ColormapEnumType::HSV;
  clmap["overunder"]=ColormapEnumType::OverUnder;

  // Look it up
  typename ColormapNameMap::iterator it = clmap.find(colormap);
  if(it == clmap.end())
    throw ConvertException("Unknown colormap %s", colormap.c_str());

  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create the filter
  typename RGBFilterType::Pointer filter = RGBFilterType::New();
  filter->SetInput(img);
  filter->SetColormap(it->second);

  // If range is not specified, use the colormap as is
  if(int_min != 0 || int_max != 0)
    {
    typename RGBFilterType::ColormapType::Pointer colormap = filter->GetModifiableColormap();
    colormap->SetMinimumInputValue(int_min);
    colormap->SetMaximumInputValue(int_max);
    filter->SetUseInputImageExtremaForScaling(false);
    }

  // Verbose
  *c->verbose << "Mapping #" << c->m_ImageStack.size() 
    << " to RGB using color map " << colormap << endl;

  // Do it!
  filter->Update();

  // Split into RGB components. 
  c->m_ImageStack.pop_back();
  for(int i = 0; i < 3; i++)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<RGBImageType, ImageType> CompFilterType;
    typename CompFilterType::Pointer fltComp = CompFilterType::New();
    fltComp->SetInput(filter->GetOutput());
    fltComp->SetIndex(i);
    fltComp->Update();
    c->m_ImageStack.push_back(fltComp->GetOutput());
    }
}

// Invocations
template class ScalarToRGB<double, 2>;
template class ScalarToRGB<double, 3>;
template class ScalarToRGB<double, 4>;
