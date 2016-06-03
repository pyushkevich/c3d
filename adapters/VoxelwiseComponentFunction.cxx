/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    VoxelwiseComponentFunction.cxx
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

#include "VoxelwiseComponentFunction.h"

#include "itkComposeImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"

namespace VoxelwiseComponentFunctionNamespace {

/** 
 * Functor for RGB to HSV conversion
 */

template <class TPixel, unsigned int VDim>
class RGB_to_HSV_Functor
{
public:
  typedef itk::FixedArray<TPixel, 3>   InputPixelType;
  typedef InputPixelType               OutputPixelType;

  OutputPixelType operator() (const InputPixelType &pix)
    {
    OutputPixelType out;

    double r = pix[0], g = pix[1], b = pix[2];
    TPixel &h = out[0], &s = out[1], &v = out[2];

    double c_min, c_max, delta;

    c_min = r < g ? r : g;
    c_min = c_min  < b ? c_min  : b;

    c_max = r > g ? r : g;
    c_max = c_max  > b ? c_max  : b;

    v = c_max;                                // v
    delta = c_max - c_min;

    if (delta < 0.00001)
      {
      s = 0;
      h = 0; // undefined, maybe nan?
      return out;
      }
    if( c_max > 0.0 ) 
      { 
      // NOTE: if Max is == 0, this divide would cause a crash
      s = (delta / c_max);                  // s
      } 
    else 
      {
      // if c_max is 0, then r = g = b = 0              
      // s = 0, v is undefined
      s = 0.0;
      h = NAN;                            // its now undefined
      return out;
      }
    if( r >= c_max )                           // > is bogus, just keeps compilor happy
      h = ( g - b ) / delta;        // between yellow & magenta
    else
      if( g >= c_max )
        h = 2.0 + ( b - r ) / delta;  // between cyan & yellow
    else
        h = 4.0 + ( r - g ) / delta;  // between magenta & cyan

    h *= 60.0;                              // degrees

    if( h < 0.0 )
        h += 360.0;

    return out;
    }
};


}; // namespace

using namespace VoxelwiseComponentFunctionNamespace;

// This function does the actual work
template <class TPixel, unsigned int VDim>
template <class TFunctor>
void 
VoxelwiseComponentFunction<TPixel, VDim>
::Apply(const TFunctor &functor, const char *func_name)
{
  // Number of components needed
  typedef typename TFunctor::InputPixelType InputPixel;
  typedef typename TFunctor::OutputPixelType OutputPixel;
  int n_comp = InputPixel::Dimension;
  int n_out_comp = OutputPixel::Dimension;
  int n_stack = c->m_ImageStack.size();

  // Get the images from the stack in forward order
  if(n_stack < n_comp)
    throw ConvertException("Too few components on the stack for VoxelwiseComponentFunction");

  // Create the filter
  typedef itk::Image<InputPixel, VDim> InputCompositeType;
  typedef itk::Image<OutputPixel, VDim> OutputCompositeType;
  typedef itk::ComposeImageFilter<ImageType, InputCompositeType> InputComposer;
  typename InputComposer::Pointer compose = InputComposer::New();

  // Take the last n_comp components
  for(int k = 0; k < n_comp; k++)
    compose->SetInput(k, c->m_ImageStack[n_stack + k - n_comp]);

  // Create the unary filter
  typedef itk::UnaryFunctorImageFilter<InputCompositeType, OutputCompositeType, TFunctor> Mapper;
  typename Mapper::Pointer mapper = Mapper::New();
  mapper->SetInput(compose->GetOutput());

  // Explain what we are doing
  *c->verbose 
    << "Applying voxelwise function " << func_name << " to #"
    << n_stack - n_comp << " - #" << n_stack - 1 << endl;

  // Update the mapper
  mapper->Update();

  // Pop the old components
  for (int i = 0; i < n_comp; i++)
    c->m_ImageStack.pop_back();

  // Push the new ones
  for (int i = 0; i < n_out_comp; i++)
    {
    typedef itk::NthElementImageAdaptor<OutputCompositeType, typename ImageType::PixelType> CompAdaptor;
    typename CompAdaptor::Pointer adaptor = CompAdaptor::New();
    adaptor->SelectNthElement(i);
    adaptor->SetImage(mapper->GetOutput());

    typedef itk::CastImageFilter<CompAdaptor, ImageType> Caster;
    typename Caster::Pointer caster = Caster::New();
    caster->SetInput(adaptor);
    caster->Update();

    c->m_ImageStack.push_back(caster->GetOutput());
    }
}


template <class TPixel, unsigned int VDim>
void
VoxelwiseComponentFunction<TPixel, VDim>
::operator() (const char *function)
{
  // Delegate to the right function
  if(!strcmp(function, "rgb2hsv"))
    {
    RGB_to_HSV_Functor<TPixel, VDim> functor;
    this->Apply(functor, function);
    }
  else
    throw ConvertException("Unknown function %s in -voxelwise-function", function);
}

// Invocations
template class VoxelwiseComponentFunction<double, 2>;
template class VoxelwiseComponentFunction<double, 3>;
template class VoxelwiseComponentFunction<double, 4>;
