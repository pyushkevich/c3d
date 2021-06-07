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
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkUnaryFunctorImageFilter.h"
#include "UnaryFunctorVectorImageFilter.h"

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

template <class TPixel, unsigned int VDim>
class SoftmaxFunctor
{
public:
  typedef itk::VariableLengthVector<TPixel>   InputPixelType;
  typedef InputPixelType                      OutputPixelType;

  static unsigned int GetNumberOfComponentsPerPixel() { return 0; }

  OutputPixelType operator() (const InputPixelType &pix)
    {
    unsigned int n = pix.GetSize();
    OutputPixelType out(n);
    double sum = 0;
    for(unsigned int i = 0; i < n; i++)
      {
      out[i] = exp(pix[i]);
      sum += out[i];
      }
    for(unsigned int i = 0; i < n; i++)
      out[i] /= sum;

    return out;
    }
};

template <class TInputPixel, class TOutputPixel, unsigned int VDim, class TFunctor> class MappingSpecialization
{
public:
  typedef itk::Image<TInputPixel, VDim> InputImageType;
  typedef itk::Image<TOutputPixel, VDim> OutputImageType;

  static unsigned int GetNumberOfInputComponents() { return TInputPixel::Dimension; }
  static unsigned int GetNumberOfOutputComponents() { return TOutputPixel::Dimension; }

  static void Map(InputImageType *input, OutputImageType *output)
    {
    // Create the unary filter
    typedef itk::UnaryFunctorImageFilter<InputImageType, OutputImageType, TFunctor> Mapper;
    typename Mapper::Pointer mapper = Mapper::New();
    mapper->SetInput(input);
    mapper->GraftOutput(output);
    mapper->Update();
    }
};

template <class TPixel, unsigned int VDim, class TFunctor>
class MappingSpecialization<
  itk::VariableLengthVector<TPixel>, itk::VariableLengthVector<TPixel>, VDim, TFunctor>
{
public:
  typedef itk::VectorImage<TPixel, VDim> InputImageType;
  typedef itk::VectorImage<TPixel, VDim> OutputImageType;

  static unsigned int GetNumberOfInputComponents() { return 0; }
  static unsigned int GetNumberOfOutputComponents() { return 0; }

  static void Map(InputImageType *input, OutputImageType *output)
    {
    // Allocate the output image
    output->CopyInformation(input);
    output->SetRegions(input->GetBufferedRegion());
    output->Allocate();

    // Create the unary filter
    typedef itk::UnaryFunctorImageFilter<InputImageType, OutputImageType, TFunctor> Mapper;
    typename Mapper::Pointer mapper = Mapper::New();
    mapper->SetInput(input);
    mapper->GraftOutput(output);
    mapper->Update();
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
  // Create a specialization based on the functor
  typedef typename TFunctor::InputPixelType InputPixel;
  typedef typename TFunctor::OutputPixelType OutputPixel;

  // Create specialization object
  typedef MappingSpecialization<InputPixel, OutputPixel, VDim, TFunctor> Specialization;

  // Number of components needed
  unsigned int n_comp = Specialization::GetNumberOfInputComponents();
  unsigned int n_out_comp = Specialization::GetNumberOfOutputComponents();
  unsigned int n_stack = c->m_ImageStack.size();

  // Zero components means use all
  if(n_comp == 0) n_comp = n_stack;
  if(n_out_comp == 0) n_out_comp = n_stack;

  // Get the images from the stack in forward order
  if(n_stack < n_comp)
    throw ConvertException("Too few components on the stack for VoxelwiseComponentFunction");

  // Create the filter
  typedef typename Specialization::InputImageType InputCompositeType;
  typedef typename Specialization::OutputImageType OutputCompositeType;
  typedef itk::ComposeImageFilter<ImageType, InputCompositeType> InputComposer;
  typename InputComposer::Pointer compose = InputComposer::New();

  // Take the last n_comp components
  for(int k = 0; k < n_comp; k++)
    compose->SetInput(k, c->m_ImageStack[n_stack + k - n_comp]);
  compose->Update();

  // Explain what we are doing
  *c->verbose 
    << "Applying voxelwise function " << func_name << " to #"
    << n_stack - n_comp << " - #" << n_stack - 1 << endl;

  // Create the unary filter
  typename OutputCompositeType::Pointer result = OutputCompositeType::New();
  Specialization::Map(compose->GetOutput(), result);

  // Pop the old components
  for (int i = 0; i < n_comp; i++)
    c->m_ImageStack.pop_back();

  // Push the new ones
  for (int i = 0; i < n_out_comp; i++)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<OutputCompositeType, ImageType> CompFilter;
    typename CompFilter::Pointer caster = CompFilter::New();
    caster->SetInput(result);
    caster->SetIndex(i);
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
  else if(!strcmp(function, "softmax"))
    {
    SoftmaxFunctor<TPixel, VDim> functor;
    this->Apply(functor, function);
    }
  else
    throw ConvertException("Unknown function %s in -voxelwise-function", function);
}

// Invocations
template class VoxelwiseComponentFunction<double, 2>;
template class VoxelwiseComponentFunction<double, 3>;
template class VoxelwiseComponentFunction<double, 4>;
