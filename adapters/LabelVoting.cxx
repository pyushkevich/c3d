/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    LabelVoting.cxx
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

#include <sstream>

#include "LabelVoting.h"
#include "itkCastImageFilter.h"
#include "itkLabelVotingImageFilter.h"

template <class TPixel, unsigned int VDim>
void
LabelVoting<TPixel, VDim>
::operator() (std::vector<int> cmd_args)
{
  // Error message prefix
  std::string err_msg  = "LabelVoting adapter - ";

  // Check and parse command arguments
  if (cmd_args.size() == 0)
  {
    err_msg += "Require at least one integer argument";
    throw ConvertException(err_msg.c_str());
  }
  if (cmd_args[0]<0)
  {
    err_msg += "Require value of the undecided pixels to be unsigned: ";
    err_msg += std::to_string(cmd_args[0]);
    throw ConvertException(err_msg.c_str());
  }
  unsigned int label_for_undecided_pixels = cmd_args[0];
  unsigned int stack_size = c->m_ImageStack.size();
  unsigned int number_of_stack_images_to_use = stack_size;
  if (cmd_args.size() >= 2)
  {
    if (cmd_args[1]<=0)
    {
      err_msg += "Require the requested number of images to be larger then zero: ";
      err_msg += std::to_string(cmd_args[1]);
      throw ConvertException(err_msg.c_str());
    }
    number_of_stack_images_to_use = cmd_args[1];
  }
  if (cmd_args.size() > 2)
  {
    err_msg += "Takes maximum two integer argument: ";
    err_msg += std::to_string(cmd_args.size());
    throw ConvertException(err_msg.c_str());
  }

  // Check if expected input images are available on the stack
  if(stack_size <= 0)
  {
    err_msg += "Require one or more multilabel images on the stack: ";
    err_msg += std::to_string(stack_size);
    throw ConvertException(err_msg.c_str());
  }
  if (number_of_stack_images_to_use > stack_size)
  {
    err_msg += "Require the requested number of images: ";
    err_msg += std::to_string(number_of_stack_images_to_use);
    err_msg += " is less than or equal to the stack size: ";
    err_msg += std::to_string(stack_size);
    throw ConvertException(err_msg.c_str());
  }

  // Define labelmap type
  typedef itk::OrientedRASImage<unsigned int, VDim> LabelImageType;

  // Create LabelVotingImageFilter
	typedef itk::LabelVotingImageFilter<LabelImageType, LabelImageType> FilterType;
	typename FilterType::Pointer filter = FilterType::New();
  unsigned int start_index = stack_size - number_of_stack_images_to_use;
  for (unsigned int i=start_index; i<stack_size; i++)
  {
    // Convert input from c3ds stack format in doubles to defined labelmap type
    using CastFilterType = itk::CastImageFilter<ImageType, LabelImageType>;
    auto castfilter = CastFilterType::New();
    castfilter->SetInput(c->m_ImageStack[i]);
    castfilter->Update();
    // Add converted input image to the LabelVotingImageFilter
		filter->PushBackInput(castfilter->GetOutput());
  }
  // SetLabelForUndecidedPixels on LabelVotingImageFilter
  filter->SetLabelForUndecidedPixels(label_for_undecided_pixels);

  // Convert output from the defined labelmap type and back to c3ds native double format
  using ReverseCastFilterType = itk::CastImageFilter<LabelImageType, ImageType>;
  auto reverse_castfilter = ReverseCastFilterType::New();
  reverse_castfilter->SetInput(filter->GetOutput());
  reverse_castfilter->Update();
  ImagePointer result = reverse_castfilter->GetOutput();

  // Write to verbose output
  *c->verbose << "LabelVoting (" << label_for_undecided_pixels;
  *c->verbose << "," << number_of_stack_images_to_use << "):";
  for (unsigned int i=start_index; i<stack_size; i++)
    *c->verbose << " #" << i;
  *c->verbose << endl;

  // Put result on stack and pop the images used
  for (unsigned int i=start_index; i<stack_size; i++)
    c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class LabelVoting<double, 2>;
template class LabelVoting<double, 3>;
template class LabelVoting<double, 4>;
