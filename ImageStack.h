/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ImageStack.h
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

#ifndef __ImageStack_h_
#define __ImageStack_h_

#include <itkSmartPointer.h>
#include <ConvertException.h>

/** A wrapper around vector<ImageType> that throws exceptions */
template <class ImageType>
class ImageStack
{
public:

  typedef itk::SmartPointer<ImageType> ImagePointer;
  typedef std::vector<ImagePointer> StackType;
  typedef typename StackType::iterator iterator;
  typedef typename StackType::const_iterator const_iterator;

  void push_back(ImageType *src)
  {
    m_Stack.push_back(src);
  }

  void pop_back()
  {
    if(m_Stack.size())
      m_Stack.pop_back();
    else throw StackAccessException();
  }

  ImageType* back() const
  {
    if(m_Stack.size())
      return m_Stack.back();
    else throw StackAccessException();
  }

  ImageType* front() const
  {
    if(m_Stack.size())
      return m_Stack.front();
    else throw StackAccessException();
  }

  ImageType * operator [] (size_t n) const
  {
    if(m_Stack.size() > n)
      return m_Stack[n];
    else throw StackAccessException();
  }

  size_t size() const 
  {
    return m_Stack.size();
  }

  void clear()
  {
    m_Stack.clear();
  }

  iterator end() { return m_Stack.end(); } 
  const_iterator end() const { return m_Stack.end(); }

  iterator begin() { return m_Stack.begin(); } 
  const_iterator begin() const { return m_Stack.begin(); }

  void insert(iterator &it, ImageType *image)
  {
    m_Stack.insert(it, image);
  }

protected:

  StackType m_Stack;
};


#endif
