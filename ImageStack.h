#ifndef __ImageStack_h_
#define __ImageStack_H_

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
