#include "SplitMultilabelImage.h"
#include "ThresholdImage.h"
#include "UpdateMetadataKey.h"
#include <set>
#include "vnl/vnl_math.h"

using namespace std;

template <class TPixel, unsigned int VDim>
void
SplitMultilabelImage<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create a list of all the finite values in the image
  set<double> sval;
  for(ConstIterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    double val = it.Get();
    if(vnl_math_isfinite(val))
      sval.insert(val);
    }
  
  // The number of finite values should be reasonable
  if(sval.size() > 255)
    throw ConvertException("Number of labels passed on to -split exceeds 255");
  else if(sval.size() == 0)
    throw ConvertException("No finite labels passed on to -split");

  // A report
  *c->verbose << "Splitting #" << c->m_ImageStack.size() << " into " 
    << sval.size() << " binary images" << endl;

  // Pop the image
  c->m_ImageStack.pop_back();

  // Clear the label set
  c->GetSplitLabelSet();

  // Generate a bunch of binary copies (unfortunately, we store them as double)
  for(set<double>::iterator it = sval.begin(); it != sval.end(); ++it)
    {
    // Push our image back on the stack
    c->m_ImageStack.push_back(img);

    // Perform thresholding
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(*it, *it, 1.0, 0.0);

    // Add label to the label set
    c->GetSplitLabelSet().push_back(*it);
    }  
}

// Invocations
template class SplitMultilabelImage<double, 2>;
template class SplitMultilabelImage<double, 3>;
