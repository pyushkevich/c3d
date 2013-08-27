#ifndef __MRFVote_h_
#define __MRFVote_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class MRFVote : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  MRFVote(Converter *c) : c(c) {}

  void operator() (double beta, size_t niter, bool flagUseSplitLabelSet);

private:
  Converter *c;

};

#endif

