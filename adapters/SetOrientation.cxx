#include "SetOrientation.h"

template <class TPixel, unsigned int VDim>
void
SetOrientation<TPixel, VDim>
::operator() (std::string rai)
{
  // Check the RAI code validity
  if(rai.length() != 3)
    throw ConvertException("Orientation code %s is not 3 characters long", rai.c_str());

  // Only valid for 3D images
  if(VDim != 3)
    throw ConvertException("Orientation codes only valid for 3D images");
  
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create a direction matrix
  vnl_matrix_fixed<double, 3, 3> eye, dm;
  eye.set_identity(); dm.set_identity();

  // RAI codes
  char codes[3][2] = { {'R', 'L'}, {'A', 'P'}, {'I', 'S'}};
  
  for(size_t i = 0; i < 3; i++)
    {
    bool matched = false;
    for(size_t j = 0; j < 3; j++)
      {
      for(size_t k = 0; k < 2; k++)
        {
        if(toupper(rai[i]) == codes[j][k])
          {
          // Set the row of the direction matrix
          dm.set_column(i, (k==0 ? 1.0 : -1.0) * eye.get_row(j));

          // Clear that code (so that we catch orientation codes like RRS)
          codes[j][0] = codes[j][1] = 'X';

          // We found a code for i
          matched = true;
          }
        }
      }

    if(!matched)
      throw ConvertException("Orientation code %s is invalid", rai.c_str());
    }

  // Explain what's being done
  *c->verbose << "Setting orientation of " << c->m_ImageStack.size() << " to " << rai << endl;

  // Set the direction in the image
  img->SetDirection(itk::Matrix<double,VDim,VDim>(dm));
}

// Invocations
template class SetOrientation<double, 2>;
template class SetOrientation<double, 3>;
