#include "TileImages.h"

#include <itkTileImageFilter.h>

template <class TPixel, unsigned int VDim>
void
TileImages<TPixel, VDim>
::operator() (const std::string &tileParam)
{
  // Create the tile filter
  typedef typename itk::TileImageFilter<ImageType, ImageType> TileFilterType;
  typename TileFilterType::Pointer fltTile = TileFilterType::New();

  // Add all of the images as input to the tile filter
  for(int i = 0; i < c->m_ImageStack.size(); i++)
    {
    fltTile->SetInput(i, c->m_ImageStack[i]);
    }

  // Set the layout. The layout can either be a letter (x, y, z) which means tiling
  // in one dimension (i.e., tile a bunch of pngs) or it can be an integer vector
  // 2x2x0 which matches the input of the tile filter
  typename TileFilterType::LayoutArrayType loArray;
  
  if(tileParam == "x" || tileParam == "X")
    {
    loArray.Fill(1);
    loArray[0] = c->m_ImageStack.size();
    }
  else if(tileParam == "y" || tileParam == "Y")
    {
    loArray.Fill(1);
    loArray[1] = c->m_ImageStack.size();
    }
  else if(tileParam == "z" || tileParam == "Z")
    {
    if(VDim < 3) throw ConvertException("Can not tile in z-dimension using c2d, use c3d");
    loArray.Fill(1);
    loArray[2] = c->m_ImageStack.size();
    }
  else
    {
    SizeType sz = c->ReadSizeVector(tileParam.c_str());
    for(int i = 0; i < VDim; i++)
      loArray[i] = sz[i];
    }

  fltTile->SetLayout(loArray);

  *c->verbose << "Tiling " << c->m_ImageStack.size() 
    << " images using layout " << loArray << endl; 

  fltTile->Update();
  
  // Put result on stack
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(fltTile->GetOutput());
}

// Invocations
template class TileImages<double, 2>;
template class TileImages<double, 3>;
