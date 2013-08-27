#include "LandmarksToSpheres.h"

template <class TPixel, unsigned int VDim>
void
LandmarksToSpheres<TPixel, VDim>
::operator() (const char *fnland, double radius)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Open the landmarks file
  FILE *f = fopen(fnland, "rt");
  if(!f)
    throw ConvertException("Unable to read file %s", fnland);

  // Define a landmark
  typedef itk::Point<double, VDim> PointType;
  typedef std::pair<PointType, double> Landmark;
  std::list<Landmark> lms;

  // Read each landmark in turn
  while(!feof(f))
    {
    PointType x; double label = 0; int rc;
    if(VDim == 2)
      rc = fscanf(f, "%lf %lf %lf\n", &x[0], &x[1], &label);
    else
      rc = fscanf(f, "%lf %lf %lf %lf\n", &x[0], &x[1], &x[2], &label);
    if (rc != VDim + 1)
      throw ConvertException("Error reading line %d in file %s", lms.size(), fnland);
    lms.push_back(make_pair(x, label));
    }
  fclose(f);

  // How many landmarks?
  if(lms.size() == 0)
    throw ConvertException("No landmarks were read from %s", fnland);
  *c->verbose << "Placing " << lms.size() << " landmarks in #" 
    << c->m_ImageStack.size() << std::endl;

  // Now we basically take the SNAP code for filling landmarks
  for(typename std::list<Landmark>::iterator it = lms.begin();
    it != lms.end(); it++)
    {
    // The landmark
    PointType x = it->first;
    double label = it->second;

    // Convert the landmark center to a continuous index
    itk::ContinuousIndex<double, VDim> idxCenter;
    img->TransformRASPhysicalPointToContinuousIndex(x, idxCenter);

    // We'll track the bounds for this bubble
    itk::ContinuousIndex<double, VDim> idxUpper = idxCenter, idxLower = idxCenter;

    // Iterate over all permutations of (-1,1) for x, y and z
    unsigned int nperm = 1 << VDim;
    for(unsigned int iperm = 0; iperm < nperm; iperm++)
      {
      // Create a point at the vertex of the cube/square with radius r
      PointType xTest = x;
      for(unsigned int i = 0; i < VDim; i++)
        {
        int sign = ((1 << i) & iperm) > 0 ? 1 : -1;
        xTest[i] += sign * radius;
        }

      // Map point to index space
      itk::ContinuousIndex<double, VDim> cTest;
      img->TransformRASPhysicalPointToContinuousIndex(xTest, cTest);

      // Update the bounding box
      for(unsigned int i = 0; i < VDim; i++)
        {
        idxLower[i] = (idxLower[i] > cTest[i]) ? cTest[i] : idxLower[i];
        idxUpper[i] = (idxUpper[i] < cTest[i]) ? cTest[i] : idxUpper[i];
        }
      }

    // Create a region in which to create the bubble
    SizeType szBubble; IndexType idxBubble;
    for(unsigned int i = 0; i < VDim; i++)
      {
      idxBubble[i] = floor(idxLower[i]);
      szBubble[i] = ceil(idxUpper[i]) - floor(idxLower[i]);
      }
    RegionType regBubble(idxBubble, szBubble);

    // Crop by the image region and skip if fully outside
    bool cropped = regBubble.Crop(img->GetBufferedRegion());
    if(!cropped)
      {
      *c->verbose << "  Landmark " << x << " fully outside image region" << std::endl;
      continue;
      }

    // Iterate over this region
    double rsq = radius * radius;
    for(Iterator itb(img, regBubble); !itb.IsAtEnd(); ++itb)
      {
      // Get the physical point for this index
      PointType pit;
      img->TransformIndexToRASPhysicalPoint(itb.GetIndex(), pit);

      // Check if it's in the bubble
      if(pit.SquaredEuclideanDistanceTo(x) <= rsq)
        itb.Set(label);
      }
    }
}

// Invocations
template class LandmarksToSpheres<double, 2>;
template class LandmarksToSpheres<double, 3>;
