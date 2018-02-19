/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertImageND.h
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

#ifndef __ConvertImageND_h_
#define __ConvertImageND_h_

#include "itkOrientedRASImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"
#include "ImageStack.h"
#include "ConvertException.h"

#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

using namespace std;

template <class TPixel, unsigned int VDim> class ConvertAdapter;

template <class TPixel, unsigned int VDim> struct ConvertAlgorithmParameters;
class Documentation;

template<class TPixel, unsigned int VDim>
class ImageConverter
{
public:

  // Image typedef
  typedef itk::OrientedRASImage<TPixel, VDim> ImageType;
  typedef itk::Image<TPixel, VDim> UnorientedImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;
  typedef std::map<TPixel, vnl_vector_fixed<double, 4> > LabelToRGBAMap;

  // Complex stuff
  typedef std::complex<TPixel> ComplexPixel;
  typedef itk::OrientedRASImage<ComplexPixel, VDim> ComplexImageType;

  // Iterators
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIterator;

  // Interpolator
  typedef itk::InterpolateImageFunction<ImageType, double> Interpolator;
  
  ImageConverter();
  ~ImageConverter();
  int ProcessCommandLine(int argc, char *argv[]);

  friend class ConvertAdapter<TPixel, VDim>;

  // Copy image on stack
  void CopyImage();

  // Get bounding box of an image
  void GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1);

  // Check if N images on the stack have the same dimensions as the top image
  // (pass zero to check all images on the stack)
  bool CheckStackSameDimensions(size_t n);

  // Read label to RGBA mapping from file (SNAP format)
  static LabelToRGBAMap ReadLabelToRGBAMap(const char *fname);

  // Print a matrix in a nice way
  void PrintMatrix(
    std::ostream &sout, vnl_matrix<double> mat, 
    const char *fmt = "%10.4f ", const char *prefix = "");

  // Set label set for split/merge operations
  typedef std::vector<double> LabelSet;
  LabelSet &GetSplitLabelSet() 
    { return m_SplitLabelSet; }

  Interpolator *GetInterpolator() const
    { return m_Interpolator; }

  void SetInterpolator(Interpolator *interp)
    { m_Interpolator = interp; }

  // Read vectors, etc from command line
  SizeType ReadSizeVector(const char *vec);
  IndexType ReadIndexVector(const char *vec);
  RealVector ReadRealVector(const char *vec, bool is_point);
  RealVector ReadRealSize(const char *vec);
  TPixel ReadIntensityValue(const char *vec);

  void PrintManual(std::ostream &out);
  void PrintCommandListing(std::ostream &out);
  void PrintCommandHelp(std::ostream &out, const char *command);

  // Add variable
  void SetVariable(std::string name, ImagePointer image);

  // Get variable
  ImageType *GetVariable(std::string name);

private:

  // Internal functions
  int ProcessCommand(int argc, char *argv[]);


  // Templated write function
  template<class TOutPixel> void TemplatedWriteImage(const char *file, double xRoundFactor);

  // Map of variable names
  typedef map<string, ImagePointer> ImageVariableMap;
  ImageVariableMap m_ImageVars;

  // Image interpolator for all interpolation commands
  itk::SmartPointer<Interpolator> m_Interpolator;

  // Implementation of the 'foreach' loop
  size_t ForEachLoop(int argc, char *argv[]);

  // Implementation of the 'foreach-comp' loop
  size_t ForEachComponentLoop(int ncomp, int argc, char *argv[]);

  // Implementation of the 'accum' loop
  size_t AccumulateLoop(int argc, char *argv[]);

  // Write multiple images (for -oo and --oomc commands)
  int WriteMultiple(int argc, char *argv[], int n_comp, const char *command);

  // Type of loop we are in currently
  enum LoopType { LOOP_NONE = 0, LOOP_FOREACH, LOOP_ACCUM };
  LoopType m_LoopType;

  // Label set for split/merge
  LabelSet m_SplitLabelSet;
  
public:

  // Stack of images from the command line
  ImageStack<ImageType> m_ImageStack;

  // Get the last image from the stack and pop it off
  ImagePointer PopImage();

  // Get the last K images from the stack and pop them off
  std::vector<ImagePointer> PopNImages(unsigned int n);

  // Push a new image to the stack
  void PushImage(ImageType *image);

  // Get the size of the stack
  int GetStackSize();

  // Print something to verbose output, but using printf syntax
  void PrintF(const char *fmt, ...);

  // Get Nth image on the stack without removing it
  ImageType *PeekImage(int k);

  // Get Nth image on the stack without removing it
  ImageType *PeekLastImage();

  // Replace an image at the end of the stack with its copy
  ImageType *PopAndPushCopy();

  // Typeid of the image to be saved
  string m_TypeId;

  // Interpolation type
  string m_Interpolation;

  // Background intensity
  double m_Background;

  // Rounding factor (not used for float and double) is 0.5 unless -noround was specified
  double m_RoundFactor;

  // Whether SPM extensions are used
  bool m_FlagSPM;

  // Whether compression is used by default
  bool m_UseCompression;

  // Whether multicomponent images are split on read
  bool m_MultiComponentSplit;

  // Number of iterations for various algorithms
  size_t m_Iterations;

  // Parameters for various algorithms
  typedef ConvertAlgorithmParameters<TPixel, VDim> ParameterType;
  ParameterType *m_Param;

  // How % is handled for intensity specs
  enum PercentIntensityMode { PIM_QUANTILE, PIM_FGQUANTILE, PIM_RANGE };
  PercentIntensityMode m_PercentIntensityMode;

  // Verbose output stream
  std::ostringstream devnull;
  std::ostream *verbose;

  // Documentation object
  Documentation *m_Documentation;
};

#endif

