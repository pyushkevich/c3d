#ifndef RANDOMFORESTCLASSIFYIMAGEFILTER_H
#define RANDOMFORESTCLASSIFYIMAGEFILTER_H

#include "itkImageToImageFilter.h"

template <class TPixel, class TLabel, int VDim> class RandomForestClassifier;

/**
 * @brief A class that takes multiple multi-component images and uses a
 * Gaussian mixture model to combine them into a single probability map.
 *
 * // TODO: derive this and the GMM filter from a common base class that
 * // simplifies working with multiple vector/scalar images
 */
template <class TInputImage, class TInputVectorImage, class TOutputImage, class TLabel>
class RandomForestClassifyImageFilter :
    public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename InputImageType::PixelType             InputPixelType;
  typedef typename InputImageType::InternalPixelType InputComponentType;
  typedef typename InputImageType::RegionType      InputImageRegionType;

  /** Define the corresponding vector image */
  typedef TInputVectorImage                        InputVectorImageType;

  /** Pixel Type of the output image */
  typedef TOutputImage                                  OutputImageType;
  typedef typename OutputImageType::PixelType           OutputPixelType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  typedef typename OutputImageType::Pointer          OutputImagePointer;

  /** Standard class typedefs. */
  typedef RandomForestClassifyImageFilter                          Self;
  typedef itk::ImageSource<OutputImageType>                  Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Classifier typedef */
  typedef TLabel                                              LabelType;
  typedef RandomForestClassifier<InputPixelType, LabelType,
                                 ImageDimension>         ClassifierType;

  /** Add a scalar input image */
  void AddScalarImage(InputImageType *image);

  /** Add a vector (multi-component) input image */
  void AddVectorImage(InputVectorImageType *image);

  /** Set the mixture model */
  void SetClassifier(ClassifierType *classifier);

  /**
   * Turn this flag on to generate class probabilities as outputs, as
   * opposed to the weighted sum of class probabilities, which is the
   * default output.
   */
  void SetGenerateClassProbabilities(bool flag);

  itkGetMacro(GenerateClassProbabilities, bool)

  /** Get the current classifier */
  itkGetMacro(Classifier, ClassifierType *);

  /** We need to override this method because of multiple input types */
  void GenerateInputRequestedRegion() ITK_OVERRIDE;

protected:

  RandomForestClassifyImageFilter();
  virtual ~RandomForestClassifyImageFilter();

  void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

  void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread,
                            itk::ThreadIdType threadId) ITK_OVERRIDE;

  void UpdateOutputs();

  ClassifierType *m_Classifier;

  bool m_GenerateClassProbabilities;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "RandomForestClassifyImageFilter.txx"
#endif



#endif // RANDOMFORESTCLASSIFYIMAGEFILTER_H
