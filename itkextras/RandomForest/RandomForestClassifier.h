#ifndef RANDOMFORESTCLASSIFIER_H
#define RANDOMFORESTCLASSIFIER_H

#include <itkDataObject.h>
#include <itkObjectFactory.h>
#include <itkSize.h>
#include "Library/classification.h"
#include <map>

template <class dataT, class labelT> class Histogram;
template <class dataT, class labelT> class AxisAlignedClassifier;
template <class HistT, class ClassT, class dataT> class DecisionForest;

/**
 * This class encapsulates a Random Forest classifier
 */
template <class TPixel, class TLabel, int VDim>
class RandomForestClassifier : public itk::DataObject
{
public:

  // Standard ITK stuff
  typedef RandomForestClassifier Self;
  typedef itk::DataObject Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self)
  
  // typedefs
  typedef TPixel GreyType;
  typedef TLabel LabelType;
  typedef Histogram<GreyType, LabelType> RFHistogramType;
  typedef AxisAlignedClassifier<GreyType, LabelType> RFAxisClassifierType;
  typedef DecisionForest<RFHistogramType, RFAxisClassifierType, GreyType> RandomForestType;
  typedef std::map<size_t, LabelType> MappingType;
  typedef itk::Size<VDim> SizeType;

  // A list of weights for each class - used to construct speed image
  typedef std::vector<double> WeightArray;

  // Reset the classifier
  void Reset()
    {
    if(m_Forest)
      delete m_Forest;

    m_Forest = new RandomForestType(true);
    m_ClassToLabelMapping.clear();
    m_BiasParameter = 0.5;
    m_PatchRadius.Fill(0);
    m_UseCoordinateFeatures = false;
    m_ClassWeights.clear();
    }

  // Get the random forest
  itkGetMacro(Forest, RandomForestType *)

  // Get the patch radius
  itkGetMacro(PatchRadius, const SizeType &)
  itkSetMacro(PatchRadius, SizeType)

  /** Whether coordinates of the voxels are used as features */
  itkGetMacro(UseCoordinateFeatures, bool)
  itkSetMacro(UseCoordinateFeatures, bool)

  // Set the bias parameter (adjusts the mapping of FG probability to speed)
  itkGetMacro(BiasParameter, double)
  itkSetMacro(BiasParameter, double)

  // Get a reference to the valid label state
  itkGetMacro(ValidLabel, bool &)

  // Get a reference to the class index to label mapping
  itkGetMacro(ClassToLabelMapping, MappingType &)

  // Get a reference to the class weights array
  itkGetMacro(ClassWeights, WeightArray &)

  // Set the weight for a class
  void SetClassWeight(size_t class_id, double weight)
    {
    m_ClassWeights[class_id] = weight;
    }

  // Test if the classifier is valid (has 2+ classes)
  bool IsValidClassifier() const
    {
    return m_ClassToLabelMapping.size() >= 2 && m_Forest->GetForestSize() > 0;
    }

protected:

  RandomForestClassifier()
    {
    m_Forest = NULL;
    this->Reset();
    }

  ~RandomForestClassifier()
    {
    if(m_Forest)
      delete m_Forest;
    }

  // The actual decision forest
  RandomForestType *m_Forest;

  // Whether the labels are valid (?)
  bool m_ValidLabel;

  // Mapping of index to label (?)
  MappingType m_ClassToLabelMapping;

  // Weight of each class
  WeightArray m_ClassWeights;

  // The patch radius
  SizeType m_PatchRadius;

  // Whether coordinate features are used
  bool m_UseCoordinateFeatures;

  // Bias parameter
  double m_BiasParameter;

  // Let the engine handle our data
  friend class RFClassificationEngine;
};

#endif // RANDOMFORESTCLASSIFIER_H
