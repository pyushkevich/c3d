#include "UpdateMetadataKey.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template <class TPixel, unsigned int VDim>
void
UpdateMetadataKey<TPixel, VDim>
::operator() (const char *key, const char *value)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Report what we are doing
  *c->verbose << "Updating metadata in #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Setting key " << key << " to value " << value << endl;

  // Update the metadata values
  itk::MetaDataDictionary &mdd = img->GetMetaDataDictionary();
  typedef itk::MetaDataObject<string> StringMetaData;
  typename StringMetaData::Pointer mdval = StringMetaData::New();
  mdval->SetMetaDataObjectValue(value);
  mdd[key] = mdval;
}

// Invocations
template class UpdateMetadataKey<double, 2>;
template class UpdateMetadataKey<double, 3>;
