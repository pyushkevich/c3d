#include "ReadImage.h"
#include "itkIOCommon.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

template <class TPixel, unsigned int VDim>
void
ReadImage<TPixel, VDim>
::operator() (const char *file)
{
  // Report
  *c->verbose << "Reading #" << (1 + c->m_ImageStack.size()) << " from " << file << endl;

  // Create an Image IO object using the ITK factory system
  typename itk::ImageIOBase::Pointer iobase = 
    itk::ImageIOFactory::CreateImageIO(file, itk::ImageIOFactory::ReadMode);
  if(!iobase)
    throw ConvertException("Unable to read image %s; IO factory can not create IO object.", file);

  // Read the image information
  iobase->SetFileName(file);
  iobase->ReadImageInformation();

  // Handle SPM origin if SPM flag is specified
  // TODO: get rid of this code, nifti is more reliable
  string ext = itksys::SystemTools::GetFilenameExtension(file);
  if((ext == ".hdr" || ext == ".img.gz" || ext == ".img") && c->m_FlagSPM)
    {
    string temp;
    if(itk::ExposeMetaData<std::string>(
        iobase->GetMetaDataDictionary(), itk::ITK_FileOriginator, temp))
      {
      // Read the SPM-style origin
      *c->verbose << "  Applying SPM origin :";
      for(size_t i=0; i < VDim; i++)
        {
        double sitk = iobase->GetSpacing(i);
        short xspm = (temp[2*i+1] << 8) + temp[2*i];
        double oitk = -sitk * xspm;
        *c->verbose << xspm << " ";
        iobase->SetOrigin(i, oitk);
        }
      *c->verbose << endl;
      }
    }

  // If the image has multiple components, we need to handle it specially
  if(iobase->GetNumberOfComponents() > 1 && c->m_MultiComponentSplit)
    {
    // Read the multi-component image
    typedef itk::VectorImage<TPixel, VDim> VectorImageType;
    typedef itk::ImageFileReader<VectorImageType> VectorReader;
    typename VectorReader::Pointer reader = VectorReader::New();
    reader->SetFileName(file);
    reader->SetImageIO(iobase);
    try { reader->Update(); }
    catch(itk::ExceptionObject &exc)
      {
      throw ConvertException("Error reading image %s\nITK exception: %s", file, exc.GetDescription());
      }

    // Report
    *c->verbose << "  Splitting " << iobase->GetNumberOfComponents() << "-component image." << endl;

    // Split the vector image into component images
    typename VectorImageType::Pointer vec = reader->GetOutput();
    size_t ncomp = vec->GetVectorLength();
    for(size_t i = 0; i < ncomp; i++)
      {
      ImagePointer icomp = ImageType::New();
      icomp->CopyInformation(vec);
      icomp->SetRegions(vec->GetBufferedRegion());
      icomp->Allocate();
      TPixel *src = vec->GetBufferPointer() + i;
      TPixel *dst = icomp->GetBufferPointer();
      TPixel *end = dst + vec->GetBufferedRegion().GetNumberOfPixels();
      for(; dst < end; dst++, src+=ncomp) *dst = *src;

      // Push the image on the stack
      c->m_ImageStack.push_back(icomp);
      }
    }
  else
    {
    // Set up the reader
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(file);
    reader->SetImageIO(iobase);

    // Try reading this file
    try { reader->Update(); }
    catch(itk::ExceptionObject &exc)
      {
      throw ConvertException("Error reading image %s\nITK exception: %s", file, exc.GetDescription());
      }
  
    ImagePointer image = reader->GetOutput();
    c->m_ImageStack.push_back(image);
    }  
}

// Invocations
template class ReadImage<double, 2>;
template class ReadImage<double, 3>;
