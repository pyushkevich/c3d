/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVoxBoCUBImageIO.cxx,v $
  Language:  C++
  Date:      $Date: 2012/11/01 09:19:14 $
  Version:   $Revision: 1.4 $  

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkVoxBoCUBImageIO.h"
#include "itkIOCommon.h"
#include "itkMacro.h"
#include "itkMetaDataObject.h"
#include "itkByteSwapper.h"
#include <iostream>
#include <list>
#include <string>
#include <math.h>

// Commented out because zlib is not available through Insight Applications
#ifdef SNAP_GZIP_SUPPORT
#include <zlib.h>
#endif

inline double myround(double x)
  { return (double) ((int) (x + 0.5)); }

#if !( (ITK_VERSION_MAJOR == 5 && ITK_VERSION_MINOR >= 1) || ITK_VERSION_MAJOR > 5 )
typedef itk::ImageIOBase IOComponentEnum;
#endif

namespace itk {


/**
 * A generic reader and writer object for VoxBo files. Basically it
 * provides uniform access to gzip and normal files
 */
class GenericCUBFileAdaptor
{
public:
  virtual unsigned char ReadByte() = 0;
  virtual void ReadData(void *data, unsigned long bytes) = 0;
  virtual void WriteData(const void *data, unsigned long bytes) = 0;
  virtual ~GenericCUBFileAdaptor() {}

  std::string ReadHeader()
    {
    // Read everything up to the \f symbol
    std::ostringstream oss;
    unsigned char byte = ReadByte();
    while(byte != '\f')
      {
      oss << byte;
      byte = ReadByte();
      }

    // Read the next byte
    unsigned char term = ReadByte();
    if(term == '\r')
      term = ReadByte();

    // Throw exception if term is not there
    if(term != '\n')
      {
      ExceptionObject exception;
      exception.SetDescription("Header is not terminated by newline.");
      throw exception;
      }

    // Return the header string
    return oss.str();
    }
};

/**
 * A reader for gzip files
 */
#ifdef SNAP_GZIP_SUPPORT

class CompressedCUBFileAdaptor : public GenericCUBFileAdaptor
{
public:
  CompressedCUBFileAdaptor(const char *file, const char *mode)
    {
    m_GzFile = ::gzopen(file, mode);
    if(m_GzFile == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be accessed");
      throw exception;
      }
    }
  
  ~CompressedCUBFileAdaptor()
    {
    if(m_GzFile)
      ::gzclose(m_GzFile);
    }

  unsigned char ReadByte()
    {
    int byte = ::gzgetc(m_GzFile);
    if(byte < 0)
      {
      std::ostringstream oss;
      oss << "Error reading byte from file at position: " << ::gztell(m_GzFile);
      ExceptionObject exception;
      exception.SetDescription(oss.str().c_str());
      throw exception;
      }
    return static_cast<unsigned char>(byte);
    }
  
  void ReadData(void *data, unsigned long bytes)
    {
    if(m_GzFile == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be read");
      throw exception;
      }

    int bread = ::gzread(m_GzFile, data, bytes);
    if(bread != bytes)
      {
      std::ostringstream oss;
      oss << "File size does not match header: " 
        << bytes << " bytes requested but only "
        << bread << " bytes available!" << std::endl
        << "At file position " << ::gztell(m_GzFile);
      ExceptionObject exception;
      exception.SetDescription(oss.str().c_str());
      throw exception;
      }
    }
  
  void WriteData(const void *data, unsigned long bytes)
    {
    if(m_GzFile == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be written");
      throw exception;
      }

    int bwritten = ::gzwrite(m_GzFile, (void *) data, bytes);
    if(bwritten != bytes)
      {
      ExceptionObject exception;
      exception.SetDescription("Could not write all bytes to file");
      throw exception;
      }
    }

private:
  ::gzFile m_GzFile;
};

#endif // SNAP_GZIP_SUPPORT

/**
 * A reader for non-gzip files
 */
class DirectCUBFileAdaptor : public GenericCUBFileAdaptor
{
public:
  DirectCUBFileAdaptor(const char *file, const char *mode)
    {
    m_File = fopen(file, mode);
    if(m_File == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be read");
      throw exception;
      }
    }
  
  virtual ~DirectCUBFileAdaptor()
    {
    if(m_File)
      fclose(m_File);
    }

  unsigned char ReadByte()
    {
    int byte = fgetc(m_File);
    if(byte == EOF)
      {
      std::ostringstream oss;
      oss << "Error reading byte from file at position: " << ::ftell(m_File);
      ExceptionObject exception;
      exception.SetDescription(oss.str().c_str());
      throw exception;
      }
    return static_cast<unsigned char>(byte);
    }
  
  void ReadData(void *data, unsigned long bytes)
    {
    if(m_File == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be read");
      throw exception;
      }

    unsigned long bread = fread(data, 1, bytes, m_File);
    if(bread != bytes)
      {
      std::ostringstream oss;
      oss << "File size does not match header: " 
        << bytes << " bytes requested but only "
        << bread << " bytes available!" << std::endl
        << "At file position " << ftell(m_File);
      ExceptionObject exception;
      exception.SetDescription(oss.str().c_str());
      throw exception;
      }
    }

  void WriteData(const void *data, unsigned long bytes)
    {
    if(m_File == NULL)
      {
      ExceptionObject exception;
      exception.SetDescription("File cannot be written");
      throw exception;
      }

    unsigned long bwritten = fwrite(data, 1, bytes, m_File);
    if(bwritten != bytes)
      {
      ExceptionObject exception;
      exception.SetDescription("Could not write all bytes to file");
      throw exception;
      }
    }
private:
  FILE *m_File;
};


/**
 * A swap helper class, used to perform swapping for any input
 * data type.
 */
template<typename TPixel> class VoxBoCUBImageIOSwapHelper
{
public:
#if (ITK_VERSION_MAJOR == 5 && ITK_VERSION_MINOR >= 1) || ITK_VERSION_MAJOR > 5
  typedef ImageIOBase::IOByteOrderEnum ByteOrder;
#else
  typedef ImageIOBase::ByteOrder ByteOrder;
#endif
  static void SwapIfNecessary(
    void *buffer, unsigned long numberOfBytes, ByteOrder order)
    {
    if ( order == ByteOrder::LittleEndian )
      {
      ByteSwapper<TPixel>::SwapRangeFromSystemToLittleEndian(
        (TPixel*)buffer, numberOfBytes / sizeof(TPixel) );
      }
    else if ( order == ByteOrder::BigEndian )
      {
      ByteSwapper<TPixel>::SwapRangeFromSystemToBigEndian(
        (TPixel *)buffer, numberOfBytes / sizeof(TPixel) );
      }
    }
};


// Strings
const char *VoxBoCUBImageIO::VB_IDENTIFIER_SYSTEM = "VB98";
const char *VoxBoCUBImageIO::VB_IDENTIFIER_FILETYPE = "CUB1";
const char *VoxBoCUBImageIO::VB_DIMENSIONS = "VoxDims(XYZ)";
const char *VoxBoCUBImageIO::VB_SPACING = "VoxSizes(XYZ)";
const char *VoxBoCUBImageIO::VB_ORIGIN = "Origin(XYZ)";
const char *VoxBoCUBImageIO::VB_DATATYPE = "DataType";
const char *VoxBoCUBImageIO::VB_BYTEORDER = "Byteorder";
const char *VoxBoCUBImageIO::VB_ORIENTATION = "Orientation";
const char *VoxBoCUBImageIO::VB_BYTEORDER_MSB = "msbfirst";
const char *VoxBoCUBImageIO::VB_BYTEORDER_LSB = "lsbfirst";
const char *VoxBoCUBImageIO::VB_DATATYPE_BYTE = "Byte";
const char *VoxBoCUBImageIO::VB_DATATYPE_INT = "Integer";
const char *VoxBoCUBImageIO::VB_DATATYPE_FLOAT = "Float";
const char *VoxBoCUBImageIO::VB_DATATYPE_DOUBLE = "Double";

/** Constructor */
VoxBoCUBImageIO::VoxBoCUBImageIO()
{
  InitializeOrientationMap();
#if (ITK_VERSION_MAJOR == 5 && ITK_VERSION_MINOR >= 1) || ITK_VERSION_MAJOR > 5
  m_ByteOrder = IOByteOrderEnum::BigEndian;
#else
  m_ByteOrder = BigEndian;
#endif
  m_Reader = NULL;
  m_Writer = NULL;
}


/** Destructor */
VoxBoCUBImageIO::~VoxBoCUBImageIO()
{
  if(m_Reader)
    delete m_Reader;
  if(m_Writer)
    delete m_Writer;
}

GenericCUBFileAdaptor *
VoxBoCUBImageIO::CreateReader(const char *filename)
{
  try
    {
    bool compressed;
    if(CheckExtension(filename, compressed))
      if(compressed)
#ifdef SNAP_GZIP_SUPPORT
          return new CompressedCUBFileAdaptor(filename, "rb");
#else
          return NULL;
#endif
      else
        return new DirectCUBFileAdaptor(filename, "rb");
    else
      return NULL;
    }
  catch(...)
    {
    return NULL;
    }
}

GenericCUBFileAdaptor *
VoxBoCUBImageIO::CreateWriter(const char *filename)
{
  try
    {
    bool compressed;
    if(CheckExtension(filename, compressed))
      if(compressed)
#ifdef SNAP_GZIP_SUPPORT
          return new CompressedCUBFileAdaptor(filename, "rb");
#else
          return NULL;
#endif        
      else
        return new DirectCUBFileAdaptor(filename, "wb");
    else
      return NULL;
    }
  catch(...)
    {
    return NULL;
    }
}

bool VoxBoCUBImageIO::CanReadFile( const char* filename ) 
{ 
  // First check if the file can be read
  GenericCUBFileAdaptor *reader = CreateReader(filename);
  if(reader == NULL)
    {
    itkDebugMacro(<<"The file is not a valid CUB file");
    return false;
    }
    
  // Now check the content
  bool iscub = true;
  try 
    {
    // Get the header
    std::istringstream iss(reader->ReadHeader());

    // Read the first two words
    std::string word;

    // Read the first line from the file
    iss >> word;
    if(word != VB_IDENTIFIER_SYSTEM)
      iscub = false;

    // Read the second line
    iss >> word;
    if(word != VB_IDENTIFIER_FILETYPE)
      iscub = false;
    }
  catch(...)
    { 
    iscub = false; 
    }

  delete reader;
  return iscub;
}

bool VoxBoCUBImageIO::CanWriteFile( const char * name )
{
  bool compressed;
  return CheckExtension(name, compressed);
}

void VoxBoCUBImageIO::Read(void* buffer)
{
  if(m_Reader == NULL)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("File cannot be read");
    throw exception;
    }

  m_Reader->ReadData(buffer, GetImageSizeInBytes());
  this->SwapBytesIfNecessary(buffer, GetImageSizeInBytes());
}

/** 
 *  Read Information about the VoxBoCUB file
 *  and put the cursor of the stream just before the first data pixel
 */
void VoxBoCUBImageIO::ReadImageInformation()
{
  // Make sure there is no other reader
  if(m_Reader)
    delete m_Reader;

  // Create a reader
  m_Reader = CreateReader(m_FileName.c_str());
  if(m_Reader == NULL)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("File cannot be read");
    throw exception;
    }

  // Set the number of dimensions to three
  SetNumberOfDimensions(3);

  // Read the file header
  std::istringstream issHeader(m_Reader->ReadHeader());

  // Read every string in the header. Parse the strings that are special
  while(issHeader.good())
    {
    // Read a line from the stream
    char linebuffer[512];
    issHeader.getline(linebuffer, 512);

    // Get the key string
    std::istringstream iss(linebuffer);
    std::string key;

    // Read the key and strip the colon from it
    iss >> key;
    if(key[key.size() - 1] == ':')
      {
      // Strip the colon off the key
      key = key.substr(0, key.size() - 1);

      // Check if this is a relevant key
      if(key == VB_DIMENSIONS)
        {
        iss >> m_Dimensions[0];
        iss >> m_Dimensions[1];
        iss >> m_Dimensions[2];
        }

      else if(key == VB_SPACING)
        {
        iss >> m_Spacing[0];
        iss >> m_Spacing[1];
        iss >> m_Spacing[2];
        }

      else if(key == VB_ORIGIN)
        {
        double ox, oy, oz;
        iss >> ox; iss >> oy; iss >> oz;
        m_Origin[0] = ox * m_Spacing[0];
        m_Origin[1] = oy * m_Spacing[1];
        m_Origin[2] = oz * m_Spacing[2];
        }

      else if(key == VB_DATATYPE)
        {
        std::string type;
        iss >> type;
#if (ITK_VERSION_MAJOR == 5 && ITK_VERSION_MINOR >= 1) || ITK_VERSION_MAJOR > 5
        m_PixelType = IOPixelEnum::SCALAR;
#else
        m_PixelType = SCALAR;
#endif
        if(type == VB_DATATYPE_BYTE)
          m_ComponentType = IOComponentEnum::UCHAR;
        else if(type == VB_DATATYPE_INT)
          m_ComponentType = IOComponentEnum::USHORT;
        else if(type == VB_DATATYPE_FLOAT)
          m_ComponentType = IOComponentEnum::FLOAT;
        else if(type == VB_DATATYPE_DOUBLE)
          m_ComponentType = IOComponentEnum::DOUBLE;
        }

      else if(key == VB_BYTEORDER)
        {
        std::string type;
        iss >> type;
        if(type == VB_BYTEORDER_MSB)
          SetByteOrderToBigEndian();
        else if(type == VB_BYTEORDER_LSB)
          SetByteOrderToLittleEndian();
        else
          {
          ExceptionObject exception(__FILE__, __LINE__);
          exception.SetDescription("Unknown byte order constant");
          throw exception;
          }
        }

      /*
      else if(key == VB_ORIENTATION)
        {
        std::string code;
        iss >> code;

        // Set the orientation code in the data dictionary
        OrientationMap::const_iterator it = m_OrientationMap.find(code);
        if(it != m_OrientationMap.end())
          {
          itk::MetaDataDictionary &dic =this->GetMetaDataDictionary();
          EncapsulateMetaData<OrientationFlags>(
            dic, ITK_CoordinateOrientation, it->second);
          }
        }
      */

      else
        {
        // Encode the right hand side of the string in the meta-data dic
        std::string word;
        std::ostringstream oss;
        while(iss >> word)
          {
          if(oss.str().size())
            oss << " ";
          oss << word;
          }
        itk::MetaDataDictionary &dic =this->GetMetaDataDictionary();
        EncapsulateMetaData<std::string>(dic, key, oss.str());
        }
      }
    }
}

void 
VoxBoCUBImageIO
::WriteImageInformation(void)
{
  // See if we have a writer already
  if(m_Writer != NULL)
    delete m_Writer;

  // First check if the file can be written to
  m_Writer = CreateWriter(m_FileName.c_str());
  if(m_Writer == NULL)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("File cannot be read");
    throw exception;
    }

  // Check that the number of dimensions is correct
  if(GetNumberOfDimensions() != 3)
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Unsupported number of dimensions");
    throw exception;
    }

  // Put together a header
  std::ostringstream header;

  // Write the identifiers
  header << VB_IDENTIFIER_SYSTEM << std::endl;
  header << VB_IDENTIFIER_FILETYPE << std::endl;

  // Write the image dimensions
  header << VB_DIMENSIONS << ": " 
    << m_Dimensions[0] << " "
    << m_Dimensions[1] << " "
    << m_Dimensions[2] << std::endl;

  // Write the spacing
  header << VB_SPACING << ": "
    << m_Spacing[0] << " "
    << m_Spacing[1] << " "
    << m_Spacing[2] << std::endl;

  // Write the origin (have to convert to bytes)
  header << VB_ORIGIN << ": "
    << (int) myround(m_Origin[0] / m_Spacing[0]) << " "
    << (int) myround(m_Origin[1] / m_Spacing[1]) << " "
    << (int) myround(m_Origin[2] / m_Spacing[2]) << std::endl;

  // Write the byte order
  header << VB_BYTEORDER << ": "
    << (ByteSwapper<char>::SystemIsBigEndian() 
      ? VB_BYTEORDER_MSB : VB_BYTEORDER_LSB) << std::endl;

  // Write the data type 
  switch(m_ComponentType) 
    {
    case IOComponentEnum::CHAR:
    case IOComponentEnum::UCHAR:
      header << VB_DATATYPE << ": " << VB_DATATYPE_BYTE << std::endl;
      break;
    case IOComponentEnum::SHORT:
    case IOComponentEnum::USHORT:
      header << VB_DATATYPE << ": " << VB_DATATYPE_INT << std::endl;
      break;
    case IOComponentEnum::FLOAT:
      header << VB_DATATYPE << ": " << VB_DATATYPE_FLOAT << std::endl;
      break;
    case IOComponentEnum::DOUBLE:
      header << VB_DATATYPE << ": " << VB_DATATYPE_DOUBLE << std::endl;
      break;
    default:
      ExceptionObject exception(__FILE__, __LINE__);
      exception.SetDescription("Unsupported pixel component type");
      throw exception;
    }

  // Write the orientation code
  /*
  MetaDataDictionary &dic = GetMetaDataDictionary();
  OrientationFlags oflag;
  if(ExposeMetaData<OrientationFlags>(dic, ITK_CoordinateOrientation, oflag))
    {
    InverseOrientationMap::const_iterator it = 
      m_InverseOrientationMap.find(oflag);
    if(it != m_InverseOrientationMap.end())
      header << VB_ORIENTATION << ": " << it->second << std::endl;
    }
  */

  // Write the terminating characters
  header << "\f\n";

  // Write the header to the file as data
  m_Writer->WriteData(header.str().c_str(), header.str().size());
}

/** The write function is not implemented */
void 
VoxBoCUBImageIO
::Write( const void* buffer) 
{
  WriteImageInformation();
  m_Writer->WriteData(buffer, GetImageSizeInBytes());
}

/** Print Self Method */
void VoxBoCUBImageIO::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "PixelType " << m_PixelType << "\n";
}


bool VoxBoCUBImageIO::CheckExtension(const char* filename, bool &isCompressed)
{
  std::string fname = filename;
  if ( fname == "" )
  {
    itkDebugMacro(<< "No filename specified.");
    return false;
  }

  bool extensionFound = false;
  isCompressed = false;

  std::string::size_type giplPos = fname.rfind(".cub");
  if ((giplPos != std::string::npos)
      && (giplPos == fname.length() - 4))
    {
      extensionFound = true;
    }

  giplPos = fname.rfind(".cub.gz");
  if ((giplPos != std::string::npos)
      && (giplPos == fname.length() - 7))
    {
    extensionFound = true;
    isCompressed = true;
    }

  return extensionFound;
}

void 
VoxBoCUBImageIO
::InitializeOrientationMap()
{
  m_OrientationMap["RIP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP;
  m_OrientationMap["LIP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP;
  m_OrientationMap["RSP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP;
  m_OrientationMap["LSP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP;
  m_OrientationMap["RIA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA;
  m_OrientationMap["LIA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA;
  m_OrientationMap["RSA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA;
  m_OrientationMap["LSA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA;
  m_OrientationMap["IRP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP;
  m_OrientationMap["ILP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP;
  m_OrientationMap["SRP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP;
  m_OrientationMap["SLP"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP;
  m_OrientationMap["IRA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA;
  m_OrientationMap["ILA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA;
  m_OrientationMap["SRA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA;
  m_OrientationMap["SLA"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA;
  m_OrientationMap["RPI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI;
  m_OrientationMap["LPI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
  m_OrientationMap["RAI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI;
  m_OrientationMap["LAI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI;
  m_OrientationMap["RPS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS;
  m_OrientationMap["LPS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
  m_OrientationMap["RAS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
  m_OrientationMap["LAS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS;
  m_OrientationMap["PRI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI;
  m_OrientationMap["PLI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI;
  m_OrientationMap["ARI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI;
  m_OrientationMap["ALI"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI;
  m_OrientationMap["PRS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS;
  m_OrientationMap["PLS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS;
  m_OrientationMap["ARS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS;
  m_OrientationMap["ALS"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS;
  m_OrientationMap["IPR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR;
  m_OrientationMap["SPR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR;
  m_OrientationMap["IAR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR;
  m_OrientationMap["SAR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR;
  m_OrientationMap["IPL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL;
  m_OrientationMap["SPL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL;
  m_OrientationMap["IAL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL;
  m_OrientationMap["SAL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL;
  m_OrientationMap["PIR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR;
  m_OrientationMap["PSR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR;
  m_OrientationMap["AIR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR;
  m_OrientationMap["ASR"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR;
  m_OrientationMap["PIL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL;
  m_OrientationMap["PSL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL;
  m_OrientationMap["AIL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL;
  m_OrientationMap["ASL"] = SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL;

  OrientationMap::const_iterator it;
  for(it = m_OrientationMap.begin(); it != m_OrientationMap.end(); ++it)
    m_InverseOrientationMap[it->second] = it->first;

}

void 
VoxBoCUBImageIO
::SwapBytesIfNecessary(void *buffer, unsigned long numberOfBytes)
{
  if(m_ComponentType == IOComponentEnum::CHAR)
    VoxBoCUBImageIOSwapHelper<char>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::UCHAR)
    VoxBoCUBImageIOSwapHelper<unsigned char>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::SHORT)
    VoxBoCUBImageIOSwapHelper<short>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::USHORT)
    VoxBoCUBImageIOSwapHelper<unsigned short>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::INT)
    VoxBoCUBImageIOSwapHelper<int>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::UINT)
    VoxBoCUBImageIOSwapHelper<unsigned int>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::LONG)
    VoxBoCUBImageIOSwapHelper<long>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::ULONG)
    VoxBoCUBImageIOSwapHelper<unsigned long>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::FLOAT)
    VoxBoCUBImageIOSwapHelper<float>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else if(m_ComponentType == IOComponentEnum::DOUBLE)
    VoxBoCUBImageIOSwapHelper<double>::SwapIfNecessary(
      buffer, numberOfBytes, m_ByteOrder);
  else 
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Pixel Type Unknown");
    throw exception;
    }
}

} // end namespace itk
