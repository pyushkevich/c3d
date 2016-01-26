/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    AffineTransformTool.cxx
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

#include "itkOrientedRASImage.h"
#include "itkImageFileReader.h"
#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFactory.h>
#include <itkByteSwapper.h>
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"
#include "vnl/vnl_inverse.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_qr.h"
#include "ConvertException.h"
#include <iostream>

#define RAS_TO_FSL 0
#define FSL_TO_RAS 1

using namespace std;

// Typedefs
typedef itk::OrientedRASImage<double, 3> ImageType;
typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
typedef std::vector<MatrixType> MatrixStack;

typedef vnl_matrix_fixed<double, 3, 3> Mat3;
typedef vnl_vector_fixed<double, 3> Vec3;


int usage()
{
  cout << 
    "RAS Affine Transform Tool"
    "Usage: \n"
    "  c3d_affine_tool [transform_files | options] \n"
    "Options: \n"
    "  -fsl2ras      Convert FSL to RAS\n"
    "  -ras2fsl      Convert RAS to FSL\n"
    "  -ref image    Set reference (fixed) image - only for -fsl2ras and ras2fsl\n"
    "  -src image    Set source (moving) image - only for -fsl2ras and -ras2fsl\n"
    "  -sform image  Read matrix from NIfTI sform\n"
    "  -o matfile    Write output matrix\n"
    "  -det          Print the determinant\n"
    "  -inv          Invert matrix\n"
    "  -mult         Multiply matrices\n"
    "  -sqrt         Matrix square root (i.e., Q s.t. A = Q * Q)\n"
    "  -itk file     Import ITK transform\n"
    "  -oitk file    Export ITK transform\n"
    "  -irtk file    Import IRTK .dof format transform\n"
    "  -oirtk file   Export IRTK .dof format transform\n"
    "  -info         Print matrix\n"
    "  -info-full    Print matrix and more detail about the transform\n"
    "  -rot theta vx vy vz:\n"
    "                Generate rotation matrix corresponding to rotation theta\n"
    "                (in degrees) around vector vx,vy,vz\n"
    "  -trans vx vy vz:\n"
    "                Generate matrix for translation by vx,vy,vz \n"
    "  -scale sx sy sz:\n"
    "                Generate matrix for scaling by sx,sy,sz\n"
    ;
  return -1;
}

void itk_read(MatrixStack &vmat, const char *fname)
{
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> MOTBType;
  typedef itk::AffineTransform<double, 3> AffTran;
  itk::TransformFactory<MOTBType>::RegisterTransform();
  itk::TransformFactory<AffTran>::RegisterTransform();

  itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
  fltReader->SetFileName(fname);
  fltReader->Update();

  itk::TransformBase *base = fltReader->GetTransformList()->front();
  typedef itk::MatrixOffsetTransformBase<double, 3, 3> MOTBType;
  MOTBType *motb = dynamic_cast<MOTBType *>(base);

  MatrixType mat;
  mat.set_identity();
  if(motb)
    {
    for(size_t r = 0; r < 3; r++)
      {
      for(size_t c = 0; c < 3; c++)
        {
        mat(r,c) = motb->GetMatrix()(r,c);
        }
      mat(r,3) = motb->GetOffset()[r];
      }
    mat(2,0) *= -1; mat(2,1) *= -1; 
    mat(0,2) *= -1; mat(1,2) *= -1;
    mat(0,3) *= -1; mat(1,3) *= -1;
    vmat.push_back(mat);
    }
  else
    throw ConvertException("Unable to read ITK transform file %s", fname);
}

void irtk_read(MatrixStack &vmat, const char *fname)
{
  // Read the DOF file
  FILE *fid = fopen(fname, "rb");
  if(!fid)
    throw ConvertException("Unable to read file %s", fname);

  // Read the magic number (or whatever that is)
  unsigned char h[12];
  if(fread(h, 1, 12, fid) < 12)
    throw ConvertException("Unable to read header from file %s", fname);

  // Check magic number
  if(h[1]!=0x0c || h[0]!=0x00 || h[3]!=0x9f || h[2]!=0x6f)
    throw ConvertException("File %s has wrong magic number for DOF file", fname);

  // Check the type of DOF (only 2/3 are supported)
  if(h[7]!=0x02 && h[7]!=0x03)
    throw ConvertException("DOF file %s is not a rigid or affine transform", fname);

  // Allocate parameter array (defaults set below)
  double p[12] = {0,0,0,0,0,0,100,100,100,0,0,0};

  // Read double data
  size_t nval = (h[7] == 0x02) ? 6 : 12;
  if(fread(p, sizeof(double), nval, fid) < nval)
    throw ConvertException("Unable to read data from file %s", fname);

  // Swap bytes if necessary
  itk::ByteSwapper<double>::SwapRangeFromSystemToBigEndian(p, nval);

  // Print the transformation parameters
  printf("DOF parameters: T=(%f, %f, %f); R = (%f, %f, %f); S = (%f, %f, %f); K = (%f, %f, %f)\n",
    p[0],p[1],p[2],p[3],p[4],p[5],
    p[6],p[7],p[8],p[9],p[10],p[11]);

  // Assign to variables
  double *t = p, *r = p+3, *s = p+6, *k = p+9;

  // Close file
  fclose(fid);

  // Initialize matrices
  MatrixType T, R[3], K, S;

  // Create rotation matrices for X, Y, Z
  for(size_t i = 0; i < 3; i++)
    {
    size_t i1 = (i + 1) % 3;
    size_t i2 = (i + 2) % 3;
    R[i].set_identity();
    R[i](i1,i1) = R[i](i2,i2) = cos(r[i] * vnl_math::pi / 180.0);
    R[i](i1,i2) = sin(r[i] * vnl_math::pi / 180.0);
    R[i](i2,i1) = -R[i](i1,i2);
    }

  // Create the translation matrix
  T.set_identity();
  for(size_t i = 0; i < 3; i++)
    T(i,3) = t[i];

  // Create the scale matrix
  S.set_identity();
  for(size_t i = 0; i < 3; i++)
    S(i,i) = 0.01 * s[i];

  // Create the skew matrix
  K.set_identity();
  K(0,1) = tan(k[0] * vnl_math::pi / 180.0);
  K(1,2) = tan(k[1] * vnl_math::pi / 180.0);
  K(0,2) = tan(k[2] * vnl_math::pi / 180.0);

  // Compute the total matrix
  MatrixType M = T * R[0] * R[1] * R[2] * K * S;

  // Push the matrix on the stack
  vmat.push_back(M);
}

void quart_print(MatrixType &mat )
{
  const double epsilon = vcl_numeric_limits<double>::epsilon();
  double m_X, m_Y, m_Z, m_W;

  // Flip the entries that must be flipped to convert to LPS
  mat(2,0) *= -1; mat(2,1) *= -1;
  mat(0,2) *= -1; mat(1,2) *= -1;
  mat(0,3) *= -1; mat(1,3) *= -1;


  typedef vnl_matrix_fixed<double, 3, 3> Mat33;
  Mat33 A = mat.extract(3,3);


  // QR decomposition
  vnl_qr<double> qr(A);
  vnl_matrix_fixed<double, 3, 3> m = qr.Q(), R = qr.R(), F;

  F.set_identity();
  for(size_t i = 0; i < 3; i++)
    if(R(i,i) < 0)
      F(i,i) = -1;

  // Scale Q and R by the flip matrix
  m = m * F;
  R = F * R;


  printf("Rotation matrix:\n");
  for(size_t i = 0; i < 3; i++)
    printf("%12.5f   %12.5f   %12.5f\n", m(i,0), m(i,1), m(i,2));

  double r[3], s[3], k[3];
  // Get the rotation angles
  r[0] = atan2(m(1,2), m(2,2)) * 180. / vnl_math::pi;
  r[1] = asin(-m(0,2)) * 180. / vnl_math::pi;
  r[2] = atan2(m(0,1), m(0,0)) * 180. / vnl_math::pi;
   
  // Get the scales
  for(size_t i = 0; i < 3; i++)
    s[i] = 100. * R(i,i);

  // Get the shears
  vnl_matrix_fixed<double, 3, 3> Sinv; Sinv.set_identity();
  for(size_t i = 0; i < 3; i++)
    Sinv(i,i) = 1.0 / R(i,i);
  vnl_matrix_fixed<double, 3, 3> K = R * Sinv;

  k[0] = atan(K(0,1)) * 180. / vnl_math::pi;
  k[1] = atan(K(1,2)) * 180. / vnl_math::pi;
  k[2] = atan(K(0,2)) * 180. / vnl_math::pi;
  
  printf("Affine parameters:  T=(%f, %f, %f); R = (%f, %f, %f); S = (%f, %f, %f); K = (%f, %f, %f)\n",
          mat(0,3), mat(1,3), mat(2,3), r[0], r[1], r[2], s[0], s[1], s[2], k[0], k[1], k[2]);
   


  double trace = m(0,0) + m(1,1) + m(2,2) + 1.0;
//std::cout << "trace: " << trace << " epsilon: " << vcl_numeric_limits<T>::epsilon() << std::endl;

  if( trace > epsilon)
    {
    const double s = 0.5 / vcl_sqrt(trace);
    m_W = 0.25 / s;
    m_X = (m(2,1) - m(1,2)) * s;
    m_Y = (m(0,2) - m(2,0)) * s;
    m_Z = (m(1,0) - m(0,1)) * s;
//std::cout << "opt 1: w " << m_W << " x " << m_X << " y " << m_Y << " z " << m_Z << std::endl;
    }
  else
    {
    if( m(0,0) > m(1,1) && m(0,0) > m(2,2) )
      {
      const double s = 2.0 * vcl_sqrt(1.0 + m(0,0) - m(1,1) - m(2,2));
      m_X = 0.25 * s;
      m_Y = (m(0,1) + m(1,0)) / s;
      m_Z = (m(0,2) + m(2,0)) / s;
      m_W = (m(1,2) - m(2,1)) / s;
//std::cout << "opt 2: w " << m_W << " x " << m_X << " y " << m_Y << " z " << m_Z << std::endl;
      }
    else
      {
      if( m(1,1) > m(2,2) )
        {
        const double s = 2.0 * vcl_sqrt(1.0 + m(1,1) - m(0,0) - m(2,2));
        m_X = (m(0,1) + m(1,0)) / s;
        m_Y = 0.25 * s;
        m_Z = (m(1,2) + m(2,1)) / s;
        m_W = (m(0,2) - m(2,0)) / s;
//std::cout << "opt 3: w " << m_W << " x " << m_X << " y " << m_Y << " z " << m_Z << std::endl;
        }
      else
        {
        const double s = 2.0 * vcl_sqrt(1.0 + m(2,2) - m(0,0) - m(1,1));
        m_X = (m(0,2) + m(2,0)) / s;
        m_Y = (m(1,2) + m(2,1)) / s;
        m_Z = 0.25 * s;
        m_W = (m(0,1) - m(1,0)) / s;
//std::cout << "opt 4: w " << m_W << " x " << m_X << " y " << m_Y << " z " << m_Z << std::endl;
        }
      }
    }
  double mag = vcl_sqrt( m_X*m_X + m_Y*m_Y + m_Z*m_Z + m_W*m_W );
  m_X /= mag;
  m_Y /= mag;
  m_Z /= mag;
  m_W /= mag;

  printf("Quaternion:\n");
  printf("%12.5f   %12.5f   %12.5f  %12.5f\n", m_X, m_Y, m_Z, m_W);

  double L = vcl_sqrt( m_X*m_X + m_Y*m_Y + m_Z*m_Z );
  double angle = (180.0/vnl_math::pi) * 2.0 * asin( L );
  double axis[3];
  axis[0] = m_X/L;
  axis[1] = m_Y/L;
  axis[2] = m_Z/L;
  
  printf("Rotation angle:\n");
  printf("%f degrees\n", angle); 

  printf("Rotation axis:\n");
  printf("%12.5f   %12.5f   %12.5f\n", axis[0], axis[1], axis[2] );

}

void irtk_write(MatrixStack &vmat, const char *fname)
{
  // Get the current matrix
  MatrixType M = vmat.back();

  // Initialize the parameters to save
  double p[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  double *t = p, *r = p+3, *s = p+6, *k = p+9;

  // Set the translation parameters
  for(size_t i = 0; i < 3; i++)
    t[i] = M(i,3);

  // Extract the 3x3 affine matrix
  typedef vnl_matrix_fixed<double, 3, 3> Mat33;
  Mat33 A = M.extract(3,3);
  
  // QR decomposition
  vnl_qr<double> qr(A);
  Mat33 Q = qr.Q(), R = qr.R();

  // Compute the flip matrix so that the scales are positive
  Mat33 F; F.set_identity();
  for(size_t i = 0; i < 3; i++)
    if(R(i,i) < 0)
      F(i,i) = -1;

  // Scale Q and R by the flip matrix
  Q = Q * F;
  R = F * R;

  // Get the rotation angles
  r[0] = atan2(Q(1,2), Q(2,2)) * 180. / vnl_math::pi;
  r[1] = asin(-Q(0,2)) * 180. / vnl_math::pi;
  r[2] = atan2(Q(0,1), Q(0,0)) * 180. / vnl_math::pi;
  
  // Get the scales
  for(size_t i = 0; i < 3; i++)
    s[i] = 100. * R(i,i);

  // Get the shears
  Mat33 Sinv; Sinv.set_identity();
  for(size_t i = 0; i < 3; i++)
    Sinv(i,i) = 1.0 / R(i,i);
  Mat33 K = R * Sinv;

  k[0] = atan(K(0,1)) * 180. / vnl_math::pi;
  k[1] = atan(K(1,2)) * 180. / vnl_math::pi;
  k[2] = atan(K(0,2)) * 180. / vnl_math::pi;


  // Print the transformation parameters
  printf("DOF parameters: T=(%f, %f, %f); R = (%f, %f, %f); S = (%f, %f, %f); K = (%f, %f, %f)\n",
    p[0],p[1],p[2],p[3],p[4],p[5],
    p[6],p[7],p[8],p[9],p[10],p[11]);

  // Write to DOF file
  const char *h = "\x00\x0c\x6f\x9f\x00\x00\x00\x03\x00\x00\x00\x0c";
  FILE *fid = fopen(fname, "wb");
  if(!fid)
    throw ConvertException("Unable to open file %s for writing", fname);

  // Write the header
  fwrite(h, 1, 12, fid);

  // Write the data
  itk::ByteSwapper<double>::SwapRangeFromSystemToBigEndian(p, 12);
  fwrite(p, sizeof(double), 12, fid);
  fclose(fid);
}

void itk_write(MatrixStack &vmat, const char *fname)
{
  // Get the current matrix
  MatrixType mat = vmat.back();
  
  // Flip the entries that must be flipped
  mat(2,0) *= -1; mat(2,1) *= -1; 
  mat(0,2) *= -1; mat(1,2) *= -1;
  mat(0,3) *= -1; mat(1,3) *= -1;

  // Create an ITK affine transform
  typedef itk::MatrixOffsetTransformBase<double, 3> AffTran;
  AffTran::Pointer atran = AffTran::New();

  // Populate its matrix
  AffTran::MatrixType amat = atran->GetMatrix();
  AffTran::OffsetType aoff = atran->GetOffset();

  for(size_t r = 0; r < 3; r++)
    {
    for(size_t c = 0; c < 3; c++)
      {
      amat(r,c) = mat(r,c);
      }
    aoff[r] = mat(r,3);
    }

  atran->SetMatrix(amat);
  atran->SetOffset(aoff);

  // Write the transform
  itk::TransformFileWriter::Pointer wrt = itk::TransformFileWriter::New();
  wrt->SetInput(atran);
  wrt->SetFileName(fname);
  wrt->Update();
}

void ras_read(MatrixStack &vmat, const char *fname)
{
  MatrixType mat;

  ifstream fin(fname);
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw "Unable to read matrix";
        }
  fin.close();

  vmat.push_back(mat);
}

void ras_write(MatrixStack &vmat, const char *fname)
{
  MatrixType mat = vmat.back();

  ofstream fout(fname);
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      fout << mat[i][j] << (j < 3 ? " " : "\n");

  fout.close();
}

void fsl_to_ras(MatrixStack &vmat, ImageType *ref, ImageType *mov, short flag)
{
  MatrixType m_fsl, m_spcref, m_spcmov, m_swpref, m_swpmov, m_ref, m_mov, m_out;
  m_fsl = vmat.back();

  // Set the ref/mov matrices
  m_ref = ref->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  m_mov = mov->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();

  // Set the swap matrices
  m_swpref.set_identity();
  if(vnl_det(m_ref) > 0)
    {
    m_swpref(0,0) = -1.0;
    m_swpref(0,3) = (ref->GetBufferedRegion().GetSize(0) - 1) * ref->GetSpacing()[0];
    }

  m_swpmov.set_identity();
  if(vnl_det(m_mov) > 0)
    {
    m_swpmov(0,0) = -1.0;
    m_swpmov(0,3) = (mov->GetBufferedRegion().GetSize(0) - 1) * mov->GetSpacing()[0];
    }

  // Set the spacing matrices
  m_spcref.set_identity();
  m_spcmov.set_identity();
  for(size_t i = 0; i < 3; i++)
    {
    m_spcref(i,i) = ref->GetSpacing()[i];
    m_spcmov(i,i) = mov->GetSpacing()[i];
    }

  // Compute the output matrix
  if (flag == FSL_TO_RAS)
    m_out = 
    m_mov * vnl_inverse(m_spcmov) * m_swpmov *
    vnl_inverse(m_fsl) * 
    m_swpref * m_spcref * vnl_inverse(m_ref);

  // NOTE: m_fsl is really m_ras here
  if (flag == RAS_TO_FSL)
    m_out =
    vnl_inverse(vnl_inverse(m_swpmov) * m_spcmov* vnl_inverse(m_mov) *
    m_fsl *
    m_ref*vnl_inverse(m_spcref)*vnl_inverse(m_swpref));

  // Put it on the stack
  vmat.pop_back();
  vmat.push_back(m_out);
}

void ras_inv(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  MatrixType minv = vnl_inverse(m);
  vmat.pop_back();
  vmat.push_back(minv);
}

void ras_sqrt(MatrixStack &vmat)
{
  MatrixType m = vmat.back();

  // Peform Denman-Beavers iteration
  MatrixType Z, Y = m;
  Z.set_identity();

  for(size_t i = 0; i < 16; i++)
    {
    MatrixType Ynext = 0.5 * (Y + vnl_inverse(Z));
    MatrixType Znext = 0.5 * (Z + vnl_inverse(Y));
    Y = Ynext;
    Z = Znext;
    }

  vmat.pop_back();
  vmat.push_back(Y);
}

void make_rotation(MatrixStack &vmat, double theta, Vec3 v)
{
  // Normalize the vector by magnitude
  v.normalize();

  // Convert the angle to radians
  double theta_rad = theta * vnl_math::pi / 180;

  // Compute the skew-symmetric matrix
  Mat3 S; S.fill(0.0);
  S(0,1) = -v[2]; S(1,0) = v[2];
  S(2,0) = -v[1]; S(0,2) = v[1];
  S(1,2) = -v[0]; S(2,1) = v[0];

  // Apply Rodriguez formula
  Mat3 R; R.set_identity();
  R += sin(theta_rad) * S;
  R += (1 - cos(theta_rad)) * (S * S);

  // Fill out the complete matrix
  MatrixType A; A.set_identity();
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      A(i,j) = R(i,j);

  vmat.push_back(A);
}


void make_translation(MatrixStack &vmat, Vec3 v)
{
  // Fill out the complete matrix
  MatrixType A; A.set_identity();
  for(int i = 0; i < 3; i++)
    A(i,3) = v[i];

  vmat.push_back(A);
}

void make_scaling(MatrixStack &vmat, Vec3 v)
{
  // Fill out the complete matrix
  MatrixType A; A.set_identity();
  for(int i = 0; i < 3; i++)
    A(i,i) = v[i];

  vmat.push_back(A);
}


void ras_det(MatrixStack &vmat)
{
  MatrixType m = vmat.back();
  double det = vnl_det(m);
  cout << "Det: " << det << endl;
}

void ras_sform(MatrixStack &vmat, const char *fname)
{
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fname);
  reader->Update();

  ImageType::Pointer img = reader->GetOutput();
 
  // Set up the directions
  MatrixType m = img->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  vmat.push_back(m);
}


void ras_print(MatrixStack &vmat, bool isfullinfo)
{
  MatrixType m = vmat.back();
  printf("Matrix #%d:\n", (int) vmat.size());
  for(size_t i = 0; i < 4; i++)
    printf("%12.5f   %12.5f   %12.5f   %12.5f\n", m(i,0), m(i,1), m(i,2), m(i,3));

  if ( isfullinfo )
    {

    quart_print(m);
    }

}

void ras_mult(MatrixStack &vmat)
{
  MatrixType A = vmat[vmat.size() - 2];
  MatrixType B = vmat[vmat.size() - 1];
  MatrixType AB = A * B;
  vmat.pop_back();
  vmat.pop_back();
  vmat.push_back(AB);
}

int main(int argc, char *argv[])
{
  // Show usage
  if(argc < 2) return usage();

  // Set up the images that might be loaded
  ImageType::Pointer ref = NULL, src = NULL;

  // Set up the matrix stack
  vector<MatrixType> vmat;

  // Parse the command line
  for(int iarg = 1; iarg < argc; iarg++)
    {
    try
      {
      string arg = argv[iarg];
      if(arg == "-ref")
        {
        // Read the reference image
        typedef itk::ImageFileReader<ImageType> ReaderType;
        ReaderType::Pointer read_ref = ReaderType::New();
        read_ref->SetFileName(argv[++iarg]);
        read_ref->Update();
        ref = read_ref->GetOutput();
        }
      else if(arg == "-src" || arg == "-mov")
        {
        typedef itk::ImageFileReader<ImageType> ReaderType;
        // Read the moving image
        ReaderType::Pointer read_src = ReaderType::New();
        read_src->SetFileName(argv[++iarg]);
        read_src->Update();
        src = read_src->GetOutput();
        }
      else if(arg == "-fsl2ras")
        {
        // Convert FSL to RAS 
        fsl_to_ras(vmat, ref, src, FSL_TO_RAS);
        }
      else if(arg == "-ras2fsl")
        {
        // Convert FSL to RAS 
        fsl_to_ras(vmat, ref, src, RAS_TO_FSL);
        }
      else if(arg == "-mult")
        {
        ras_mult(vmat);
        }
      else if(arg == "-det")
        {
        ras_det(vmat);
        }
      else if(arg == "-inv")
        {
        ras_inv(vmat);
        }
      else if(arg == "-sqrt")
        {
        ras_sqrt(vmat);
        }
      else if(arg == "-info")
        {
        ras_print(vmat, false);
        }
      else if(arg == "-info-full")
        {
        ras_print(vmat, true);
        }
      else if(arg == "-sform")
        {
        ras_sform(vmat, argv[++iarg]);
        }
      else if(arg == "-itk")
        {
        itk_read(vmat, argv[++iarg]);
        }
      else if(arg == "-irtk")
        {
        irtk_read(vmat, argv[++iarg]);
        }
      else if(arg == "-o")
        {
        ras_write(vmat, argv[++iarg]);
        }
      else if(arg == "-oitk")
        {
        itk_write(vmat, argv[++iarg]);
        }
      else if(arg == "-oirtk")
        {
        irtk_write(vmat, argv[++iarg]);
        }
      else if(arg == "-rot")
        {
        double theta = atof(argv[++iarg]);
        Vec3 v;
        v[0] = atof(argv[++iarg]);
        v[1] = atof(argv[++iarg]);
        v[2] = atof(argv[++iarg]);
        make_rotation(vmat, theta, v);
        }
      else if(arg == "-tran")
        {
        Vec3 v;
        v[0] = atof(argv[++iarg]);
        v[1] = atof(argv[++iarg]);
        v[2] = atof(argv[++iarg]);
        make_translation(vmat, v);
        }
      else if(arg == "-scale")
        {
        Vec3 v;
        v[0] = atof(argv[++iarg]);
        v[1] = atof(argv[++iarg]);
        v[2] = atof(argv[++iarg]);
        make_scaling(vmat, v);
        }
      else if(arg[0] != '-')
        {
        ras_read(vmat, arg.c_str());
        }
      else
        {
        cerr << "Unknown option " << arg << endl;
        return usage();
        }
      }
    catch(std::exception &exc)
      {
      cerr << "Exception raised during processing:" << endl;
      cerr << exc.what() << endl;
      return -1;
      }
    }
  
}
