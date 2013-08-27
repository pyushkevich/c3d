#ifndef __itkOrientedRASImage_h_
#define __itkOrientedRASImage_h_

#include "itkImage.h"

namespace itk {

/** 
 * Oriented image with RAS physical coordinates (as opposed to LPS)
 */
template <class TPixel, unsigned int VImageDimension>
class ITK_EXPORT OrientedRASImage : public Image<TPixel, VImageDimension>
{
public:
  /** Standard class typedefs */
  typedef OrientedRASImage               Self;
  typedef Image<TPixel, VImageDimension>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef WeakPointer<const Self>  ConstWeakPointer;
  typedef Matrix<double, VImageDimension+1, VImageDimension+1> TransformMatrixType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OrientedRASImage, Image);

  /** Index typedef support. An index is used to access pixel values. */
  typedef typename Superclass::IndexType  IndexType;

  /** Direction typedef support. The direction cosines of the image. */
  typedef typename Superclass::DirectionType  DirectionType;

  /** Spacing typedef support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
  typedef typename Superclass::SpacingType SpacingType;

  typedef typename Superclass::AccessorType        AccessorType;
  typedef typename Superclass::AccessorFunctorType AccessorFunctorType;
  typedef typename Superclass::IOPixelType         IOPixelType;

  /** Tyepdef for the functor used to access a neighborhood of pixel pointers.*/
  typedef NeighborhoodAccessorFunctor< Self > 
                                            NeighborhoodAccessorFunctorType;

  /** Return the NeighborhoodAccessor functor. This method is called by the 
   * neighborhood iterators. */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() 
    { return NeighborhoodAccessorFunctorType(); }
  
  /** Return the NeighborhoodAccessor functor. This method is called by the 
   * neighborhood iterators. */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
    { return NeighborhoodAccessorFunctorType(); }
  

  /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
  template<class TCoordRep>
  bool TransformRASPhysicalPointToContinuousIndex(
              const Point<TCoordRep, VImageDimension>& point,
              ContinuousIndex<TCoordRep, VImageDimension>& index   ) const
    {
    Point<TCoordRep, VImageDimension> p_lps = point;
    p_lps[0] = -point[0]; p_lps[1] = -point[1];
    return Superclass::TransformPhysicalPointToContinuousIndex(p_lps, index);
    }

  /** Get the index (discrete) from a physical point.
   * Floating point index results are truncated to integers.
   * Returns true if the resulting index is within the image, false otherwise
   * \sa Transform */
  template<class TCoordRep>
  bool TransformRASPhysicalPointToIndex(
            const Point<TCoordRep, VImageDimension>& point,
            IndexType & index                                ) const
    {
    Point<TCoordRep, VImageDimension> p_lps = point;
    p_lps[0] = -point[0]; p_lps[1] = -point[1];
    return Superclass::TransformPhysicalPointToIndex(p_lps, index);
    }

  /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
  template<class TCoordRep>
  void TransformContinuousIndexToRASPhysicalPoint(
            const ContinuousIndex<TCoordRep, VImageDimension>& index,
            Point<TCoordRep, VImageDimension>& point        ) const
    {
    Superclass::TransformContinuousIndexToPhysicalPoint(index, point);
    point[0] = -point[0];
    point[1] = -point[1];
    }

  /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
  template<class TCoordRep>
  void TransformIndexToRASPhysicalPoint(
                      const IndexType & index,
                      Point<TCoordRep, VImageDimension>& point ) const
    {
    Superclass::TransformIndexToPhysicalPoint(index, point);
    point[0] = -point[0];
    point[1] = -point[1];
    }

  /** Take a vector or covariant vector that has been computed in the
   * coordinate system parallel to the image grid and rotate it by the
   * direction cosines in order to get it in terms of the coordinate system of
   * the image acquisition device.  This implementation in the Image
   * multiply the array (vector or covariant vector) by the matrix of Direction
   * Cosines. The arguments of the method are of type FixedArray to make
   * possible to use this method with both Vector and CovariantVector.
   * The Method is implemented differently in the itk::Image.
   *
   * \sa Image
   */ 
  template<class TCoordRep>
  void TransformLocalVectorToRASPhysicalVector(
    const FixedArray<TCoordRep, VImageDimension> & inputGradient,
          FixedArray<TCoordRep, VImageDimension> & outputGradient ) const
    {
    Superclass::TransformLocalVectorToPhysicalVector(inputGradient, outputGradient);
    outputGradient[0] = -outputGradient[0];
    outputGradient[1] = -outputGradient[1];
    }

  /** 
   * Get a matrix that maps points voxel coordinates to RAS coordinates
   */
  TransformMatrixType GetVoxelSpaceToRASPhysicalSpaceMatrix()
    {
    // Generate intermediate terms
    vnl_matrix<double> m_dir, m_ras_matrix;
    vnl_diag_matrix<double> m_scale, m_lps_to_ras;
    vnl_vector<double> v_origin, v_ras_offset;

    // Compute the matrix
    m_dir = this->GetDirection().GetVnlMatrix();
    m_scale.set(this->GetSpacing().GetVnlVector());
    m_lps_to_ras.set(vnl_vector<double>(VImageDimension, 1.0));
    m_lps_to_ras[0] = -1;
    m_lps_to_ras[1] = -1;
    m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

    // Compute the vector
    v_origin = this->GetOrigin().GetVnlVector();
    v_ras_offset = m_lps_to_ras * v_origin;

    // Create the larger matrix
    TransformMatrixType mat;
    vnl_vector<double> vcol(VImageDimension+1, 1.0);
    vcol.update(v_ras_offset);
    mat.SetIdentity();
    mat.GetVnlMatrix().update(m_ras_matrix);
    mat.GetVnlMatrix().set_column(VImageDimension, vcol);

    return mat;
    };

  /** 
   * Set a matrix that maps points voxel coordinates to RAS coordinates
   */
  void SetVoxelSpaceToRASPhysicalSpaceMatrix(vnl_matrix<double> mat)
    {
    // Generate intermediate terms
    vnl_matrix<double> m_dir, m_ras_matrix, m_dist;
    vnl_diag_matrix<double> m_ras_to_lps, m_scale;
    vnl_vector<double> v_origin ;
    vnl_vector<double> m_spacing(VImageDimension, 0.0);

    // Get the dim x dim submatrix from mat
    vnl_matrix<double> smat(VImageDimension,VImageDimension,0.0);
    for (size_t i=0; i< VImageDimension; i++)
      for (size_t j=0; j< VImageDimension; j++)
        smat[i][j] = mat[i][j];
    //smat = mat.get_n_rows(0, VImageDimension).get_n_columns(0, VImageDimension);
    // Get the origin
    m_ras_to_lps.set(vnl_vector<double>(VImageDimension, 1.0));
    m_ras_to_lps[0] = -1;
    m_ras_to_lps[1] = -1;
    vnl_vector<double> v_ras_offset(VImageDimension,0.0);
    v_ras_offset.fill(0.0);
    for (size_t i=0; i< VImageDimension; i++)
      v_ras_offset[i] = mat[i][VImageDimension];
    v_origin = m_ras_to_lps * v_ras_offset;

    // Get the Spacing
    // First, create a matrix of the form [1 0 0; 0 1 0; 0 0 1; 0 0 0] to get distances between consecutive voxels
    // along each axis. When RAS mat will be applied to this matrix, we'll have 3 distance vectors
    vnl_diag_matrix<double> offsetmat(VImageDimension+1, VImageDimension);
    offsetmat.fill(0.0);
    for (size_t i=0; i < VImageDimension+1; i++)
      offsetmat[i]=1.0;
    m_dist = mat * offsetmat;
    // Then compute magnitude of the distance vectors, that's our spacing
    for (size_t i=0; i< VImageDimension; i++)
    {
      vnl_vector<double> distcol(m_dist.get_column(i));
      m_spacing[i] = distcol.magnitude();
    }
    m_scale.set(m_spacing);
    
    // Get the direction
    m_scale.invert_in_place();
    m_dir = m_ras_to_lps * smat * m_scale;
    
    // Set everything
    itk::Matrix<double, VImageDimension, VImageDimension> dir;
    dir.SetIdentity();
    for (size_t i=0; i< VImageDimension; i++)
      for (size_t j=0; j< VImageDimension; j++)
        dir[i][j] = m_dir[i][j];
    this->SetDirection(dir);
    double origin[VImageDimension];
    for (size_t i=0; i< VImageDimension; i++)
      origin[i] = v_origin[i];
    this->SetOrigin(origin);
    double spacing[VImageDimension];
    for (size_t i=0; i< VImageDimension; i++)
      spacing[i] = m_spacing[i];
    this->SetSpacing(spacing);
     
    };

  /** 
   * Get a matrix that maps points in the x * spacing + origin space to
   * the RAS space
   */
  TransformMatrixType GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix()
    {
    TransformMatrixType mat;
    mat.SetIdentity();

    for(size_t i = 0; i < VImageDimension; i++)
      {
      double ras_flip = (i < 2) ? -1 : 1;
      mat[i][VImageDimension] = ras_flip * this->GetOrigin()[i];
      for(size_t j = 0; j < VImageDimension; j++)
        {
        mat[i][j] = ras_flip * this->GetDirection()(i,j) * this->GetSpacing()[i];
        mat[i][VImageDimension] -= ras_flip * this->GetDirection()(i,j) * this->GetOrigin()[i];
        }      
      }

    return mat;
    }



protected:
  OrientedRASImage() {};
  virtual ~OrientedRASImage() {};

private:
  OrientedRASImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} //namespace itk

#endif

