/*
 * cl_FEM_Geometry_Interpolator.hpp
 *
 *  Created on: Jan 31, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_

#include "typedefs.hpp" //MRS/COR/src

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"              //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src

#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp" //LINALG/src
#include "op_times.hpp"  //LINALG/src
#include "op_plus.hpp"   //LINALG/src
#include "op_minus.hpp"  //LINALG/src
#include "fn_trans.hpp"  //LINALG/src
#include "fn_det.hpp"    //LINALG/src
#include "fn_inv.hpp"    //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    /**
     * \brief a special interpolation class for geometry
     */
    class Geometry_Interpolator
    {
        // pointer to space interpolation function object
        Interpolation_Function_Base * mSpaceInterpolation = nullptr;

        // pointer to time interpolation function object
        Interpolation_Function_Base * mTimeInterpolation = nullptr;

        // number of space bases, number of physical and parametric dimensions
        uint mNumSpaceBases;
        uint mNumSpaceDim;
        uint mNumSpaceParamDim;

        // number of time bases and dimensions
        uint mNumTimeBases;
        uint mNumTimeDim;

        // matrix of space coefficients xHat
        Matrix < DDRMat > mXHat;

        // matrix of time coefficients tHat
        Matrix < DDRMat > mTHat;

        // matrix of element parametric coordinates
        Matrix< DDRMat > mParamCoords;

        // matrix of element physical coordinates
        Matrix< DDRMat > mPhysCoords;

        // element geometry type
        mtk::Geometry_Type mGeometryType;

        // boolean true if side interpolation
        bool mSpaceSideset = false;

        // side geometry type
        mtk::Geometry_Type mSideGeometryType;

        // pointer to side space interpolation function object
        Interpolation_Function_Base * mSideSpaceInterpolation = nullptr;

        // number of space bases, number of physical and parametric dimensions for the side
        uint mSideNumSpaceBases;
        uint mSideNumSpaceDim;
        uint mSideNumSpaceParamDim;

        // Vertices ordinals for each face of the parent element
        moris::Cell< moris::Cell< moris::moris_index > > mVerticesOrdinalsPerFace;

        // pointer to function for second derivative
        void ( * mSecondDerivativeMatricesSpace )( const Matrix< DDRMat > & aJt,
                                                         Matrix< DDRMat > & aKt,
                                                         Matrix< DDRMat > & aLt,
                                                   const Matrix< DDRMat > & ad2NdXi2,
                                                   const Matrix< DDRMat > & aXHat );

        void ( * mSecondDerivativeMatricesTime )( const Matrix< DDRMat > & aJt,
                                                        Matrix< DDRMat > & aKt,
                                                        Matrix< DDRMat > & aLt,
                                                  const Matrix< DDRMat > & ad2NdTau2,
                                                  const Matrix< DDRMat > & aTHat );

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] interpolation rule for geometry
         * @param[ in ] flag true if side interpolation
         */
        Geometry_Interpolator( const Interpolation_Rule & aInterpolationRule,
                               const bool                 aSpaceSideset = false );

//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Geometry_Interpolator();

//------------------------------------------------------------------------------
        /**
         * returns the order of the space interpolation
         */
        mtk::Interpolation_Order get_space_interpolation_order() const
        {
            return mSpaceInterpolation->get_interpolation_order();
        }

//------------------------------------------------------------------------------
        /**
         * returns the order of the time interpolation
         */
        mtk::Interpolation_Order get_time_interpolation_order() const
        {
            return mTimeInterpolation->get_interpolation_order();
        }

//------------------------------------------------------------------------------
        /**
         * returns the type of the space interpolation
         */
        Interpolation_Type get_space_interpolation_type() const
        {
            return mSpaceInterpolation->get_interpolation_type();
        }

//------------------------------------------------------------------------------
        /**
         * returns the type of the time interpolation
         */
        Interpolation_Type get_time_interpolation_type() const
        {
            return mTimeInterpolation->get_interpolation_type();
        }

//------------------------------------------------------------------------------
        /**
         * returns the number of space dimensions
         */
         uint get_number_of_space_dimensions() const
         {
             return mSpaceInterpolation->get_number_of_dimensions();
         }

//------------------------------------------------------------------------------
         /**
          * returns the number of time dimensions
          */
          uint get_number_of_time_dimensions() const
          {
              return mTimeInterpolation->get_number_of_dimensions();
          }

//------------------------------------------------------------------------------
         /**
          * returns the number of bases for the space function
          */
          uint get_number_of_space_bases() const
          {
              return mSpaceInterpolation->get_number_of_bases();
          }

//------------------------------------------------------------------------------
          /**
           * returns the number of bases for the space function
           */
           uint get_number_of_time_bases() const
           {
               return mTimeInterpolation->get_number_of_bases();
           }

//------------------------------------------------------------------------------
          /**
           * returns the side geometry type
           */
           mtk::Geometry_Type get_side_geometry_type() const
           {
               return mSideGeometryType;
           }

//------------------------------------------------------------------------------
        /**
         * set the space and time coefficients of the geometry field xHat, tHat
         * @param[ in ] space coefficients
         * @param[ in ] time coefficients
         */
        void set_coeff( const Matrix< DDRMat > & aXHat,
                        const Matrix< DDRMat > & aTHat );

//------------------------------------------------------------------------------
        /**
         * set the space coefficients of the geometry field xHat
         * @param[ in ] space coefficients
         */
        void set_space_coeff( const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------
        /**
         * set the time coefficients of the geometry field tHat
         * @param[ in ] time coefficients
         */
        void set_time_coeff( const Matrix< DDRMat > & aTHat );

//------------------------------------------------------------------------------
         /**
          * get the space coefficients of the geometry field xHat
          */
          Matrix< DDRMat > get_space_coeff() const
          {
              return mXHat;
          }

//------------------------------------------------------------------------------
          /**
           * get the time coefficients of the geometry field tHat
           */
           Matrix< DDRMat > get_time_coeff() const
           {
               return mTHat;
           }

//------------------------------------------------------------------------------
        /**
         * get the space time parametric coordinates
         */
        void get_space_time_param_coords();

//------------------------------------------------------------------------------
        /**
         * get the space time physical coordinates
         */
        void get_space_time_phys_coords();

//------------------------------------------------------------------------------
        /**
         * get the vertices ordinals of each face of the parent element
         */
        void get_face_vertices_ordinals();

//------------------------------------------------------------------------------
        /**
         * get the parametric coordinates of a space side
         * @param[ in ] a space face ordinal
         */
        Matrix< DDRMat > get_space_side_param_coords( const moris_index aSpaceOrdinal );

//------------------------------------------------------------------------------
        /**
         * get the parametric coordinates of a space side
         * @param[ in ] a space face ordinal
         */
        Matrix< DDRMat > get_space_side_phys_coords( const moris_index aSpaceOrdinal );

//------------------------------------------------------------------------------
        /**
         * get the parametric coordinates of a time side
         * @param[ in ] a time face ordinal
         */
        Matrix< DDRMat > get_time_side_param_coords( const moris_index aTimeOrdinal );

//------------------------------------------------------------------------------
        /**
         * evaluates the space shape functions at a given evaluation point
         * @param[ out ] aNXi shape functions ( 1 x <number of nodes> )
         * @param[ in ]  aXi  evaluation point ( <number of dimensions> x 1 )
         */
        Matrix < DDRMat > NXi( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the time shape functions at a given evaluation point
         * @param[ out ] aNTau shape functions ( 1 x <number of nodes> )
         * @param[ in ]  aTau  evaluation point ( <number of dimensions> x 1 )
         */
        Matrix < DDRMat > NTau( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the first derivatives of the space shape functions
         * wrt parametric coordinates at a given evaluation point
         * @param[ out ] adNdXi derivatives
         *                      ( <number of dimensions> x <number of nodes> )
         * @param[ in ]  aXi    evaluation point
         *                      ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > dNdXi( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the first derivatives of the time shape functions
         * wrt parametric coordinates at a given evaluation point
         * @param[ out ] adNdTau first order derivatives
         *                       ( <number of dimensions> x <number of nodes> )
         * @param[ in ] aTau     evaluation point
         *                       ( <number of dimensions>  x 1 )
         */
         Matrix< DDRMat > dNdTau( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the second derivatives of the space shape functions
         *  wrt parametric coordinates at a given evaluation point
         * @param[ out ] ad2NdXi2 second order derivatives
         *                        ( <number of dimensions> x <number of nodes> )
         * @param[ in ] aXi       evaluation point
         *                        ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > d2NdXi2 ( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the second derivatives of the time shape functions
         *  wrt parametric coordinates at a given evaluation point
         * @param[ out ] ad2NdTau2 second order derivatives
         *                         ( <number of dimensions> x <number of nodes> )
         * @param[ in ] aTau       evaluation point
         *                         ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > d2NdTau2 ( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian in space
         * @param[ out ] tJt    transposed of geometry Jacobian in space
         * @param[ in ]  adNdXi first derivatives of space shape functions in
         *                      parameter space
         */
        Matrix< DDRMat > space_jacobian( const Matrix< DDRMat > & adNdXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian in time
         * @param[ out ] tJt     transposed of geometry Jacobian in time
         * @param[ in ]  adNdTau first derivatives of time shape functions in
         *                       parameter space
         */
        Matrix< DDRMat > time_jacobian( const Matrix< DDRMat > & adNdTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the determinant of the Jacobian mapping
         * at given space and time evaluation point
         * @param[ in ]  aParamPoint evaluation point
         */
        real det_J( const Matrix< DDRMat > & aParamPoint );

//------------------------------------------------------------------------------
        /**
         * evaluates the determinant of the Jacobian mapping
         * in the case of a time side interpolation
         * at given space and time evaluation point
         * @param[ in ]  aSideParamPoint evaluation point on the side
         * @param[ in ]  aTimeOrdinal    a time face ordinal
         */
         real time_surf_det_J( const Matrix< DDRMat > & aSideParamPoint,
                               const moris_index      & aTimeOrdinal );

//------------------------------------------------------------------------------
         /**
          * get the parametric coordinates of a point
          * in the side parametric space
          * in the parent parametric space
          */
         Matrix< DDRMat > surf_val( const Matrix< DDRMat > & aSideParamPoint,
                                    const moris_index      & aSideOrdinal );

//------------------------------------------------------------------------------
         /**
          * get the parametric coordinates of a point
          * in the side parametric space
          * in the parent parametric space
          */
         Matrix < DDRMat > time_surf_val( const Matrix< DDRMat > & aSideParamPoint,
                                          const moris_index      & aTimeOrdinal );

//------------------------------------------------------------------------------
         /**
          * evaluates the determinant of the Jacobian mapping
          * in the case of a space side interpolation
          * at given space and time evaluation point
          * @param[ out ] aTimeSurfDetJ   determinant of the Jacobian
          * @param[ out ] aNormal         normal to the face
          * @param[ in ]  aSideParamPoint evaluation point on the face
          * @param[ in ]  aSpaceOrdinal   a space face ordinal
          */
         void surf_det_J(       real             & aSurfDetJ,
                                Matrix< DDRMat > & aNormal,
                          const Matrix< DDRMat > & aSideParamPoint,
                          const moris_index      & aSpaceOrdinal );

//------------------------------------------------------------------------------
         /**
          * evaluates the determinant of the Jacobian mapping
          * in the case of a space side interpolation
          * at given space and time evaluation point
          * @param[ in ]  aSideParamPoint evaluation point on the face
          * @param[ in ]  aSpaceOrdinal   a space face ordinal
          */
         real surf_det_J_new( const Matrix< DDRMat > & aSideParamPoint,
                              const moris_index      & aSpaceOrdinal );

//------------------------------------------------------------------------------
         /**
          * evaluates the normal to the side
          * in the case of a space side interpolation
          * at given space and time evaluation point
          * @param[ in ]  aSideParamPoint evaluation point on the face
          * @param[ in ]  aSpaceOrdinal   a space face ordinal
          */
         Matrix< DDRMat > surf_normal( const Matrix< DDRMat > & aSideParamPoint,
                                       const moris_index      & aSpaceOrdinal );

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derivatives wrt to space in physical space
         * @param[ in ]  aJt      transposed of geometry Jacobian
         * @param[ out ] aKt      matrix for second derivatives in physical space
         * @param[ out ] aLt      matrix for second derivatives in physical space
         * @param[ in ]  adNdXi   first derivatives of N in parameter space
         * @param[ in ]  ad2NdXi2 second derivatives of N in parameter space
         *
         */
        void space_jacobian_and_matrices_for_second_derivatives(       Matrix< DDRMat > & aJt,
                                                                       Matrix< DDRMat > & aKt,
                                                                       Matrix< DDRMat > & aLt,
                                                                 const Matrix< DDRMat > & adNdXi,
                                                                 const Matrix< DDRMat > & ad2NdXi2 ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derivatives wrt to time in physical space
         * @param[ in ]  aJt       transposed of geometry Jacobian
         * @param[ out ] aKt       matrix for second derivatives in physical space
         * @param[ out ] aLt       matrix for second derivatives in physical space
         * @param[ in ]  adNdTau   first derivatives of N in parameter space
         * @param[ in ]  ad2NdTau2 second derivatives of N in parameter space
         *
         */
        void time_jacobian_and_matrices_for_second_derivatives(       Matrix< DDRMat > & aJt,
                                                                      Matrix< DDRMat > & aKt,
                                                                      Matrix< DDRMat > & aLt,
                                                                const Matrix< DDRMat > & adNdTau,
                                                                const Matrix< DDRMat > & ad2NdTau2 ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the space geometry field at a given evaluation point in space
         * @param[ out ] x   location in space
         * @param[ in ]  aXi evaluation point in space
         */
        Matrix< DDRMat > valx( const Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------
        /**
         * evaluates the space geometry field at a given evaluation point in time
         * @param[ out ] t   location in time
         * @param[ in ]  aTau evaluation point in time
         */
        Matrix< DDRMat > valt( const Matrix< DDRMat > & aTau );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
        /**
         * sets the function pointers for 2D and 3D. Called during construction.
         */
        void set_function_pointers();

//------------------------------------------------------------------------------
        /**
         * evaluates matrices that are needed for the second derivative
         * in space, 2D version. It is
         *
         * \f[
         *      \mathbf{L}^T \, \mathbf{\frac{\partial^2 N}{\partial x^2}}
         *      = \mathbf{\frac{\partial^2 N}{\partial \xi^2}}
         *      - K^T \, mathbf{\frac{\partial N}{\partial x}}
         * \f]
         *
         * @param[ in ]  aJt          transposed of geometry Jacobian
         * @param[ out ] aKt          transposed help matrix K
         * @param[ out ] aLt          transposed help matrix L
         * @param[ in ]  adNdXi       first derivatives in parameter space
         * @param[ in ]  ad2NdX2i     second derivatives in parameter space
         *
         */
        static void eval_matrices_for_second_derivative_1d( const Matrix< DDRMat > & aJt,
                                                                  Matrix< DDRMat > & aKt,
                                                                  Matrix< DDRMat > & aLt,
                                                            const Matrix< DDRMat > & ad2NdXi2,
                                                            const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------
        /**
         * evaluates matrices that are needed for the second derivative
         * in space, 2D version. It is
         *
         * \f[
         *      \mathbf{L}^T \, \mathbf{\frac{\partian^2 N}{\partial x^2}}
         *      = \mathbf{\frac{\partian^2 N}{\partial \xi^2}}
         *      - K^T \, mathbf{\frac{\partian N}{\partial x}}
         * \f]
         *
         * @param[ in ]  aJt          transposed of geometry Jacobian
         * @param[ out ] aKt          transposed help matrix K
         * @param[ out ] aLt          transposed help matrix L
         * @param[ in ]  adNdXi       first derivatives in parameter space
         * @param[ in ]  ad2NdX2i     second derivatives in parameter space
         *
         */
        static void eval_matrices_for_second_derivative_2d( const Matrix< DDRMat > & aJt,
                                                                  Matrix< DDRMat > & aKt,
                                                                  Matrix< DDRMat > & aLt,
                                                            const Matrix< DDRMat > & ad2NdXi2,
                                                            const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------
        /**
         * evaluates matrices that are needed for the second derivative
         * in space, 3D version. It is
         *
         * \f[
         *      \mathbf{L}^T \, \mathbf{\frac{\partian^2 N}{\partial x^2}}
         *      = \mathbf{\frac{\partian^2 N}{\partial \xi^2}}
         *      - K^T \, mathbf{\frac{\partian N}{\partial x}}
         * \f]
         *
         * @param[ in ]  aJt          transposed of geometry Jacobian
         * @param[ out ] aKt          transposed help matrix K
         * @param[ out ] aLt          transposed help matrix L
         * @param[ in ]  adNdXi       first derivatives in parameter space
         * @param[ in ]  ad2NdX2i     second derivatives in parameter space
         *
         */
        static void eval_matrices_for_second_derivative_3d( const Matrix< DDRMat > & aJt,
                                                                  Matrix< DDRMat > & aKt,
                                                                  Matrix< DDRMat > & aLt,
                                                            const Matrix< DDRMat > & ad2NdXi2,
                                                            const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------
        /**
         * get the geometry type of a side
         */
        void get_auto_side_geometry_type();

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */
