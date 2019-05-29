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
        // and matrix of time coefficients tHat
        Matrix < DDRMat > mXHat;
        Matrix < DDRMat > mTHat;

        // matrix of space param coefficients xiHat in the interpolation param space
        // and matrix of time param coefficients tauHat in the interpolation param space
        Matrix < DDRMat > mXiHat;
        Matrix < DDRMat > mTauHat;

        // element geometry type
        mtk::Geometry_Type mGeometryType;
        mtk::Geometry_Type mTimeGeometryType;

        // boolean true if side interpolation
        bool mSpaceSideset = false;

        //fixme to remove when GI on the side
        // side geometry type
        mtk::Geometry_Type mSideGeometryType;

        // pointer to side space interpolation function object
        Interpolation_Function_Base * mSideSpaceInterpolation = nullptr;

        // number of space bases, number of physical and parametric dimensions for the side
        uint mSideNumSpaceBases;
        uint mSideNumSpaceDim;
        uint mSideNumSpaceParamDim;

        // matrix of space coefficients for the side in the physical space
        // and matrix of time coefficients for the side in the physical space
        Matrix< DDRMat > mSideXHat;
        Matrix< DDRMat > mSideTHat;

        // matrix of space coefficients for the side in the interpolation parametric space
        // matrix of time coefficients for the side in the interpolation parametric space
        Matrix< DDRMat > mSideXiHat;
        Matrix< DDRMat > mSideTauHat;
        //end fixme to remove when GI on the side

        // Vertices ordinals for each face of the parent element
        moris::Cell< moris::Cell< moris::moris_index > > mVerticesPerFace;

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
              // check that mXHat is set
              MORIS_ASSERT( mXHat.numel()>0, "Geometry_Interpolator::get_space_coeff - mXHat is not set." );

              return mXHat;
          }

//------------------------------------------------------------------------------
          /**
           * get the time coefficients of the geometry field tHat
           */
           Matrix< DDRMat > get_time_coeff() const
           {
               // check that mTHat is set
               MORIS_ASSERT( mTHat.numel()>0, "Geometry_Interpolator::get_time_coeff - mTHat is not set." );

               return mTHat;
           }

//------------------------------------------------------------------------------
        /**
         * set the space param coefficients of the geometry field xiHat
         * @param[ in ] space coefficients
         */
           // fixme param coords form integration mesh
        void set_param_coeff();

//------------------------------------------------------------------------------
        /**
         * set the space param coefficients of the geometry field xiHat
         * @param[ in ] space coefficients
         */
        void set_space_param_coeff( const Matrix< DDRMat > & aXiHat );

//------------------------------------------------------------------------------
        /**
         * set the time param coefficients of the geometry field tHat
         * @param[ in ] time coefficients
         */
        void set_time_param_coeff( const Matrix< DDRMat > & aTauHat );

//------------------------------------------------------------------------------
        /**
         * set the space param coefficients for the side
         * @param[ in ] side space coefficients
         */
        void set_side_space_param_coeff( const Matrix< DDRMat > & aSideXiHat );

//------------------------------------------------------------------------------
        /**
         * set the time param coefficients for the side
         * @param[ in ] side time coefficients
         */
        void set_side_time_param_coeff( const Matrix< DDRMat > & aSideTauHat );

//------------------------------------------------------------------------------
         /**
          * get the space parametric coefficients of the geometry field xiHat
          */
          Matrix< DDRMat > get_space_param_coeff() const
          {
              // check that mXiHat is set
              MORIS_ASSERT( mXiHat.numel()>0, "Geometry_Interpolator::get_space_param_coeff - mXiHat is not set." );

              return mXiHat;
          }

//------------------------------------------------------------------------------
          /**
           * get the time parametric coefficients of the geometry field tauHat
           */
           Matrix< DDRMat > get_time_param_coeff() const
           {
               // check that mTauHat is set
               MORIS_ASSERT( mTauHat.numel()>0, "Geometry_Interpolator::get_time_param_coeff - mTauHat is not set." );

               return mTauHat;
           }

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
        Matrix< DDRMat > extract_space_side_space_param_coeff( moris_index aSpaceOrdinal );

        Matrix< DDRMat > extract_space_param_coeff();

//------------------------------------------------------------------------------
        /**
         * get the parametric coordinates of a space side
         */
        void build_time_side_time_param_coeff( moris_index aTimeOrdinal )
        {
            mSideTauHat = {{ mTauHat( aTimeOrdinal ) }};
        }

//------------------------------------------------------------------------------
        /**
         * build the time physical coordinates for a time side
         * @param[ in ] a time side ordinal
         */
        void build_time_side_time_phys_coeff( moris_index aTimeOrdinal )
        {
            mSideTHat = {{ mTHat( aTimeOrdinal ) }};
        }

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
         */
         real time_surf_det_J( const Matrix< DDRMat > & aSideParamPoint );

//------------------------------------------------------------------------------
         /**
          * get the parametric coordinates of a point
          * in the side parametric space
          * in the parent parametric space
          */
         Matrix < DDRMat > time_surf_val( const Matrix< DDRMat > & aSideParamPoint );

//------------------------------------------------------------------------------
         /**
          * evaluates the normal to the side
          * in the case of a space side interpolation
          * at given space and time evaluation point
          * @param[ in ]  aSideParamPoint evaluation point on the face
          */
         Matrix< DDRMat > get_normal( const Matrix< DDRMat > & aSideParamPoint );

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
        /**
         * map an integration point from local param coords to global param coords
         * @param[ out ] aLocalParamPoint param coords in local parametric space
         * @param[ in ]  aParamPoint      param coords in global parametric space
         */
        Matrix< DDRMat > map_integration_point( const Matrix< DDRMat > & aLocalParamPoint );

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

////------------------------------------------------------------------------------
//         /**
//          * evaluates the determinant of the Jacobian mapping
//          * in the case of a space side interpolation
//          * at given space and time evaluation point
//          * @param[ out ] aTimeSurfDetJ   determinant of the Jacobian
//          * @param[ out ] aNormal         normal to the face
//          * @param[ in ]  aSideParamPoint evaluation point on the face
//          * @param[ in ]  aSpaceOrdinal   a space face ordinal
//          */
//         void surf_det_J(       real             & aSurfDetJ,
//                                Matrix< DDRMat > & aNormal,
//                          const Matrix< DDRMat > & aSideParamPoint,
//                                moris_index        aSpaceOrdinal );
//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */
