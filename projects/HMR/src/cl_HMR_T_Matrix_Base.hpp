/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_Base.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_T_MATRIX_BASE_HPP_
#define SRC_HMR_CL_HMR_T_MATRIX_BASE_HPP_

#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Cell.hpp" //CNT/src

namespace moris::hmr
{
    /**
     * Base T-matrix class
     */
    class T_Matrix_Base
    {
    public:
        //! ref to settings container
        const Parameters   * mParameters;

        //! ref to Lagrange Mesh
        Lagrange_Mesh_Base * mLagrangeMesh;

        //! ref to B-Spline Mesh
        BSpline_Mesh_Base  * mBSplineMesh;

        //! matrix containing the ijk positions of reference element
        Matrix< DDUMat > mBSplineIJK;

        //! matrix containing the ijk positions of reference element / ijk position of nodes an element
        Matrix< DDUMat > mLagrangeIJK;

        //! ordering scheme for Elements
        Matrix< DDUMat > mBasisIndex;

        // unity matrix
        Matrix< DDRMat > mEye;

        // zero vector
        Matrix< DDRMat > mZero;

        // cell containing child matrices ( transposed )
        Cell< Matrix< DDRMat > > mChild;
        Cell< Cell< Matrix< DDRMat > > >mChildMultiplied;

        //! parameter coordinates for lagrange element / natural coordinates
        Matrix< DDRMat > mLagrangeParam;

        //! T-Matrix for B-Spline to Lagrange conversion
        Matrix< DDRMat > mTMatrixLagrange;

        //! weights needed for truncation
        Matrix< DDRMat > mTruncationWeights;

        //! Lagrange coefficients for interpolation
        Matrix< DDRMat > mLagrangeCoefficients;

        //! order of B-Spline Mesh
        uint        mBSplineOrder;

        //! order of Lagrange Mesh
        uint        mLagrangeOrder;

        //! number of nodes per Lagrange element
        uint        mNumberOfNodes;

        //! matrices for refining Lagrange node values
        Cell< Matrix< DDRMat > > mLagrangeRefinementMatrix;

        //! matrices for changing the order of a Lagrange mesh
        Cell< Matrix< DDRMat > > mLagrangeChangeOrderMatrix;

        //! pointer to T-Matrix calculation function
        //! points to either calculate_untruncated_t_matrix
        //! or calculate_truncated_t_matrix
        void ( T_Matrix_Base:: *mTMatrixFunction )( luint aMemoryIndex,
                                                       Matrix< DDRMat > & aTMatrixTransposed,
                                                       Cell< Basis* >   & aDOFs );

        //! pointer to function for geometry interpolation
        void ( * mEvalNGeo )( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN );

        //! pointer to corner node function
        void ( * mGetCorners )(  uint aChildIndex, Matrix< DDRMat > & aXi );

    public:

        /**
         * Constructor initializing Lagrange coefficients
         *
         * @param aParameters HMR Parameters
         * @param aLagrangeMesh Lagrange mesh pointer
         * @param aBSplineMesh B-spline Mesh pointer
         */
        T_Matrix_Base(
                const Parameters*   aParameters,
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh = nullptr );

        // destructor
        virtual ~T_Matrix_Base();

        //-------------------------------------------------------------------------------

        const Matrix< DDRMat > & get_lagrange_matrix()
        {
            return mTMatrixLagrange;
        }

        //-------------------------------------------------------------------------------

        const Matrix< DDRMat > & get_refinement_matrix( uint aChildIndex ) const
        {
            return mLagrangeRefinementMatrix( aChildIndex );
        }

        //-------------------------------------------------------------------------------

        const Matrix< DDRMat > & get_change_order_matrix( uint aOrder ) const
        {
            return mLagrangeChangeOrderMatrix( aOrder );
        }

        //-------------------------------------------------------------------------------

        void calculate_t_matrix(
                luint             aMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

        //-------------------------------------------------------------------------------

        void calculate_untruncated_t_matrix(
                luint             aMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

        //-------------------------------------------------------------------------------

        void calculate_truncated_t_matrix(
                luint             aMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

        virtual void evaluate( uint aBSplineMeshIndex,
                               bool aBool = true);

        void evaluate_trivial( uint aBSplineMeshIndex,
                               bool aBool );

    private:

        /**
         * Initializes lagrange coefficients for Lagrange interpolation
         */
        void init_lagrange_coefficients();

        const Matrix< DDRMat> & get_child_matrix( const Cell< uint > & aChildIndices );


        //------------------------------------------------------------------------------

        /**
         * 1D shape function
         */
        static real b_spline_shape_1d( uint aOrder,
                                uint aK,
                                real aXi );

        //------------------------------------------------------------------------------

        /**
         * 1D shape function
         */
        real lagrange_shape_1d( uint aBasisNumber,
                                real aXi ) const;

        //------------------------------------------------------------------------------

        /**
         * 2D shape function
         */
        void lagrange_shape_2d( const Matrix< DDRMat > & aXi,
                                      Matrix< DDRMat > & aN ) const;

        //------------------------------------------------------------------------------

        /**
         * 3D shape function
         */
        void lagrange_shape_3d( const Matrix< DDRMat > & aXi,
                                      Matrix< DDRMat > & aN ) const;
    };
}

#endif /* SRC_HMR_CL_HMR_T_MATRIX_BASE_HPP_ */

