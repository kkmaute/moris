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
    protected:

        //! ref to Lagrange Mesh
        Lagrange_Mesh_Base* mLagrangeMesh;

        //! ref to B-Spline Mesh
        BSpline_Mesh_Base* mBSplineMesh;

        //! Whether or not to truncate B-splines
        const bool mTruncate;

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

        //! Modified parameter coordinates for lagrange element / natural coordinates
        Matrix< DDRMat > mLagrangeParamModified;

        //! T-Matrix for B-Spline to Lagrange conversion
        Matrix< DDRMat > mTMatrixLagrange;

        //! Modified T-Matrix for B-Spline to Lagrange conversion
        Matrix< DDRMat > mTMatrixLagrangeModified;

        //! weights needed for truncation
        Matrix< DDRMat > mTruncationWeights;

        //! Lagrange coefficients for interpolation
        Matrix< DDRMat > mLagrangeCoefficients;

        //! order of Lagrange Mesh
        uint        mLagrangeOrder;

        //! number of nodes per Lagrange element
        uint        mNumberOfNodes;

        //! matrices for refining Lagrange node values
        Cell< Matrix< DDRMat > > mLagrangeRefinementMatrix;

        //! matrices for changing the order of a Lagrange mesh
        Cell< Matrix< DDRMat > > mLagrangeChangeOrderMatrix;

    public:

        /**
         * Constructor initializing Lagrange coefficients
         *
         * @param aLagrangeMesh Lagrange mesh pointer
         * @param aBSplineMesh B-spline Mesh pointer
         * @param aTruncate Whether or not to truncate B-splines
         */
        explicit T_Matrix_Base(
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh = nullptr,
                bool                aTruncate = true );

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

        /**
         * Evaluates an extended T-matrix
         *
         * @param aBsplineElement B-spline element
         * @param aLagrangeElement Lagrange element
         * @param aBsplineBasis B-spline basis
         * @param aWeights Weights
         */
        void evaluate_extended_t_matrix(
                Element*                                    aBsplineElement,
                Element*                                    aLagrangeElement,
                moris::Cell< moris::Cell< mtk::Vertex* > >& aBsplineBasis,
                moris::Cell< Matrix< DDRMat > >&            aWeights );

        virtual void evaluate( uint aBSplineMeshIndex,
                               bool aBool = true);

        void evaluate_trivial( uint aBSplineMeshIndex,
                               bool aBool );

        /**
         * Calculates the truncated or untruncated T-matrix for a B-spline element, based on truncation parameter.
         *
         * @param aElementMemoryIndex Memory index of the B-spline element for T-matrix computation
         * @param aTMatrixTransposed Transposed T-matrix
         * @param aDOFs Active bases on the element
         */
        void calculate_t_matrix(
                luint             aElementMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

    private:

        /**
         * Calculates the untruncated T-matrix for a B-spline element.
         *
         * @param aElementMemoryIndex Memory index of the B-spline element for T-matrix computation
         * @param aTMatrixTransposed Transposed T-matrix
         * @param aDOFs Active bases on the element
         */
        void calculate_untruncated_t_matrix(
                luint             aElementMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

        /**
         * Calculates the truncated T-matrix for a B-spline element.
         *
         * @param aElementMemoryIndex Memory index of the B-spline element for T-matrix computation
         * @param aTMatrixTransposed Transposed T-matrix
         * @param aDOFs Active bases on the element
         */
        void calculate_truncated_t_matrix(
                luint             aElementMemoryIndex,
                Matrix< DDRMat >& aTMatrixTransposed,
                Cell< Basis* >&   aDOFs );

        /**
         * Initializes lagrange coefficients for Lagrange interpolation
         */
        void init_lagrange_coefficients();

        /**
         * Initializes the modified lagrange parameter coordinates
         *
         * @param aLagrangeElement Lagrange element
         * @param aBSplineElement B-spline element
         */
        virtual void init_modified_lagrange_parameter_coordinates(
                Element* aLagrangeElement,
                Element* aBSplineElement ) = 0;

        const Matrix< DDRMat> & get_child_matrix( const Cell< uint > & aChildIndices );

        /**
         * Recompute the Lagrange matrix for extended T-matrices
         */
        virtual void recompute_lagrange_matrix() = 0;

    };
}

#endif /* SRC_HMR_CL_HMR_T_MATRIX_BASE_HPP_ */
