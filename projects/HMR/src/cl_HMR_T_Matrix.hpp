/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_T_MATRIX_HPP_
#define SRC_HMR_CL_HMR_T_MATRIX_HPP_

#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Cell.hpp" //CON/src

namespace moris
{
    namespace hmr
    {
        /**
         *  \brief a T-Matrix Generator
         */
//-------------------------------------------------------------------------------

        class T_Matrix
        {
            protected:
            //! ref to settings container
            const Parameters   * mParameters;

            //! ref to B-Spline Mesh
            BSpline_Mesh_Base  * mBSplineMesh;

            //! ref to Lagrange Mesh
            Lagrange_Mesh_Base * mLagrangeMesh;

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

            Cell< Cell< Matrix< DDRMat > > >mChild_2;

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

            //! container for gauss points in one direction
            //Matrix< DDRMat > mGaussPoints;

            //! container for gauss weights in one direction
            //Matrix< DDRMat > mGaussWeights;

            //! mass matrix for L2-projection
            // Matrix< DDRMat > mBSplineMass;
            //Matrix< DDRMat > mLagrangeMass;

            //! pointer to T-Matrix calculation function
            //! points to either calculate_untruncated_t_matrix
            //! or calculate_truncated_t_matrix
            void ( T_Matrix:: * mTMatrixFunction )( const luint            & aMemoryIndex,
                                                           Matrix< DDRMat > & aTMatrixTransposed,
                                                           Cell< Basis* >   & aDOFs );

            //! pointer to function for geometry interpolation
            void ( * mEvalNGeo )( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN );

            //! pointer to corner node function
            void ( * mGetCorners )(  const uint & aChildindex, Matrix< DDRMat > & aXi );

            //! pointer to shape function
            void ( T_Matrix :: * mEvalN )( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN ) const;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            // constructor
            T_Matrix( const Parameters         * aParameters,
                            BSpline_Mesh_Base  * aBSplineMesh,
                            Lagrange_Mesh_Base * aLagrangeMesh );

            T_Matrix( const Parameters         * aParameters,
                            Lagrange_Mesh_Base * aLagrangeMesh );

//-------------------------------------------------------------------------------

            // destructor
            virtual ~T_Matrix();

//-------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_lagrange_matrix()
            {
                return mTMatrixLagrange;
            }

//-------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_refinement_matrix( const uint & aChildIndex ) const
            {
                return mLagrangeRefinementMatrix( aChildIndex );
            }

//-------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_change_order_matrix( const uint & aOrder ) const
            {
                return mLagrangeChangeOrderMatrix( aOrder );
            }

//-------------------------------------------------------------------------------

            void calculate_t_matrix( const luint            & aMemoryIndex,
                                           Matrix< DDRMat > & aTMatrixTransposed,
                                           Cell< Basis* >   & aDOFs );

//-------------------------------------------------------------------------------

            void calculate_untruncated_t_matrix( const luint            & aMemoryIndex,
                                                       Matrix< DDRMat > & aTMatrixTransposed,
                                                       Cell< Basis* >   & aDOFs );

//-------------------------------------------------------------------------------

            void calculate_truncated_t_matrix( const luint            & aMemoryIndex,
                                                     Matrix< DDRMat > & aTMatrixTransposed,
                                                     Cell< Basis* >   & aDOFs );

//-------------------------------------------------------------------------------

            const Matrix< DDRMat > & get_child_matrix( const uint aChildIndex ) const
            {
                return mChild( aChildIndex );
            }

//-------------------------------------------------------------------------------
            virtual void evaluate( const uint aBSplineMeshIndex,
                           const bool aBool = true);

//-------------------------------------------------------------------------------

            void evaluate_trivial( const uint aBSplineMeshIndex,
                                   const bool aBool );

//-------------------------------------------------------------------------------
       private:
//-------------------------------------------------------------------------------

            /**
             * determins the sorting order of the nodes
             */
            void init_basis_index();

//-------------------------------------------------------------------------------

            /**
             * initializes the container for the unity matrix
             */
            void init_unity_matrix();

//-------------------------------------------------------------------------------
            /**
             * pre-calculates the child relation matrices
             */
            void init_child_matrices();

//-------------------------------------------------------------------------------

            void child_multiplication();

//-------------------------------------------------------------------------------

            const Matrix< DDRMat> & get_child_matrix_1( const Cell< uint > & aChildIndices );

//------------------------------------------------------------------------------

            void init_truncation_weights();

//------------------------------------------------------------------------------

            void init_lagrange_parameter_coordinates();

//------------------------------------------------------------------------------

            /**
             * calculates the matrix that converts B-Spline DOFs per element to Lagrange DOFs.
             */
            void init_lagrange_matrix();

//------------------------------------------------------------------------------

            void init_lagrange_refinement_matrices();

//-------------------------------------------------------------------------------

            void init_lagrange_change_order_matrices();

//-------------------------------------------------------------------------------

            /**
             * initializes interpolation coefficiens for Lagrange interpolation
             */
            void init_lagrange_coefficients();

//-------------------------------------------------------------------------------

            /**
             * initializes gauss points
             */
            //void
            //init_gauss_points();

 //-------------------------------------------------------------------------------

            /**
             * calculates the B-Spline mass matrix
             */
            //void
            //init_mass_matrices();

//-------------------------------------------------------------------------------

            Background_Element_Base * create_background_element();

//------------------------------------------------------------------------------

            /**
             * 1D shape function
             */
            real b_spline_shape_1d( const uint & aOrder,
                                    const uint & aK,
                                    const real & aXi ) const;

//------------------------------------------------------------------------------

            /**
             * 2D shape function
             */
            void b_spline_shape( const real             & aXi,
                                 const real             & aEta,
                                       Matrix< DDRMat > & aN ) const;

//------------------------------------------------------------------------------

            /**
             * 3D shape function
             */
            void b_spline_shape( const real             & aXi,
                                 const real             & aEta,
                                 const real             & aZeta,
                                       Matrix< DDRMat > & aN ) const;

//------------------------------------------------------------------------------

            /**
             * 1D shape function
             */
            real lagrange_shape_1d( const uint & aBasisNumber,
                                    const real & aXi ) const;

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

//------------------------------------------------------------------------------

            /**
             * calculates the Legendre polynomials
             *
             * @param[ in ] number of polynomial - 1
             * @param[ in ] x value of supporting point
             * @param[ in ] value of polynomial
             * @param[ in ] first derivative of polynomial
             */
            void legendre( const uint         & aIndex,
                           const long double  & aX,
                                 long double  & aP,
                                 long double  & adPdX ) const;

//------------------------------------------------------------------------------

            /**
             * returns the corner nodes of a child and dimension
             */
            static void get_child_corner_nodes_2d( const uint & aChildIndex, Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------

            /**
             * returns the corner nodes of a child and dimension
             */
            static void get_child_corner_nodes_3d(  const uint & aChildIndex, Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------

            /**
             * quam4 shape function
             */
            static void N_quad4( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN );

//------------------------------------------------------------------------------

            /**
             * quam4 shape function
             */
            static void N_hex8( const Matrix< DDRMat > & aXi, Matrix< DDRMat > & aN );

//------------------------------------------------------------------------------

            Matrix< DDRMat > get_supporting_points( const uint aDimension, const uint aOrder );

        };

//------------------------------------------------------------------------------
    }
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_T_MATRIX_HPP_ */

