/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Factory.hpp
 *
 */

#pragma once

#include "cl_HMR_Background_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"       //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"      //HMR/src
#include "cl_HMR_T_Matrix_Advanced.hpp"
#include "cl_HMR_Parameters.hpp"              //HMR/src
#include "moris_typedefs.hpp"                       //COR/src
#include "cl_Matrix.hpp"                      //LINALG/src

namespace moris::hmr
{
    // Forward declare T-matrix template
    template< uint N >
    class T_Matrix;

    // ----------------------------------------------------------------------------
    /**
     * \brief factory class that generates pointers to templated meshes
     */
    class Factory
    {
        // ----------------------------------------------------------------------------

      private:
        // Stored parameters for building classes
        const Parameters* mParameters;

        // ----------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------

        /**
         * Constructor which takes in the HMR parameters needed to construct classes
         *
         * @param aParameters HMR parameters to use
         */
        explicit Factory( const Parameters* aParameters );

        // ----------------------------------------------------------------------------

        /**
         * creates a background mesh depending on the number of dimensions set
         *
         * @return Background_Mesh_Base*   pointer to new background mesh
         */
        Background_Mesh_Base* create_background_mesh();

        // ----------------------------------------------------------------------------

        /**
         * @brief creates a Lagrange mesh depending on the number of dimensions set
         *
         * @param[in] aBackgroundMesh       pointer to background mesh
         * @param[in] aBSplineMeshes
         * @param[in] aActivationPattern
         * @param[in] aPolynomialDegree     degree of Lagrange mesh
         * @param[in] aMeshIndex
         *
         * @return Lagrange_Mesh_Base* pointer to new Lagrange mesh
         */
        Lagrange_Mesh_Base*
        create_lagrange_mesh(
                Background_Mesh_Base*      aBackgroundMesh,
                Vector< BSpline_Mesh_Base* > aBSplineMeshes,
                uint                       aActivationPattern,
                luint                      aPolynomialDegree,
                uint                       aMeshIndex = MORIS_UINT_MAX );

        // ----------------------------------------------------------------------------

        /**
         * @brief creates a Lagrange mesh depending on the number of dimensions set
         *
         * @param[in] aBackgroundMesh       pointer to background mesh
         * @param[in] aActivationPattern    pattern (i.e. grid) used for building the mesh
         * @param[in] aPolynomialDegree     degree of Lagrange mesh
         * @param[in] aMeshIndex            index to be assigned to the mesh
         *
         * @return BSpline_Mesh_Base* pointer to new Lagrange mesh
         */
        BSpline_Mesh_Base*
        create_bspline_mesh(
                Background_Mesh_Base* aBackgroundMesh,
                uint                  aActivationPattern,
                luint                 aPolynomialDegree,
                uint                  aMeshIndex = MORIS_UINT_MAX );

        // ----------------------------------------------------------------------------

        /**
         * Creates a T-matrix specified by the given parameters
         *
         * @param aLagrangeMesh Lagrange mesh pointer
         * @param aBSplineMesh B-spline mesh pointer
         * @param aLagrangeMeshFine Pointer to finer Lagrange mesh, only needed for creating advanced T-matrices
         * @return Base T-matrix pointer
         */
        T_Matrix_Base* create_t_matrix(
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh        = nullptr,
                Lagrange_Mesh_Base* aLagrangeMeshCoarse = nullptr );

        // ----------------------------------------------------------------------------

        /**
         * Creates a T-matrix for a specific dimension, if desired and known at compile-time
         *
         * @param N Spatial dimension
         * @param aLagrangeMesh Lagrange mesh pointer
         * @param aBSplineMesh B-spline mesh pointer
         * @param aLagrangeMeshFine Pointer to finer Lagrange mesh, only needed for creating advanced T-matrices
         * @return T-matrix pointer
         */
        template< uint N >
        T_Matrix< N >*
        create_t_matrix(
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh      = nullptr,
                Lagrange_Mesh_Base* aLagrangeMeshFine = nullptr )
        {
            // Use Advanced T-matrices
            if ( mParameters->use_advanced_t_matrices() and aLagrangeMeshFine )
            {
                return new T_Matrix_Advanced< N >( aLagrangeMeshFine, aBSplineMesh, aLagrangeMesh, mParameters->truncate_bsplines() );
            }

            // Use regular T-matrices
            else
            {
                return new T_Matrix< N >( aLagrangeMesh, aBSplineMesh, mParameters->truncate_bsplines() );
            }
        }

        // ----------------------------------------------------------------------------

    };    // end class: Factory

}    // namespace moris::hmr
