/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Mesh.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_

#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp"
#include "cl_HMR_BSpline_Mesh_Base.hpp"
#include "cl_HMR_Lagrange_Edge.hpp"
#include "cl_HMR_Lagrange_Edge2.hpp"
#include "cl_HMR_Lagrange_Edge3.hpp"
#include "cl_HMR_Lagrange_Edge4.hpp"
#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_HMR_Lagrange_Facet.hpp"
#include "cl_HMR_Lagrange_Facet_Line2.hpp"
#include "cl_HMR_Lagrange_Facet_Line3.hpp"
#include "cl_HMR_Lagrange_Facet_Line4.hpp"
#include "cl_HMR_Lagrange_Facet_Quad16.hpp"
#include "cl_HMR_Lagrange_Facet_Quad4.hpp"
#include "cl_HMR_Lagrange_Facet_Quad9.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR_T_Matrix.hpp"
#include "cl_HMR_T_Matrix_Advanced.hpp"
#include "HMR_Globals.hpp"
#include "moris_typedefs.hpp"
#include "cl_Stopwatch.hpp"
// #include "cl_Map.hpp"
#include "cl_Tracer.hpp"
#include "fn_stringify.hpp"

namespace moris::hmr
{
    /**
     * Lagrange mesh class
     *
     * @param N Number of dimensions
     * @param P Polynomial order
     */
    template< uint N, uint P >
    class Lagrange_Mesh : public Lagrange_Mesh_Base
    {
        // ----------------------------------------------------------------------------

        //! Lookup table containing offset for node IDs
        luint mNodeLevelOffset[ gMaxNumberOfLevels ];

        //! Lookup table containing number of elements per dimension for each level
        luint mNumberOfElementsPerDimensionIncludingAura[ gMaxNumberOfLevels ][ N ];

        //! Lookup table for node IDs
        luint mMySubdomainOffset[ gMaxNumberOfLevels ][ N ];

        //! calculation object that calculates the T-Matrices
        Vector< T_Matrix< N >* > mTMatrix;

        // ----------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------

        /**
         * Constructor for Lagrange Mesh
         *
         * @param[in] aParameters       ref to container of user defined settings
         * @param[in] aBackgroundMesh pointer to background mesh
         *
         */
        Lagrange_Mesh(
                const Parameters*           aParameters,
                Background_Mesh_Base*       aBackgroundMesh,
                Vector< BSpline_Mesh_Base* >& aBSplineMeshes,
                uint                        aActivationPattern,
                uint                        aMeshIndex )
                : Lagrange_Mesh_Base(
                        aParameters,
                        aBackgroundMesh,
                        aBSplineMeshes,
                        P,
                        aActivationPattern )
        {
            // trace this operation
            Tracer tTracer( "HMR", "Lagrange Mesh #" + std::to_string( aMeshIndex ), "Create" );

            // collect the B-spline mesh indices which this Lagrange is associated with
            uint                tNumBspMeshes = aBSplineMeshes.size();
            Vector< moris_index > tBsplineMeshIndices( tNumBspMeshes );
            for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get access to the B-spline mesh
                BSpline_Mesh_Base* tBsplineMeshPtr = aBSplineMeshes( iBspMesh );

                // copy its index, but only if it isn't a nullptr
                if ( tBsplineMeshPtr )
                {
                    tBsplineMeshIndices( iBspMesh ) = tBsplineMeshPtr->get_index();
                }
                else    // B-spline mesh not constructed as the Lagrange mesh instead uses its own basis functions (probably)
                {
                    tBsplineMeshIndices( iBspMesh ) = -1;
                }
            }

            // report based on what this lagrange mesh is created on
            MORIS_LOG_INFO(
                    "Creating Lagrange mesh index #%i with polynomial order p=%i on pattern #%i, associated with B-spline meshes %s.",
                    aMeshIndex,
                    mOrder,
                    aActivationPattern,
                    ios::stringify_cell( tBsplineMeshIndices.data() ).c_str() );

            // assign and store the mesh index
            this->set_index( aMeshIndex );

            // ask background mesh for number of elements per ijk-direction
            this->get_number_of_elements_per_dimension();

            // calculate lookup table mNodeLevelOffset
            this->calculate_level_offset();

            // find out coordinate of first point on proc subdomain
            this->calculate_subdomain_offset();

            // calculate any value that can change after refinement
            this->update_mesh();

            // Initialize T-matrices
            this->init_t_matrices();
        }

        // ----------------------------------------------------------------------------

        /**
         * Default destructor.
         */
        ~Lagrange_Mesh() override
        {
            this->delete_t_matrices();
            this->delete_pointers();
            this->delete_facets();

            if ( N == 3 )
            {
                this->delete_edges();
            }

            this->delete_t_matrix_lagrange_mesh();
        }

        // ----------------------------------------------------------------------------

        /**
         * Creates a Lagrange element and links it to corresponding element
         * on background mesh.
         *
         * @param[in] aElement  pointer to element on background mesh
         *
         * @return Element*  new Lagrange element
         */
        Element* create_element( Background_Element_Base* aElement ) override;

        // ----------------------------------------------------------------------------

      protected:
        // ----------------------------------------------------------------------------

        Facet* create_facet( Background_Facet* aFacet ) override;

        // ----------------------------------------------------------------------------

        Edge* create_edge( Background_Edge* aEdge ) override;

        // ----------------------------------------------------------------------------

      private:
        // ----------------------------------------------------------------------------
        /**
         * Initializes T-matrices
         */
        void
        init_t_matrices()
        {
            // report on this operation
            MORIS_LOG_INFO( "Initializing T-matrices on HMR Lagrange mesh" );

            // Resize for T-matrices
            mTMatrix.resize( mNumBSplineMeshes, nullptr );
            mLagrangeMeshForTMatrix.resize( mNumBSplineMeshes, nullptr );

            // initialize T-matrices for every B-spline mesh associated with the current Lagrange mesh
            for ( uint iBspMesh = 0; iBspMesh < mNumBSplineMeshes; iBspMesh++ )
            {
                // create factory object
                Factory tFactory( mParameters );

                // Get B-spline mesh and order
                BSpline_Mesh_Base* tMesh = mBSplineMeshes( iBspMesh );

                // Check if B-spline mesh exists
                if ( tMesh )
                {
                    // Check if Lagrange order is less than B-spline order for advanced T-matrices
                    uint tBSplineOrder = tMesh->get_min_order();
                    if ( P < tBSplineOrder and mParameters->use_advanced_t_matrices() )
                    {
                        // get the activation pattern associated with the current mesh
                        uint tActivationPattern = this->get_activation_pattern();

                        // create Lagrange mesh object
                        mLagrangeMeshForTMatrix( iBspMesh ) = tFactory.create_lagrange_mesh(
                                mBackgroundMesh,
                                mBSplineMeshes,
                                tActivationPattern,
                                tBSplineOrder );
                    }
                }

                // Create T-matrix
                // Note: the last input is a pointer to a finer Lagrange Mesh if advanced T-matrix scheme is used, otherwise it's just a nullptr
                mTMatrix( iBspMesh ) = tFactory.create_t_matrix< N >(
                        this,
                        tMesh,
                        mLagrangeMeshForTMatrix( iBspMesh ) );
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique node ID (1D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                uint  aLevel,
                luint aI ) override
        {
            if ( aLevel < gMaxNumberOfLevels && N == 1 )
            {
                return aI + mMySubdomainOffset[ aLevel ][ 0 ]
                     + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique node ID (2D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                uint  aLevel,
                luint aI,
                luint aJ ) override
        {
            if ( aLevel < gMaxNumberOfLevels && N == 2 )
            {
                return aI + mMySubdomainOffset[ aLevel ][ 0 ]
                     + ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] )
                               * ( P * mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 0 ] + 1 )
                     + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique node ID (3D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of node
         * @param[in]  aI        proc local i-position of node
         * @param[in]  aJ        proc local j-position of node
         * @param[in]  aK        proc local k-position of node
         * @return uint          domain wide unique ID
         */
        luint
        calculate_node_id(
                uint  aLevel,
                luint aI,
                luint aJ,
                luint aK ) override
        {
            if ( aLevel < gMaxNumberOfLevels && N == 3 )
            {
                return aI + mMySubdomainOffset[ aLevel ][ 0 ]
                     + ( P * mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 0 ] + 1 )
                               * ( ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] )
                                       + ( aK + mMySubdomainOffset[ aLevel ][ 2 ] )
                                                 * ( P * mNumberOfElementsPerDimensionIncludingAura[ aLevel ][ 1 ] + 1 ) )
                     + mNodeLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Internal function. Asks the background mesh for number of elements
         * per direction, stores result in  mAuraNumberOfElementsPerDimension
         *
         * @return void
         *
         */
        void
        get_number_of_elements_per_dimension()
        {
            // get elements per level from background mesh
            Matrix< DDLUMat > tMat =
                    mBackgroundMesh->get_number_of_elements_per_direction();

            // convert matrix to fixed size array
            for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < N; ++k )
                {
                    mNumberOfElementsPerDimensionIncludingAura[ l ][ k ] = tMat( k, l );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         *  Private function, creates the mNodeLevelOffset lookup table.
         *
         *  @return void
         */
        void
        calculate_level_offset()
        {
            // calculate node level offset
            mNodeLevelOffset[ 0 ] = 0;

            for ( uint l = 1; l < gMaxNumberOfLevels; ++l )
            {
                // calculate number of nodes on this level
                luint tNumberOfNodes = 1;
                for ( uint k = 0; k < N; ++k )
                {
                    tNumberOfNodes *= P * mNumberOfElementsPerDimensionIncludingAura[ l - 1 ][ k ] + 1;
                }

                // add number of nodes to offset table
                mNodeLevelOffset[ l ] = mNodeLevelOffset[ l - 1 ] + tNumberOfNodes;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Private function calculates the mMySubdomainOffset lookup table
         *
         * @return void
         */
        void
        calculate_subdomain_offset()
        {
            Matrix< DDLUMat > tIJK = mBackgroundMesh->get_subdomain_offset_of_proc();

            for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < N; ++k )
                {
                    mMySubdomainOffset[ l ][ k ] = P * tIJK( k, l );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates XZY coordinates for each node
         *
         * @return void
         */
        void
        calculate_node_coordinates() override
        {
            // get domain dimensions from settings
            const Matrix< DDRMat >& tDomainDimensions = mParameters->get_domain_dimensions();

            // get number of elements on coarsest level from settings
            const Vector< luint >& tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

            // calculate step width
            real tDeltaX[ gMaxNumberOfLevels ][ N ];

            // calculate width for first level
            for ( uint k = 0; k < N; ++k )
            {
                tDeltaX[ 0 ][ k ] = tDomainDimensions( k ) / ( (real)( P * tNumberOfElements( k ) ) );
            }

            // loop over all higher levels
            for ( uint l = 1; l < gMaxNumberOfLevels; ++l )
            {
                for ( uint k = 0; k < N; ++k )
                {
                    tDeltaX[ l ][ k ] = 0.5 * tDeltaX[ l - 1 ][ k ];
                }
            }

            // get domain offset
            Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

            // domain offset
            real tOffset[ N ];

            // get coords from background mesh
            Matrix< DDRMat > tOffsetCoords = mBackgroundMesh->get_domain_offset();

            // un-flatten coords to a normal array
            for ( uint k = 0; k < N; ++k )
            {
                tOffset[ k ] = tOffsetCoords( k );
            }

            // loop over all nodes
            for ( auto tNode : mAllBasisOnProc )
            {
                // get ijk position of node
                const luint* tIJK = tNode->get_ijk();

                // get level of node
                luint tLevel = tNode->get_level();

                // array containing coordinate
                real tXYZ[ N ];

                // loop over all dimensions
                for ( uint k = 0; k < N; ++k )
                {
                    tXYZ[ k ] = ( (real)( tIJK[ k ]
                                          + mMySubdomainOffset[ tLevel ][ k ] ) )
                                      * tDeltaX[ tLevel ][ k ]
                              + tOffset[ k ];
                }

                // write XYZ coordinate into node
                tNode->set_xyz( tXYZ );
            }
        }

        // ----------------------------------------------------------------------------

        void
        delete_t_matrices()
        {
            for ( T_Matrix< N >* tTMatrix : mTMatrix )
            {
                delete tTMatrix;
            }
        }

        // ----------------------------------------------------------------------------

        void
        calculate_t_matrices( bool aBool ) override
        {
            // log & trace this operation
            Tracer tTracer( "HMR", "Lagrange Mesh #" + std::to_string( this->get_index() ), "Compute T-matrices" );

            // evaluate T-matrices for each B-spline mesh associated with this Lagrange mesh
            for ( uint iBspMesh = 0; iBspMesh < mNumBSplineMeshes; iBspMesh++ )
            {
                // access the B-spline mesh
                BSpline_Mesh_Base* tMesh = mBSplineMeshes( iBspMesh );

                if ( tMesh != nullptr )
                {
                    uint tBSplineOrder = tMesh->get_min_order();

                    if ( P < tBSplineOrder and mParameters->use_advanced_t_matrices() )
                    {
                        MORIS_ERROR( mLagrangeMeshForTMatrix( iBspMesh ) != nullptr,
                                "Lagrange_Mesh_Base::calculate_t_matrices(), Higher order Lagrange mesh for T-Matrices does not exist." );

                        mLagrangeMeshForTMatrix( iBspMesh )->update_mesh();
                    }

                    // compute the actual T-matrix
                    mTMatrix( iBspMesh )->evaluate( iBspMesh, aBool );
                }
                else
                {
                    mTMatrix( iBspMesh )->evaluate_trivial( iBspMesh, aBool );
                }
            }

            // delete temporary Lagrange meshes for T-matrix evaluation
            this->delete_t_matrix_lagrange_mesh();

        }    // end function: hmr::Lagrange_Mesh::calculate_t_matrices()

        //------------------------------------------------------------------------------

        void
        get_extended_t_matrix(
                moris_index                       aDiscretizationMeshIndex,
                moris_index                       aBSplineCellIndex,
                Element&                          aLagrangeCell,
                Vector< Vector< mtk::Vertex* > >& aBsplineBasis,
                Vector< Matrix< DDRMat > >&       aWeights ) override
        {
            // get B-Spline pattern of this mesh
            uint tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

            // get Lagrange pattern of this mesh
            uint tLagrangePattern = this->get_activation_pattern();

            // set the activation pattern to that of the B-spline mesh
            mBackgroundMesh->set_activation_pattern( tBSplinePattern );

            // get pointer to b-spline and background elements
            Element* tBsplineElement = mBSplineMeshes( 0 )->get_element_including_aura( aBSplineCellIndex );

            // Evaluate extended T-matrix
            mTMatrix( aDiscretizationMeshIndex )->evaluate_extended_t_matrix( tBsplineElement, &aLagrangeCell, aBsplineBasis, aWeights );

            // set the activation pattern back to the pattern of the Lagrange mesh
            mBackgroundMesh->set_activation_pattern( tLagrangePattern );
        }

        //------------------------------------------------------------------------------

        void
        get_L2_projection_matrix(
                moris_index                             aDiscretizationMeshIndex,
                const Element*                          aRootBSplineCell,
                const Element*                          aExtendedBSplineCell,
                Vector< Vector< const mtk::Vertex* > >& aRootBsplineBasis,
                Vector< const mtk::Vertex* >&           aExtendedBsplineBasis,
                Vector< Matrix< DDRMat > >&             aWeights ) override
        {
            // ask the t-matrix object to compute the weights
            mTMatrix( aDiscretizationMeshIndex )->evaluate_L2_projection( aRootBSplineCell, aExtendedBSplineCell, aRootBsplineBasis, aExtendedBsplineBasis, aWeights );
        }

        //------------------------------------------------------------------------------
    };

    template< uint N, uint P >
    inline Element*
    Lagrange_Mesh< N, P >::create_element(
            Background_Element_Base* aElement )
    {
        MORIS_ERROR( false, "Don't know how to create Lagrange element." );
        return nullptr;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 2, 1 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 2, 4 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 2, 2 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 2, 9 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 2, 3 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 2, 16 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 2, 4 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 2, 25 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 2, 5 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 2, 36 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 3, 1 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 3, 8 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 3, 2 >::create_element(
            Background_Element_Base* aElement )
    {

        return new Lagrange_Element< 3, 27 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 3, 3 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 3, 64 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 3, 4 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 3, 125 >( aElement, mActivationPattern );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Element*
    Lagrange_Mesh< 3, 5 >::create_element(
            Background_Element_Base* aElement )
    {
        return new Lagrange_Element< 3, 216 >( aElement, mActivationPattern );
    }

    // ----------------------------------------------------------------------------

    template< uint N, uint P >
    inline Facet*
    Lagrange_Mesh< N, P >::create_facet(
            Background_Facet* aFacet )
    {
        MORIS_ERROR( false, "Don't know how to create Lagrange facet." );
        return nullptr;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 2, 1 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 2, 2 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 2, 2 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 2, 3 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 2, 3 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 2, 4 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 2, 4 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 2, 5 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 2, 5 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 2, 6 >( this, aFacet );
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 3, 1 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 3, 4 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 3, 2 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 3, 9 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 3, 3 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 3, 16 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 3, 4 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 3, 25 >( this, aFacet );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Facet*
    Lagrange_Mesh< 3, 5 >::create_facet(
            Background_Facet* aFacet )
    {
        return new Lagrange_Facet< 3, 36 >( this, aFacet );
    }

    // ----------------------------------------------------------------------------

    template< uint N, uint P >
    inline Edge*
    Lagrange_Mesh< N, P >::create_edge(
            Background_Edge* aEdge )
    {
        MORIS_ERROR( false, "Don't know how to create Lagrange edge." );
        return nullptr;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Edge*
    Lagrange_Mesh< 3, 1 >::create_edge(
            Background_Edge* aEdge )
    {
        return new Lagrange_Edge< 2 >( this, aEdge );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Edge*
    Lagrange_Mesh< 3, 2 >::create_edge(
            Background_Edge* aEdge )
    {
        return new Lagrange_Edge< 3 >( this, aEdge );
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<>
    inline Edge*
    Lagrange_Mesh< 3, 3 >::create_edge(
            Background_Edge* aEdge )
    {
        return new Lagrange_Edge< 4 >( this, aEdge );
    }

}    // namespace moris::hmr

#endif /* SRC_HMR_CL_HMR_LAGRANGE_MESH_HPP_ */
