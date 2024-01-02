/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh.cpp
 *
 */

#include "cl_HMR_Mesh.hpp"    //HMR/src

#include <string>

#include "cl_HMR.hpp"    //HMR/src
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Background_Element.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_HMR_refinement_transition_locations.hpp"
#include "MTK_Tools.hpp"
#include "fn_sort.hpp"
#include "fn_unique.hpp"

namespace moris::hmr
{
    //-----------------------------------------------------------------------------

    Mesh::Mesh(
            std::shared_ptr< Database > aDatabase,
            uint                        aLagrangeOrder,
            uint                        aLagrangePattern )
    {
        // copy database pointer
        mDatabase = aDatabase;

        bool tMeshFound = false;

        // get number of meshes
        uint tNumberOfMeshes = mDatabase->get_number_of_lagrange_meshes();

        // find correct block
        for ( uint k = 0; k < tNumberOfMeshes; ++k )
        {
            auto tMesh = mDatabase->get_lagrange_mesh_by_index( k );

            // test if mesh uses active pattern
            if ( tMesh->get_activation_pattern() == aLagrangePattern &&    // FIXME
                    tMesh->get_order() == aLagrangeOrder )
            {
                mMesh      = tMesh;
                tMeshFound = true;

                break;
            }
        }

        if ( !tMeshFound )
        {
            tNumberOfMeshes = mDatabase->get_number_of_additional_lagrange_meshes( aLagrangePattern );

            for ( uint k = 0; k < tNumberOfMeshes; ++k )
            {
                auto tMesh = mDatabase->get_additional_lagrange_mesh_by_index( k, aLagrangePattern );

                // test if mesh uses active pattern
                if ( tMesh->get_activation_pattern() == aLagrangePattern &&    // FIXME
                        tMesh->get_order() == aLagrangeOrder )
                {
                    mMesh      = tMesh;
                    tMeshFound = true;
                    break;
                }
            }
        }

        MORIS_ERROR( mMesh != nullptr, "Mesh::Mesh(), Mesh not found" );

        if ( mDatabase->is_finalized() )
        {
            // setup_glb_to_local_maps();
        }

        // MORIS_ERROR( mMesh != nullptr, "Could not find mesh, do you parameters for lagrange_orders contain the provided aLagrangeOrder?" );
    }

    //-----------------------------------------------------------------------------

    Mesh::Mesh(
            std::shared_ptr< Database > aDatabase,
            uint                        aLagrangeMeshIndex )
    {
        // copy database pointer
        mDatabase = aDatabase;

        MORIS_ASSERT( aLagrangeMeshIndex < mDatabase->get_number_of_lagrange_meshes(),
                "HMR::Mesh::Mesh() - Could not find mesh, Lagrange mesh index %-5i exceeds number of Lagrange meshes. Check input file.",
                aLagrangeMeshIndex );

        mMesh = mDatabase->get_lagrange_mesh_by_index( aLagrangeMeshIndex );

        MORIS_ERROR( mMesh != nullptr, "HMR::Mesh::Mesh() - Pointer to Lagrange mesh %-5i is nullptr", aLagrangeMeshIndex );

        // mDatabase->get_background_mesh()->set_activation_pattern( aLagrangeMeshIndex );
        if ( mDatabase->is_finalized() )
        {
            setup_glb_to_local_maps();
        }
    }

    //-----------------------------------------------------------------------------

    Mesh::Mesh(
            std::shared_ptr< Database > aDatabase,
            uint                        aOrder,
            uint                        aLagrangePattern,
            uint                        aBsplinePattern )
    {
        // copy database pointer
        mDatabase = aDatabase;

        mDummyBSplineMeshes.resize( 3, nullptr );

        // Create factory
        Factory tFactory( mDatabase->get_parameters() );

        for ( uint iBspMesh = 0; iBspMesh < 3; iBspMesh++ )
        {
            // FIXME only one mesh
            mDummyBSplineMeshes( iBspMesh ) = tFactory.create_bspline_mesh(
                    mDatabase->get_background_mesh(),
                    aBsplinePattern,
                    aOrder,
                    MORIS_UINT_MAX );
        }

        // Create Lagrange mesh
        mMesh = tFactory.create_lagrange_mesh(
                mDatabase->get_background_mesh(),
                mDummyBSplineMeshes,
                aLagrangePattern,
                aOrder );

        // FIXME dont add
        mDatabase->add_lagrange_mesh( mMesh, aLagrangePattern );

        // remember active pattern
        auto tActivePattern = mDatabase->get_background_mesh()->get_activation_pattern();

        // activate output pattern
        mDatabase->get_background_mesh()->set_activation_pattern( mMesh->get_activation_pattern() );

        mMesh->update_mesh();

        mMesh->calculate_node_indices();
        mMesh->calculate_node_sharing();
        mMesh->calculate_t_matrices();

        for ( auto tMesh : mDummyBSplineMeshes )    // FIXME only one mesh
        {
            tMesh->calculate_basis_indices( mDatabase->get_communication_table() );
        }

        // reset active pattern
        if ( mDatabase->get_background_mesh()->get_activation_pattern() != tActivePattern )
        {
            mDatabase->get_background_mesh()->set_activation_pattern( tActivePattern );
        }
    }

    //-----------------------------------------------------------------------------

    Mesh::Mesh(
            std::shared_ptr< Database > aDatabase,
            uint                        aLagrangeOrder,
            uint                        aLagrangePattern,
            uint                        aBSplineOrder,
            uint                        aBSplinePattern )
    {
        // copy database pointer
        mDatabase = aDatabase;

        mDummyBSplineMeshes.resize( 1, nullptr );
        Factory tFactory( mDatabase->get_parameters() );

        for ( uint iBspMesh = 0; iBspMesh < 1; iBspMesh++ )
        {
            // FIXME only one mesh
            mDummyBSplineMeshes( iBspMesh ) = tFactory.create_bspline_mesh(
                    mDatabase->get_background_mesh(),
                    aBSplinePattern,
                    aBSplineOrder,
                    MORIS_UINT_MAX );
        }

        // Create Lagrange mesh
        mMesh = tFactory.create_lagrange_mesh(
                mDatabase->get_background_mesh(),
                mDummyBSplineMeshes,
                aLagrangePattern,
                aLagrangeOrder );

        // FIXME dont add
        mDatabase->add_lagrange_mesh( mMesh, aLagrangePattern );

        // remember active pattern
        auto tActivePattern = mDatabase->get_background_mesh()->get_activation_pattern();

        // activate output pattern
        mDatabase->get_background_mesh()->set_activation_pattern( mMesh->get_activation_pattern() );

        // mMesh->update_mesh();

        mMesh->calculate_node_indices();
        mMesh->calculate_node_sharing();
        mMesh->calculate_t_matrices();

        for ( auto tMesh : mDummyBSplineMeshes )    // FIXME only one mesh
        {
            tMesh->calculate_basis_indices( mDatabase->get_communication_table() );
        }

        // reset active pattern
        if ( mDatabase->get_background_mesh()->get_activation_pattern() != tActivePattern )
        {
            mDatabase->get_background_mesh()->set_activation_pattern( tActivePattern );
        }
    }

    //-----------------------------------------------------------------------------

    Mesh::~Mesh()
    {
        // Delete dummy B-spline meshes
        for ( auto tMesh : mDummyBSplineMeshes )
        {
            delete tMesh;
        }

        // Clear dummy B-spline meshes
        mDummyBSplineMeshes.clear();
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_communication_table() const
    {
        return mDatabase->get_communication_table();
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_proc_neighbors() const
    {
        return mDatabase->get_proc_neighbors();
    }

    //-----------------------------------------------------------------------------

    std::shared_ptr< Field >
    Mesh::create_field(
            const std::string& aLabel,
            uint               aBSplineIndex )
    {
        // fixme: this is not the best solution. See also
        // https://forum.libcinder.org/topic/solution-calling-shared-from-this-in-the-constructor

        // create temporary weak pointer so that shared from this works
        auto tWptr = std::shared_ptr< Mesh >( this, []( Mesh* ) {} );

        // create field
        return std::make_shared< Field >(
                aLabel,
                this->shared_from_this(),
                aBSplineIndex,
                mDatabase,
                mMesh );
    }

    //-----------------------------------------------------------------------------
    // MTK
    //-----------------------------------------------------------------------------

    uint
    Mesh::get_spatial_dim() const
    {
        return mDatabase->get_parameters()->get_number_of_dimensions();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_order()
    {
        return mMesh->get_order();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_discretization_order( uint aDiscretizationIndex )
    {
        return mMesh->get_bspline_order( aDiscretizationIndex );
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_entities(
            const mtk::EntityRank aEntityRank,
            const moris_index     aIndex ) const
    {
        if ( mMesh->get_activation_pattern() != mDatabase->get_background_mesh()->get_activation_pattern() )
        {
            mMesh->select_activation_pattern();
        }
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return this->get_num_nodes();
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                return this->get_num_edges();
                break;
            }
            case mtk::EntityRank::FACE:
            {
                return this->get_num_faces();
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return this->get_num_elems();
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                return this->get_max_num_coeffs_on_proc( aIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "unknown entity rank" );
                return 0;
                break;
            }
        }
    }
    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_elements_including_aura() const
    {
        if ( mMesh->get_activation_pattern() != mDatabase->get_background_mesh()->get_activation_pattern() )
        {
            mMesh->select_activation_pattern();
        }
        return mMesh->get_number_of_elements_including_aura();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_nodes() const
    {
        return mMesh->get_number_of_nodes_on_proc();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_edges() const
    {
        if ( mMesh->get_activation_pattern() != mDatabase->get_background_mesh()->get_activation_pattern() )
        {
            mMesh->select_activation_pattern();
        }
        return mMesh->get_number_of_edges();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_faces() const
    {
        if ( mMesh->get_activation_pattern() != mDatabase->get_background_mesh()->get_activation_pattern() )
        {
            mMesh->select_activation_pattern();
        }
        return mMesh->get_number_of_facets();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_elems() const
    {
        if ( mMesh->get_activation_pattern() != mDatabase->get_background_mesh()->get_activation_pattern() )
        {
            mMesh->select_activation_pattern();
        }
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            return this->get_num_elements_including_aura();
        }
        else
        {
            return mMesh->get_number_of_elements();
        }
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_element_indices_in_block_set( uint aSetIndex )
    {
        // Get number of elements
        uint tNumElements = this->get_num_elems();

        // Initialize element indices
        Matrix< IndexMat > tElementIndices( tNumElements, 1 );

        // Fill element indices
        if ( aSetIndex == 0 )
        {
            for ( uint tElementIndex = 0; tElementIndex < tNumElements; tElementIndex++ )
            {
                tElementIndices( tElementIndex ) = tElementIndex;
            }
        }
        else    // Should not be called for a set index other than 0, but just in case:
        {
            tElementIndices.set_size( 0, 0 );
        }

        return tElementIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_element_ids_in_block_set( uint aSetIndex )
    {
        // Get number of elements
        uint tNumElements = this->get_num_elems();

        // Initialize element IDs
        Matrix< IdMat > tElementIDs( tNumElements, 1 );

        // Fill element IDs
        if ( aSetIndex == 0 )
        {
            for ( uint tElementIndex = 0; tElementIndex < tNumElements; tElementIndex++ )
            {
                if ( mDatabase->get_parameters()->use_number_aura() and    //
                        mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
                {
                    tElementIDs( tElementIndex ) = mMesh->get_element_including_aura( tElementIndex )->get_hmr_id();
                }
                else
                {
                    tElementIDs( tElementIndex ) = mMesh->get_element( tElementIndex )->get_hmr_id();
                }
            }
        }
        else    // Should not be called for a set index other than 0, but just in case:
        {
            tElementIDs.set_size( 0, 0 );
        }

        return tElementIDs;
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_max_num_coeffs_on_proc( const uint aBsplineMeshIndex ) const
    {
        return mMesh->get_number_of_bsplines_on_proc( aBsplineMeshIndex );
    }

    //-----------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Mesh::get_t_matrix_of_node_loc_ind(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        return *mMesh->get_node_by_index( aNodeIndex )->get_interpolation( aBSplineMeshIndex )->get_weights();
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_coefficient_indices_of_node(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_interpolation( aBSplineMeshIndex )->get_indices();
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_coefficient_IDs_of_node(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_interpolation( aBSplineMeshIndex )->get_ids();
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_coefficient_owners_of_node(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_interpolation( aBSplineMeshIndex )->get_owners();
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Mesh::get_coefficient_ijkl_IDs_of_node(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_interpolation( aBSplineMeshIndex )->get_ijkl_id();
    }

    //-----------------------------------------------------------------------------
    Matrix< IndexMat >
    Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aOutputEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                switch ( aInputEntityRank )
                {
                    case mtk::EntityRank::NODE:
                    {
                        return this->get_nodes_connected_to_node_loc_inds( aEntityIndex );
                        break;
                    }
                    case mtk::EntityRank::EDGE:
                    {
                        Matrix< IndexMat > tNodeToEdge = this->get_nodes_connected_to_edge_loc_inds( aEntityIndex );
                        return tNodeToEdge;

                        break;
                    }
                    case mtk::EntityRank::FACE:
                    {
                        Matrix< IndexMat > tNodeToFace = this->get_nodes_connected_to_face_loc_inds( aEntityIndex );
                        return tNodeToFace;
                        break;
                    }
                    case mtk::EntityRank::ELEMENT:
                    {
                        Matrix< IndexMat > tNodeToElement = this->get_nodes_connected_to_element_loc_inds( aEntityIndex );
                        return tNodeToElement;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                }
                break;
            }    // end output rank node
            case mtk::EntityRank::EDGE:
            {
                switch ( aInputEntityRank )
                {
                    case mtk::EntityRank::NODE:
                    {
                        return this->get_edges_connected_to_node_loc_inds( aEntityIndex );
                        break;
                    }
                    case mtk::EntityRank::EDGE:
                    {
                        MORIS_ERROR( false, "HMR does not provide edge to edge connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                    case mtk::EntityRank::FACE:
                    {
                        MORIS_ERROR( false, "HMR does not provide edge to face connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                    case mtk::EntityRank::ELEMENT:
                    {
                        if ( this->get_spatial_dim() == 3 )
                        {
                            return this->get_edges_connected_to_element_loc_inds( aEntityIndex );
                        }
                        else
                        {
                            return this->get_faces_connected_to_element_loc_inds( aEntityIndex );
                        }
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                }
                break;
            }    // end output rank edge
            case mtk::EntityRank::FACE:
            {
                switch ( aInputEntityRank )
                {
                    case mtk::EntityRank::NODE:
                    {
                        return this->get_faces_connected_to_node_loc_inds( aEntityIndex );
                        break;
                    }
                    case mtk::EntityRank::EDGE:
                    {
                        MORIS_ERROR( false, "HMR does not provide face to edge connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                    case mtk::EntityRank::FACE:
                    {
                        MORIS_ERROR( false, "HMR does not provide face to face connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                    case mtk::EntityRank::ELEMENT:
                    {
                        return this->get_faces_connected_to_element_loc_inds( aEntityIndex );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                }
                break;
            }    // end output rank face
            case mtk::EntityRank::ELEMENT:
            {
                switch ( aInputEntityRank )
                {
                    case mtk::EntityRank::NODE:
                    {
                        return this->get_elements_connected_to_node_loc_inds( aEntityIndex );
                        break;
                    }
                    case mtk::EntityRank::EDGE:
                    {
                        MORIS_ERROR( false, "HMR does not provide element to edge connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                    case mtk::EntityRank::FACE:
                    {
                        return this->get_elements_connected_to_face_loc_inds( aEntityIndex );
                        break;
                    }
                    case mtk::EntityRank::ELEMENT:
                    {
                        return this->get_elements_connected_to_element_and_face_ind_loc_inds( aEntityIndex );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                }
                break;
            }    // end output rank element
            case mtk::EntityRank::BSPLINE:
            {
                switch ( aInputEntityRank )
                {
                    case mtk::EntityRank::ELEMENT:
                    {
                        return this->get_inds_of_active_elements_connected_to_basis( mMesh->get_bspline_mesh( aIndex )
                                                                                             ->get_basis_by_index( aEntityIndex ) );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                        return Matrix< IndexMat >( 0, 0 );
                        break;
                    }
                }
            }    // end output rank linear bspline
            default:
            {
                MORIS_ERROR( false, "HMR does not provide the requested connectivity" );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
        }
    }

    //-----------------------------------------------------------------------------
    Matrix< IndexMat >
    Mesh::get_entity_connected_to_entity_glob_ids(
            moris_index       aEntityId,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {

        moris_index tEntityIndex = this->get_loc_entity_ind_from_entity_glb_id( aEntityId, aInputEntityRank );

        Matrix< IndexMat > tEntityToEntity = this->get_entity_connected_to_entity_loc_inds( tEntityIndex, aInputEntityRank, aOutputEntityRank, aIndex );

        for ( uint i = 0; i < tEntityToEntity.numel(); i++ )
        {
            tEntityToEntity( i ) = this->get_glb_entity_id_from_entity_loc_index( tEntityToEntity( i ), aOutputEntityRank, aIndex );
        }

        return tEntityToEntity;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_nodes_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        // get pointer to basis
        const Basis* tBasis = mMesh->get_node_by_index( aNodeIndex );

        // get number of connected elements
        uint tNumberOfElements = tBasis->get_element_counter();

        // get number of nodes connected to element
        uint tCount = 0;
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            tCount += ( tBasis->get_element( e )->get_number_of_vertices() - 1 );
        }

        // allocate temporary Matrix
        Matrix< IndexMat > tNodeIndices( tCount, 1 );

        // reset counter
        tCount = 0;

        // get ID of this basis
        auto tMyID = tBasis->get_hmr_id();

        // loop over all elements
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            // get pointer to element
            const Element* tElement = tBasis->get_element( e );

            // ask element about number of nodes
            uint tNumberOfVertices = tElement->get_number_of_vertices();

            // loop over all connected vertices
            for ( uint k = 0; k < tNumberOfVertices; ++k )
            {
                // test if this vertex is not myself
                if ( tElement->get_basis( k )->get_hmr_id() != tMyID )
                {
                    // add basis index to Indices
                    tNodeIndices( tCount++ ) = tElement->get_basis( k )->get_index();
                }
            }
        }

        // make result unique
        Matrix< IndexMat > aNodeIndices;
        unique( tNodeIndices, aNodeIndices );
        return aNodeIndices;
    }
    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_nodes_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const
    {
        return mMesh->get_edge( aEdgeIndex )->get_vertex_inds();
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_nodes_connected_to_face_loc_inds( moris_index aFaceIndex ) const
    {
        return mMesh->get_facet( aFaceIndex )->get_vertex_inds();
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const
    {
        Element* tElement;

        // get pointer to element
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            tElement = mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            tElement = mMesh->get_element( aElementIndex );
        }

        // get number of nodes
        uint tNumberOfNodes = tElement->get_number_of_vertices();

        // allocate output
        Matrix< IndexMat > aIndices( tNumberOfNodes, 1 );

        // populate output
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            aIndices( k ) = tElement->get_basis( k )->get_index();
        }

        return aIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        // get pointer to basis
        Basis* tBasis = mMesh->get_node_by_index( aNodeIndex );

        uint tNumberOfEdges = tBasis->get_edge_counter();

        Matrix< IndexMat > tEdgeIndex( tNumberOfEdges, 1 );

        for ( uint k = 0; k < tNumberOfEdges; ++k )
        {
            tEdgeIndex( k ) = tBasis->get_edge( k )->get_index();
        }

        Matrix< IndexMat > aEdgeIndex;
        sort( tEdgeIndex, aEdgeIndex );
        return aEdgeIndex;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const
    {
        // get pointer to element
        Element* tElement;

        // get pointer to element
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            tElement = mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            tElement = mMesh->get_element( aElementIndex );
        }

        // will be 0 for 2D, 12 for 3D
        uint tNumberOfEdges = tElement->get_background_element()->get_number_of_edges();

        // allocate output
        Matrix< IndexMat > aIndices( tNumberOfEdges, 1 );

        // populate output
        for ( uint e = 0; e < tNumberOfEdges; ++e )
        {
            aIndices( e ) = tElement->get_hmr_edge( e )->get_index();
        }

        return aIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        // get pointer to basis
        Basis* tBasis = mMesh->get_node_by_index( aNodeIndex );

        uint tNumberOfFacets = tBasis->get_facet_counter();

        Matrix< IndexMat > tFaceIndex( tNumberOfFacets, 1 );

        for ( uint k = 0; k < tNumberOfFacets; ++k )
        {
            tFaceIndex( k ) = tBasis->get_facet( k )->get_index();
        }

        Matrix< IndexMat > aFaceIndex;
        sort( tFaceIndex, aFaceIndex );
        return aFaceIndex;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const
    {
        // get pointer to element
        Element* tElement;

        // get pointer to element
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            tElement = mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            tElement = mMesh->get_element( aElementIndex );
        }

        uint tNumberOfFaces = tElement->get_background_element()->get_number_of_facets();

        Matrix< IndexMat > aIndices( tNumberOfFaces, 1 );

        for ( uint f = 0; f < tNumberOfFaces; ++f )
        {
            aIndices( f ) = tElement->get_hmr_facet( f )->get_index();
        }

        return aIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        // collect memory indices of active elements
        Matrix< DDLUMat > tMemoryIndices;
        this->collect_memory_indices_of_active_elements_connected_to_node( aNodeIndex,
                tMemoryIndices );

        Matrix< IndexMat > aIndices;
        this->get_element_indices_from_memory_indices( tMemoryIndices,
                aIndices );

        return aIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const
    {
        // get pointer to facet
        Facet* tFacet = mMesh->get_facet( aFaceIndex );

        Matrix< IndexMat > aIndices;

        if ( tFacet->get_follower() == nullptr )
        {
            if ( tFacet->get_hmr_leader()->is_active() )
            {
                aIndices.set_size( 1, 1 );
                aIndices( 0 ) = tFacet->get_hmr_leader()->get_index();
            }
        }
        else
        {
            Element* tLeader   = tFacet->get_hmr_leader();
            Element* tFollower = tFacet->get_hmr_follower();

            if ( tLeader->is_active() && tFollower->is_active() )
            {
                aIndices.set_size( 2, 1 );
                aIndices( 0 ) = tLeader->get_index();
                aIndices( 1 ) = tFollower->get_index();
            }
            else if ( tLeader->is_active() )
            {
                aIndices.set_size( 1, 1 );
                aIndices( 0 ) = tLeader->get_index();
            }
            else if ( tFollower->is_active() )
            {
                aIndices.set_size( 1, 1 );
                aIndices( 0 ) = tFollower->get_index();
            }
        }

        return aIndices;
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_elements_in_support_of_basis(
            const uint          aMeshIndex,
            const uint          aBasisIndex,
            Matrix< IndexMat >& aElementIndices )
    {
        mMesh->get_my_elements_in_basis_support( aMeshIndex, aBasisIndex, aElementIndices );
    }

    //-------------------------------------------------------------------------------

    uint
    Mesh::get_num_basis_functions( const uint aMeshIndex )
    {
        BSpline_Mesh_Base* tMesh = mMesh->get_bspline_mesh( aMeshIndex );

        if ( tMesh != nullptr )
        {
            return tMesh->get_number_of_indexed_basis();
        }
        else
        {
            return this->get_num_nodes();
        }
    }

    //-------------------------------------------------------------------------------

    void
    Mesh::get_nodes_indices_in_bounding_box(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aBoundingBoxSize,
            Matrix< IndexMat >&     aNodeIndices )
    {
        mMesh->calculate_nodes_indices_in_bounding_box( aPoint, aBoundingBoxSize, aNodeIndices );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_elements_in_bspline_element(
            moris_index const          aBspElementIndex,
            moris_index const          aDiscretizationMeshIndex,
            Vector< mtk::Cell* >& aCells )
    {
        mMesh->get_elements_in_bspline_element(
                aBspElementIndex,
                aDiscretizationMeshIndex,
                aCells );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_lagrange_elements_in_bspline_elements(
            moris_index const                          aDiscretizationMeshIndex,
            Vector< Vector< mtk::Cell* > >&  aCells,
            Vector< Vector< moris_index > >& aCellIndices,
            Vector< moris_index >&                aLagToBspCellIndices,
            Vector< uint >&                       aBspCellRefineLevels,
            Vector< mtk::Cell* >&                 aBspCells )
    {
        mMesh->get_lagrange_elements_in_bspline_elements(
                aDiscretizationMeshIndex,
                aCells,
                aCellIndices,
                aLagToBspCellIndices,
                aBspCellRefineLevels,
                aBspCells );
    }

    // ----------------------------------------------------------------------------

    const luint*
    Mesh::get_bspline_element_ijk_level(
            moris_index      aDiscretizationMeshIndex,
            const mtk::Cell* aBsplineElement,
            uint&            aLevel )
    {
        // get the HMR cell
        const Element* aHMRCell = dynamic_cast< const Element* >( aBsplineElement );

        // get level and ijk
        aLevel = aHMRCell->get_level();
        return aHMRCell->get_ijk();
    }

    // ----------------------------------------------------------------------------

    void
    Mesh::get_extended_t_matrix(
            moris_index                                 aDiscretizationMeshIndex,
            moris_index                                 aBSplineCellIndex,
            moris::mtk::Cell&                           aLagrangeCell,
            Vector< Vector< mtk::Vertex* > >& tBsplineBasis,
            Vector< Matrix< DDRMat > >&            tWeights )
    {
        Element& aHMRLagrangeCell = dynamic_cast< Element& >( aLagrangeCell );
        mMesh->get_extended_t_matrix( aDiscretizationMeshIndex,
                aBSplineCellIndex,
                aHMRLagrangeCell,
                tBsplineBasis,
                tWeights );
    }

    // ----------------------------------------------------------------------------

    void
    Mesh::get_L2_projection_matrix(
            moris_index                                       aDiscretizationMeshIndex,
            const mtk::Cell*                                  aRootBSplineCell,
            const mtk::Cell*                                  aExtendedBSplineCell,
            Vector< Vector< const mtk::Vertex* > >& tRootBsplineBasis,
            Vector< const mtk::Vertex* >&                tExtendedBsplineBasis,
            Vector< Matrix< DDRMat > >&                  tWeights )
    {
        const Element* aHMRRootCell     = dynamic_cast< const Element* >( aRootBSplineCell );
        const Element* aHMRExtendedCell = dynamic_cast< const Element* >( aExtendedBSplineCell );

        mMesh->get_L2_projection_matrix( aDiscretizationMeshIndex,
                aHMRRootCell,
                aHMRExtendedCell,
                tRootBsplineBasis,
                tExtendedBsplineBasis,
                tWeights );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_elements_in_interpolation_cluster(
            moris_index                aElementIndex,
            moris_index                aDiscretizationMeshIndex,
            Vector< mtk::Cell* >& aCells )
    {
        mMesh->get_elements_in_interpolation_cluster(
                aElementIndex,
                aDiscretizationMeshIndex,
                aCells );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_elements_in_bspline_element_and_side_ordinal(
            moris_index const          aBsplineElementIndex,
            moris_index const          aDiscretizationMeshIndex,
            moris_index const          aSideOrdinal,
            Vector< mtk::Cell* >& aCells )
    {
        mMesh->get_elements_in_bspline_element_and_side_ordinal(
                aBsplineElementIndex,
                aDiscretizationMeshIndex,
                aSideOrdinal,
                aCells );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_elements_in_interpolation_cluster_and_side_ordinal(
            moris_index const          aElementIndex,
            moris_index const          aDiscretizationMeshIndex,
            moris_index const          aSideOrdinal,
            Vector< mtk::Cell* >& aCells )
    {
        mMesh->get_elements_in_interpolation_cluster_and_side_ordinal(
                aElementIndex,
                aDiscretizationMeshIndex,
                aSideOrdinal,
                aCells );
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        // collect memory indices of active neighbors
        Matrix< DDLUMat > tMemoryIndices;
        Matrix< DDLUMat > tThisCellFacetsOrds;
        Matrix< DDLUMat > tNeighborCellFacetsOrds;
        Matrix< DDLUMat > tTransitionNeighborCellLocation;
        luint             tNumberOfNeighbors;

        this->collect_memory_indices_of_active_element_neighbors( aElementIndex,
                tMemoryIndices,
                tThisCellFacetsOrds,
                tNeighborCellFacetsOrds,
                tTransitionNeighborCellLocation,
                tNumberOfNeighbors );
        // reserve memory for output matrix
        Matrix< IndexMat > tIndices( 2, tNumberOfNeighbors );

        // my cell
        Element* tMyCell = nullptr;
        if ( mDatabase->get_parameters()->use_number_aura()    //
                and mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            tMyCell = mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            tMyCell = mMesh->get_element( aElementIndex );
        }

        // copy indices from pointers
        for ( luint k = 0; k < tNumberOfNeighbors; ++k )
        {

            Element* tOtherCell = mMesh->get_element_by_memory_index( tMemoryIndices( k ) );

            MORIS_ASSERT( tOtherCell->get_index() != MORIS_INDEX_MAX, "A Neighbor Cell with max index is about to be outputted." );
            tIndices( 0, k ) = tOtherCell->get_index();

            // verify the facet neighborhood makes sense
            Facet* tFacet = tMyCell->get_hmr_facet( tThisCellFacetsOrds( k ) );
            MORIS_ASSERT( tFacet != nullptr, "get_elements_connected_to_element_and_face_ind_loc_inds(), facet is nullptr" );

            tIndices( 1, k ) = tFacet->get_index();
        }
        return tIndices;
    }

    //-----------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const
    {
        // collect memory indices of active neighbors
        Matrix< DDLUMat > tMemoryIndices;
        Matrix< DDLUMat > tThisCellFacetsOrds;
        Matrix< DDLUMat > tNeighborCellFacetsOrds;
        Matrix< DDLUMat > tTransitionNeighborCellLocation;
        luint             tNumberOfNeighbors;

        this->collect_memory_indices_of_active_element_neighbors( aElementIndex,
                tMemoryIndices,
                tThisCellFacetsOrds,
                tNeighborCellFacetsOrds,
                tTransitionNeighborCellLocation,
                tNumberOfNeighbors );

        // reserve memory for output matrix
        Matrix< IndexMat > tNeighborIndicesAndSideOrds( 4, tNumberOfNeighbors );

        // copy indices from pointers
        for ( luint k = 0; k < tNumberOfNeighbors; ++k )
        {

            Element* tOtherCell = mMesh->get_element_by_memory_index( tMemoryIndices( k ) );

            MORIS_ASSERT( tOtherCell->get_index() != MORIS_INDEX_MAX, "A Neighbor Cell with max index is about to be outputted." );
            //                MORIS_ASSERT(tOtherCell->get_level() == mMesh->get_element( aElementIndex )->get_level()     ||
            //                             tOtherCell->get_level() == mMesh->get_element( aElementIndex )->get_level() - 1 ||
            //                             tOtherCell->get_level() == mMesh->get_element( aElementIndex )->get_level() + 1   ,"Neighboring cells must be on the within 1 level of the provided cell.");

            tNeighborIndicesAndSideOrds( 0, k ) = tOtherCell->get_index();
            tNeighborIndicesAndSideOrds( 1, k ) = tThisCellFacetsOrds( k );
            tNeighborIndicesAndSideOrds( 2, k ) = tNeighborCellFacetsOrds( k );
            tNeighborIndicesAndSideOrds( 3, k ) = tTransitionNeighborCellLocation( k );
        }
        return tNeighborIndicesAndSideOrds;
    }

    //-----------------------------------------------------------------------------

    bool
    Mesh::get_elements_connected_to_element_through_face_ord(
            moris_index                 aElementIndex,
            moris_index                 aSideOrdinal,
            moris_index&                aMyRefineLevel,
            Vector< moris_index >& aNeighborElements,
            Vector< moris_index >& aNeighborSideOrdinals,
            Vector< moris_index >& aTransitionLocations,
            Vector< moris_index >& aNeighborRefinementLevels ) const
    {
        // get the current element's refinement level
        Element* tElement = mMesh->get_element_including_aura( aElementIndex );
        aMyRefineLevel    = tElement->get_level();

        // initialize output list for active neighbor search
        Matrix< DDLUMat > tMemoryIndices;
        Matrix< DDLUMat > tThisCellFacetsOrds;
        Matrix< DDLUMat > tNeighborCellFacetsOrds;
        Matrix< DDLUMat > tTransitionNeighborCellLocation;
        luint             tNumberOfNeighbors;

        // collect memory indices of active neighbors
        this->collect_memory_indices_of_active_element_neighbors(
                aElementIndex,
                tMemoryIndices,
                tThisCellFacetsOrds,
                tNeighborCellFacetsOrds,
                tTransitionNeighborCellLocation,
                tNumberOfNeighbors );

        // initialize counter
        uint                tNumNeighborsOnSideOrd = 0;
        Vector< uint > tValidNeighbors( 0 );

        // find the number of neighbors on the given side ordinal
        for ( uint iNeighbor = 0; iNeighbor < tNumberOfNeighbors; iNeighbor++ )
        {
            if ( tThisCellFacetsOrds( iNeighbor ) == (uint)aSideOrdinal )
            {
                tNumNeighborsOnSideOrd++;
                tValidNeighbors.push_back( iNeighbor );
            }
        }

        // check that there are actually any neighbors
        MORIS_ASSERT( tNumNeighborsOnSideOrd > 0,
                "Mesh::get_elements_connected_to_element_through_face_ord() - No neighbors found on side ordinal." );

        // mark transition as trivial or non-trivial (non-trivial case 1: transition from big to small)
        bool tNonTrivialFlag = ( tNumNeighborsOnSideOrd > 1 );

        // initialize output lists with correct size
        aNeighborElements.resize( tNumNeighborsOnSideOrd );
        aNeighborSideOrdinals.resize( tNumNeighborsOnSideOrd );
        aTransitionLocations.resize( tNumNeighborsOnSideOrd );
        aNeighborRefinementLevels.resize( tNumNeighborsOnSideOrd );

        // find all neighbors on the given side ordinal
        for ( uint iNeighbor = 0; iNeighbor < tNumNeighborsOnSideOrd; iNeighbor++ )
        {
            // get local neighbor index
            uint tLocalNeighborIndex = tValidNeighbors( iNeighbor );

            // get access to the neighbor element
            Element* tOtherCell = mMesh->get_element_by_memory_index( tMemoryIndices( tLocalNeighborIndex ) );

            // check validity of neighbor element
            MORIS_ASSERT( tOtherCell->get_index() != MORIS_INDEX_MAX,
                    "Mesh::get_elements_connected_to_element_through_face_ord() - A Neighbor Cell with max index is about to be outputted." );

            aNeighborElements( iNeighbor )         = tOtherCell->get_index();
            aNeighborSideOrdinals( iNeighbor )     = tNeighborCellFacetsOrds( tLocalNeighborIndex );
            aTransitionLocations( iNeighbor )      = tTransitionNeighborCellLocation( tLocalNeighborIndex );
            aNeighborRefinementLevels( iNeighbor ) = tOtherCell->get_level();
        }

        // check whether the current element is finer than its neighbor, if so use the neighbors transition location to the current element
        if ( aMyRefineLevel > aNeighborRefinementLevels( 0 ) )
        {
            // run some checks that the data makes sense
            MORIS_ASSERT( aTransitionLocations.size() == 1,
                    "hmr::Mesh::get_elements_connected_to_element_through_face_ord() - "
                    "Neighbor element coarser, but multiple neighbor elements listed. This does not make sense." );

            // mark transition as non-trivial (case 2: transition from small to big)
            tNonTrivialFlag = true;

            // get the neighbor's side ordinal
            moris_index tBigElementsSideOrdinal = aNeighborSideOrdinals( 0 );

            // find the child ordinal that the current element sits at, if it has one
            Background_Element_Base* tBackElement      = tElement->get_background_element();
            moris_index              tMyOctreePosition = tElement->get_background_element()->get_child_index();

            // get the number of spatial dimensions
            uint tNumSpatialDims = tBackElement->get_number_of_facets() / 2;
            MORIS_ASSERT( tNumSpatialDims == 2 || tNumSpatialDims == 3,
                    "hmr::Mesh::get_elements_connected_to_element_through_face_ord() - "
                    "Number of spatial dimensions not 2 or 3." );

            // get the transition location from the coarser neighbor element's point of view
            aTransitionLocations( 0 ) = get_refinement_transition_location_for_neighbor_child_ordinal(
                    tBigElementsSideOrdinal,
                    tMyOctreePosition,
                    tNumSpatialDims );
        }

        // return whether element transition is non-trivial
        return tNonTrivialFlag;
    }

    //-----------------------------------------------------------------------------
    // private
    //-----------------------------------------------------------------------------

    void
    Mesh::get_element_indices_from_memory_indices(
            const Matrix< DDLUMat >& aMemoryIndices,
            Matrix< IndexMat >&      aIndices ) const
    {
        // get length of vector
        uint tNumberOfElements = aMemoryIndices.length();

        // allocate memory
        aIndices.set_size( tNumberOfElements, 1 );

        // copy indices from pointers
        for ( uint k = 0; k < tNumberOfElements; ++k )
        {
            aIndices( k ) = mMesh->get_element_by_memory_index( aMemoryIndices( k ) )->get_index();
        }
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::collect_memory_indices_of_active_element_neighbors(
            const moris_index  aElementIndex,
            Matrix< DDLUMat >& aMemoryIndices,
            Matrix< DDLUMat >& aThisCellFacetOrds,
            Matrix< DDLUMat >& aNeighborCellFacetOrds,
            Matrix< DDLUMat >& aTransitionNeighborCellLocation,
            luint&             aCounter ) const
    {
        // get active index of this mesh
        uint tPattern = mMesh->get_activation_pattern();

        // get pointer to element
        Element* tElement;

        // get pointer to element
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            tElement = mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            tElement = mMesh->get_element( aElementIndex );
        }

        // get pointer to background element
        Background_Element_Base* tBackElement = tElement->get_background_element();

        // count elements that are not padding
        aCounter = 0;

        uint tNumberOfFacets = tBackElement->get_number_of_facets();

        for ( uint k = 0; k < tNumberOfFacets; ++k )
        {
            // get pointer to neighbor
            Background_Element_Base* tNeighbor = tBackElement->get_neighbor( k );

            if ( tNeighbor != nullptr )
            {
                if ( !tNeighbor->is_padding() )
                {
                    if ( tNeighbor->is_refined( tPattern ) )
                    {
                        tNeighbor->get_number_of_active_descendants( tPattern, aCounter );
                    }
                    else
                    {
                        // increment counter
                        ++aCounter;
                    }
                }
            }
        }

        // allocate matrix with memory indices
        aMemoryIndices.set_size( aCounter, 1 );
        aMemoryIndices.fill( MORIS_INDEX_MAX );

        aThisCellFacetOrds.set_size( aCounter, 1 );
        aThisCellFacetOrds.fill( MORIS_INDEX_MAX );

        aNeighborCellFacetOrds.set_size( aCounter, 1 );
        aNeighborCellFacetOrds.fill( MORIS_INDEX_MAX );

        aTransitionNeighborCellLocation.set_size( aCounter, 1 );
        aTransitionNeighborCellLocation.fill( MORIS_INDEX_MAX );

        // reset counter
        aCounter = 0;

        // number of descendance
        uint tStart = 0;
        uint tEnd   = 0;

        // loop over all neighbors
        for ( uint k = 0; k < tNumberOfFacets; ++k )
        {
            // get pointer to neighbor
            Background_Element_Base* tNeighbor = tBackElement->get_neighbor( k );

            // get the neighbor side ordinal
            int tNeighborSideOrd = tBackElement->get_neighbor_side_ordinal( k );

            if ( tNeighbor != nullptr )
            {
                if ( !tNeighbor->is_padding() )
                {
                    // if there's multiple neighbors on this facet
                    if ( tNeighbor->is_refined( tPattern ) )
                    {
                        // get the neighbor child cell ordinal
                        // int tNeighborSideOrd = tBackElement->get_neighbor_side_ordinal(k);

                        // get my child cell ordinals that would be on this side (for use in XTK ghost)
                        Matrix< IndexMat > tMyChildOrds;
                        tBackElement->get_child_cell_ordinals_on_side( k, tMyChildOrds );

                        tStart = aCounter;

                        tNeighbor->collect_active_descendants_by_memory_index(
                                tPattern,
                                aMemoryIndices,
                                aCounter,
                                k );

                        // mark facets that we share with other element
                        tEnd = aCounter;

                        uint tZeroStart = 0;
                        for ( uint i = tStart; i < tEnd; i++ )
                        {
                            aThisCellFacetOrds( i )              = k;
                            aNeighborCellFacetOrds( i )          = tNeighborSideOrd;
                            aTransitionNeighborCellLocation( i ) = tMyChildOrds( tZeroStart++ );
                        }
                    }
                    // if there's only one neighbor on this facet
                    else
                    {
                        while ( !tNeighbor->is_active( tPattern ) )
                        {
                            tNeighbor = tNeighbor->get_parent();
                        }

                        MORIS_ASSERT( tNeighbor->get_memory_index() != MORIS_INDEX_MAX, "A Neighbor with max index should not be outputted." );

                        aThisCellFacetOrds( aCounter )     = k;
                        aNeighborCellFacetOrds( aCounter ) = tNeighborSideOrd;
                        aMemoryIndices( aCounter++ )       = tNeighbor->get_memory_index();
                    }
                }
            }
        }

        aMemoryIndices.resize( aCounter, 1 );
        aThisCellFacetOrds.resize( aCounter, 1 );
        aNeighborCellFacetOrds.resize( aCounter, 1 );
    }

    //-----------------------------------------------------------------------------

    void
    Mesh::collect_memory_indices_of_active_elements_connected_to_node(
            const moris_index  aNodeIndex,
            Matrix< DDLUMat >& aMemoryIndices ) const
    {
        // get active index of this mesh
        uint tPattern = mMesh->get_activation_pattern();

        // get pointer to node
        Basis* tNode = mMesh->get_node_by_index( aNodeIndex );

        // get number of elements
        luint tElementCounter = tNode->get_element_counter();

        // element counter with active elements
        luint tCount = 0;

        for ( uint k = 0; k < tElementCounter; ++k )
        {
            if ( tNode->get_element( k )->get_background_element()->is_active( tPattern ) )
            {
                // increment counter
                ++tCount;
            }
        }

        // allocate matrix with memory indices
        aMemoryIndices.set_size( tCount, 1 );

        // reset counter
        tCount = 0;

        for ( uint k = 0; k < tElementCounter; ++k )
        {
            // get pointer to element
            Background_Element_Base* tBackElement = tNode->get_element( k )
                                                            ->get_background_element();

            if ( tBackElement->is_active( tPattern ) )
            {
                // increment counter
                aMemoryIndices( tCount++ ) = tBackElement->get_memory_index();
            }
        }
    }

    //-----------------------------------------------------------------------------

    Matrix< DDRMat >
    Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_coords();
    }

    //-----------------------------------------------------------------------------

    moris_id
    Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_node_by_index( aEntityIndex )->get_id();
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                return mMesh->get_edge( aEntityIndex )->get_id();
                break;
            }
            case mtk::EntityRank::FACE:
            {
                return mMesh->get_facet( aEntityIndex )->get_id();
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                if ( mDatabase->get_parameters()->use_number_aura() and    //
                        mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
                {
                    return mMesh->get_element_including_aura( aEntityIndex )->get_id();
                }
                else
                {
                    return mMesh->get_element( aEntityIndex )->get_id();
                }
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                if ( !mMesh->get_bspline_mesh_is_trivial_interpolation( aIndex ) )
                {
                    return mMesh->get_bspline( aIndex, aEntityIndex )->get_id();
                }
                else
                {
                    return this->get_glb_entity_id_from_entity_loc_index( aEntityIndex, mtk::EntityRank::NODE, aIndex );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "unknown entity rank" );
                return 0;
            }
        }
    }

    //-----------------------------------------------------------------------------

    moris_id
    Mesh::get_glb_element_id_from_element_loc_index( moris_index aEntityIndex ) const
    {
        return mMesh->get_element_including_aura( aEntityIndex )->get_id();
    }

    //-----------------------------------------------------------------------------

    moris_index
    Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id          aEntityId,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        auto tIter = mEntityGlobalToLocalMap( (uint)aEntityRank + aIndex ).find( aEntityId );

        MORIS_ERROR( tIter != mEntityGlobalToLocalMap( (uint)aEntityRank + aIndex ).end(),
                "HMR::Mesh::get_loc_entity_ind_from_entity_glb_id() - "
                "Provided Entity Id is not in the map. aEntityId = %u mtk::EntityRank = %u on process %u. Has the map been initialized?; Size of map: %zu",
                aEntityId,
                (uint)aEntityRank,
                par_rank(),
                mEntityGlobalToLocalMap.size() );

        return tIter->second;
    }

    //-----------------------------------------------------------------------------
    moris_id
    Mesh::get_max_entity_id(
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_max_node_id();
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                return mMesh->get_max_edge_id();
                break;
            }
            case mtk::EntityRank::FACE:
            {
                return mMesh->get_max_facet_id();
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return mMesh->get_max_element_id();
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                uint tNumEntities = this->get_num_entities( aEntityRank, aIndex );

                moris_id tLocalMaxId = 0;

                for ( uint i = 0; i < tNumEntities; i++ )
                {
                    moris_id tId = this->get_glb_entity_id_from_entity_loc_index( i, aEntityRank, aIndex );

                    if ( tId > tLocalMaxId )
                    {
                        tLocalMaxId = tId;
                    }
                }

                moris_id tGlobalMaxId = moris::max_all( tLocalMaxId );
                return tGlobalMaxId;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "unknown entity rank" );
                return 0;
            }
        }
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_node_owner( moris_index aNodeIndex ) const
    {
        return mMesh->get_node_by_index( aNodeIndex )->get_owner();
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_element_owner( moris_index aElementIndex ) const
    {
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            return mMesh->get_element_including_aura( aElementIndex )->get_owner();
        }
        else
        {
            return mMesh->get_element( aElementIndex )->get_owner();
        }
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_entity_owner(
            moris_index       aEntityIndex,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_node_by_index( aEntityIndex )->get_owner();
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                if ( this->get_spatial_dim() == 3 )
                {
                    return mMesh->get_edge( aEntityIndex )->get_owner();
                }
                else
                {
                    return mMesh->get_facet( aEntityIndex )->get_owner();
                }
                break;
            }
            case mtk::EntityRank::FACE:
            {
                return mMesh->get_facet( aEntityIndex )->get_owner();
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                if ( mDatabase->get_parameters()->use_number_aura() and    //
                        mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
                {
                    return mMesh->get_element_including_aura( aEntityIndex )->get_owner();
                }
                else
                {
                    return mMesh->get_element( aEntityIndex )->get_owner();
                }
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                if ( !mMesh->get_bspline_mesh_is_trivial_interpolation( aIndex ) )
                {
                    return mMesh->get_bspline( aIndex, aEntityIndex )->get_owner();
                }
                else
                {
                    return this->get_entity_owner( aEntityIndex, mtk::EntityRank::NODE, aIndex );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "unknown entity rank" );
                return 0;
                break;
            }
        }
    }

    //-----------------------------------------------------------------------------
    void
    Mesh::get_processors_whom_share_entity(
            moris_index      aEntityIndex,
            mtk::EntityRank  aEntityRank,
            Matrix< IdMat >& aProcsWhomShareEntity ) const
    {
        MORIS_ASSERT( par_size() == 1, "Not implemented in HMR (pending completion of entity sharing info in HMR" );

        aProcsWhomShareEntity.set_size( 0, 0 );
    }

    //-----------------------------------------------------------------------------

    mtk::EntityRank
    Mesh::get_facet_rank() const
    {
        return mtk::EntityRank::FACE;
    }

    //-----------------------------------------------------------------------------

    Vector< mtk::Vertex const * >
    Mesh::get_all_vertices() const
    {
        uint tNumVertices = this->get_num_entities( mtk::EntityRank::NODE );

        Vector< mtk::Vertex const * > tVertices( tNumVertices );

        for ( uint i = 0; i < tNumVertices; i++ )
        {
            tVertices( i ) = &get_mtk_vertex( (moris_index)i );
        }

        return tVertices;
    }

    ////-----------------------------------------------------------------------------
    //
    //        Vector<mtk::Vertex const *>
    //        Mesh::get_all_vertices_including_aura() const
    //        {
    //            uint tNumVertices = this->get_num_nodes_including_aura();
    //
    //            Vector<mtk::Vertex const *> tVertices (tNumVertices);
    //
    //            for(uint  i = 0; i < tNumVertices; i++)
    //            {
    //                tVertices(i) = &get_mtk_vertex_including_aura( ( moris_index ) i );
    //            }
    //
    //            return tVertices;
    //        }

    //-----------------------------------------------------------------------------

    void
    Mesh::get_adof_map(
            const uint                    aBSplineIndex,
            map< moris_id, moris_index >& aAdofMap ) const
    {
        aAdofMap.clear();

        moris_index tNumberOfBSplines = mMesh->get_number_of_bsplines_on_proc( aBSplineIndex );

        for ( moris_index k = 0; k < tNumberOfBSplines; ++k )
        {
            Basis* tBasis                = mMesh->get_bspline( aBSplineIndex, k );
            aAdofMap[ tBasis->get_id() ] = tBasis->get_index();
        }
    }

    //-----------------------------------------------------------------------------

    uint
    Mesh::get_num_fields(
            const mtk::EntityRank aEntityRank,
            const moris_index     aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_number_of_real_scalar_fields();
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_number_of_real_scalar_fields();
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Entity not supported in hmr::Mesh::get_num_fields()" );
                return 0;
            }
        }
    }

    //-------------------------------------------------------------------------------

    real&
    Mesh::get_value_of_scalar_field(
            const moris_index     aFieldIndex,
            const mtk::EntityRank aEntityRank,
            const uint            aEntityIndex,
            const moris_index     aIndex )
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_real_scalar_field_data( aFieldIndex )( aEntityIndex );
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_real_scalar_field_coeffs( aFieldIndex )( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Entity not supported in hmr::Mesh::get_value_of_scalar_field()" );
                return mDummyReal;
                break;
            }
        }
    }

    //-------------------------------------------------------------------------------

    const real&
    Mesh::get_value_of_scalar_field(
            const moris_index     aFieldIndex,
            const mtk::EntityRank aEntityRank,
            const uint            aEntityIndex,
            const moris_index     aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_real_scalar_field_data( aFieldIndex )( aEntityIndex );
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_real_scalar_field_coeffs( aFieldIndex )( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Entity not supported in hmr::Mesh::get_value_of_scalar_field()" );
                return mDummyReal;
            }
        }
    }

    //-------------------------------------------------------------------------------

    Matrix< DDRMat >&
    Mesh::get_field(
            const moris_index     aFieldIndex,
            const mtk::EntityRank aEntityRank,
            const moris_index     aIndex )
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_real_scalar_field_data( aFieldIndex );
                break;
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_real_scalar_field_coeffs( aFieldIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Entity not supported in hmr::Mesh::get_field()" );
                return mDummyMatrix;
            }
        }
    }

    //-------------------------------------------------------------------------------

    moris_index
    Mesh::get_field_ind(
            const std::string&    aFieldLabel,
            const mtk::EntityRank aEntityRank ) const
    {
        if ( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::BSPLINE )
        {
            // fixme: this is not a good method. A map would be better
            moris_index aIndex          = gNoIndex;
            moris_index tNumberOfFields = mMesh->get_number_of_real_scalar_fields();
            for ( moris_index k = 0; k < tNumberOfFields; ++k )
            {
                if ( mMesh->get_real_scalar_field_label( k ) == aFieldLabel )
                {
                    aIndex = k;
                    break;
                }
            }
            return aIndex;
        }
        else
        {
            MORIS_ERROR( false, "Entity not supported in hmr::Mesh::get_field_ind()" );
            return gNoIndex;
        }
    }

    //-------------------------------------------------------------------------------

    void
    Mesh::get_sideset_elems_loc_inds_and_ords(
            const std::string&  aSetName,
            Matrix< IndexMat >& aElemIndices,
            Matrix< IndexMat >& aSideOrdinals ) const
    {
        Matrix< DDUMat > tPatternList = mDatabase->create_output_pattern_list();

        if ( mMesh->get_activation_pattern() == tPatternList( 0 ) )    // FIXME loop for more pattern
        {
            // get ref to set
            const Side_Set& tSet = mDatabase->get_output_side_set( aSetName );

            if ( tSet.mElemIdsAndSideOrds.n_rows() > 0 )
            {
                // copy indices into output
                aElemIndices = tSet.mElemIndices;

                // get side id of this set
                uint tSide = tSet.mElemIdsAndSideOrds( 0, 1 );

                // initialize ordinals
                aSideOrdinals.set_size( aElemIndices.length(), 1, tSide );
            }
            else
            {
                // sideset does not exist on this proc
                aElemIndices  = Matrix< IndexMat >( 0, 1 );
                aSideOrdinals = Matrix< IndexMat >( 0, 1 );
            }
        }
        else
        {
            MORIS_ERROR( false, "HMR only generates sidesets for meshes that are linked to the output pattern" );
        }
    }

    //-------------------------------------------------------------------------------

    Vector< std::string >
    Mesh::get_set_names( mtk::EntityRank aSetEntityRank ) const
    {
        if ( aSetEntityRank == mtk::EntityRank::ELEMENT )
        {
            std::string tDummy = "HMR_dummy";

            Vector< std::string > tSetNames( 1, tDummy );

            return tSetNames;
        }
        else if ( aSetEntityRank == mtk::EntityRank::FACE )
        {
            Matrix< DDUMat > tPatternList = mDatabase->create_output_pattern_list();

            if ( mMesh->get_activation_pattern() == tPatternList( 0 ) )
            {
                uint tNumSideSets = mDatabase->get_side_sets().size();

                Vector< std::string > tSetNames( tNumSideSets );

                for ( uint iEntity = 0; iEntity < tNumSideSets; ++iEntity )
                {
                    tSetNames( iEntity ) = mDatabase->get_side_sets()( iEntity ).mInfo.mSideSetName;
                }

                return tSetNames;
            }
            else
            {
                return Vector< std::string >( 0 );
            }
        }
        else if ( aSetEntityRank == mtk::EntityRank::NODE )
        {
            return Vector< std::string >( 0 );
        }
        else
        {
            MORIS_ERROR( false, "Mesh::get_set_names(), only mtk::EntityRank::ELEMENT/FACE is implemented for HMR. Rest can be implemented by you." );
        }

        return Vector< std::string >( 0 );
    }

    //-------------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_set_entity_loc_inds(
            mtk::EntityRank aSetEntityRank,
            std::string     aSetName ) const
    {
        if ( aSetEntityRank == mtk::EntityRank::ELEMENT )
        {
            uint tNumEntities = mMesh->get_number_of_elements();

            Matrix< IndexMat > tOutputEntityInds( tNumEntities, 1 );

            for ( uint iEntity = 0; iEntity < tNumEntities; ++iEntity )
            {
                tOutputEntityInds( iEntity ) = mMesh->get_element( iEntity )->get_index();
            }

            return tOutputEntityInds;
        }

        if ( aSetEntityRank == mtk::EntityRank::FACE )
        {
            // todo: fix this
            Matrix< IndexMat > tSideSetElementInd = mDatabase->get_output_side_set( aSetName ).mElemIndices;

            return tSideSetElementInd;
        }
        else
        {
            MORIS_ERROR( false, "Mesh::get_set_entity_loc_inds(), only mtk::EntityRank::ELEMENT/FACE is implemented for HMR. Rest can be implemented by you." );
        }

        return Matrix< IndexMat >( 0, 0 );
    }

    //-------------------------------------------------------------------------------

    uint
    Mesh::get_level_of_entity_loc_ind(
            const mtk::EntityRank aEntityRank,
            const uint            aEntityIndex,
            const moris_index     aIndex )
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::ELEMENT:
            {
                return mMesh->get_element( aEntityIndex )->get_level();
            }
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_node_by_index( aEntityIndex )->get_level();
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_bspline( aIndex, aEntityIndex )->get_level();
            }
            default:
            {
                MORIS_ERROR( false, "get_level_of_entity_loc_ind: invalid entity rank" );
                return 0;
            }
        }
    }

    //-------------------------------------------------------------------------------

    uint
    Mesh::get_max_level_of_entity(
            const mtk::EntityRank aEntityRank,
            const moris_index     aIndex )
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mMesh->get_max_level();
            }
            case mtk::EntityRank::BSPLINE:
            {
                return mMesh->get_bspline_mesh( aIndex )->get_max_level();
            }
            default:
            {
                MORIS_ERROR( false, "get_level_of_entity_loc_ind: invalid entity rank" );
                return 0;
            }
        }
    }

    //-------------------------------------------------------------------------------

    mtk::Cell&
    Mesh::get_mtk_cell( moris_index aElementIndex )
    {
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            return *mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            return *mMesh->get_element( aElementIndex );
        }
    }

    //-------------------------------------------------------------------------------

    mtk::Cell&
    Mesh::get_writable_mtk_cell( moris_index aElementIndex )
    {
        return this->get_mtk_cell( aElementIndex );
    }

    //-------------------------------------------------------------------------------

    mtk::Cell const &
    Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        if ( mDatabase->get_parameters()->use_number_aura() and    //
                mDatabase->get_parameters()->is_output_mesh( mMesh->get_index() ) )
        {
            return *mMesh->get_element_including_aura( aElementIndex );
        }
        else
        {
            return *mMesh->get_element( aElementIndex );
        }
    }

    //-------------------------------------------------------------------------------

    Matrix< IndexMat >
    Mesh::get_inds_of_active_elements_connected_to_basis( const Basis* aBasis ) const
    {
        // ask basis for number of elements
        luint tNumberOfElements = aBasis->get_element_counter();

        // get activation pattern from this mesh
        uint tPattern = mMesh->get_activation_pattern();

        // counter for active elements
        luint tCount = 0;

        // count active elements
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            aBasis->get_element( e )->get_background_element()->    //
                    get_number_of_active_descendants( tPattern, tCount );
        }

        // allocate cell with connected elements
        Vector< const Background_Element_Base* > tBackElements( tCount, nullptr );

        // reset counter
        tCount = 0;

        // collect background elements
        for ( luint e = 0; e < tNumberOfElements; ++e )
        {
            aBasis->get_element( e )->get_background_element()->    //
                    collect_active_descendants(
                            tPattern,
                            tBackElements,
                            tCount );
        }

        // allocate output matrix
        Matrix< IndexMat > aElementIndices( tCount, 1 );

        // loop over all background elements
        for ( luint e = 0; e < tCount; ++e )
        {
            // grab element by memory index and copy moris index into matrix
            aElementIndices( e ) = mMesh->get_element_by_memory_index( tBackElements( e )->get_memory_index() )->get_index();
        }

        return aElementIndices;
    }

    //-------------------------------------------------------------------------------

    void
    Mesh::setup_glb_to_local_maps()
    {
        // report on this operation
        MORIS_LOG_INFO( "Setting up global to local maps on HMR mesh" );

        // Initialize global to local map
        mEntityGlobalToLocalMap = Vector< std::unordered_map< moris_id, moris_index > >( 4 + mMesh->get_number_of_bspline_meshes() );

        // count number of mesh entities
        uint tCounter = 0;

        setup_entity_global_to_local_map( mtk::EntityRank::NODE, tCounter );
        setup_entity_global_to_local_map( mtk::EntityRank::EDGE, tCounter );
        setup_entity_global_to_local_map( mtk::EntityRank::FACE, tCounter );
        setup_entity_global_to_local_map( mtk::EntityRank::ELEMENT, tCounter );

        for ( uint iBspMesh = 0; iBspMesh < mMesh->get_number_of_bspline_meshes(); iBspMesh++ )
        {
            if ( mMesh->get_bspline_mesh( iBspMesh ) != nullptr )
            {
                setup_entity_global_to_local_map( mtk::EntityRank::BSPLINE, tCounter, iBspMesh );
            }
            else
            {
                // trivial case when all t-matrices are 1
                setup_entity_global_to_local_map( mtk::EntityRank::NODE, tCounter );
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Mesh::setup_entity_global_to_local_map(
            mtk::EntityRank   aEntityRank,
            uint&             aCounter,
            const moris_index aIndex )
    {
        uint tNumEntities = this->get_num_entities( aEntityRank, aIndex );

        moris_index tCount = 0;

        for ( uint i = 0; i < tNumEntities; i++ )
        {
            moris_id tEntityId = this->get_glb_entity_id_from_entity_loc_index(
                    i,
                    aEntityRank,
                    aIndex );

            // MORIS_ASSERT(mEntityGlobalToLocalMap(aCounter).find(tEntityId) == mEntityGlobalToLocalMap(aCounter).end(),"ID already exists.");

            MORIS_ASSERT( tEntityId >= 0, "EntityID received is smaller than 0." );

            mEntityGlobalToLocalMap( aCounter )[ tEntityId ] = tCount;

            tCount++;
        }

        aCounter++;
    }

    //-------------------------------------------------------------------------------

     mtk::CellTopology
    Mesh::get_blockset_topology( const std::string& aSetName )
    {
        uint tNumberOfDimensions = mDatabase->get_number_of_dimensions();

        uint tOrder = this->get_order();

         mtk::CellTopology tCellTopology = mtk::CellTopology::UNDEFINED;

        switch ( tNumberOfDimensions )
        {
            case 1:
            {
                MORIS_ERROR( false, "1D not implemented" );
                break;
            }
            case 2:
            {
                switch ( tOrder )
                {
                    case 1:
                    {
                        tCellTopology = mtk::CellTopology::QUAD4;
                        break;
                    }
                    case 2:
                    {
                        tCellTopology = mtk::CellTopology::QUAD9;
                        break;
                    }
                    case 3:
                    {
                        tCellTopology = mtk::CellTopology::QUAD16;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Order not implemented" );
                        break;
                    }
                }
                break;
            }
            case 3:
            {
                switch ( tOrder )
                {
                    case 1:
                    {
                        tCellTopology = mtk::CellTopology::HEX8;
                        break;
                    }
                    case 2:
                    {
                        tCellTopology = mtk::CellTopology::HEX27;
                        break;
                    }
                    case 3:
                    {
                        tCellTopology = mtk::CellTopology::HEX64;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Order not implemented" );
                        break;
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Number of dimensions not implemented" );
                break;
            }
        }
        return tCellTopology;
    }

    //-------------------------------------------------------------------------------

    mtk::CellShape
    Mesh::get_IG_blockset_shape( const std::string& aSetName )
    {
        return mtk::CellShape::RECTANGULAR;
    }

    //-------------------------------------------------------------------------------

    mtk::CellShape
    Mesh::get_IP_blockset_shape( const std::string& aSetName )
    {
        return mtk::CellShape::RECTANGULAR;
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
