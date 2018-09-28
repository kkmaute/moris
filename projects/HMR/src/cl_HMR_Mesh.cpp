#include <string>
#include "cl_MTK_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Mesh.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR.hpp" //HMR/src
#include "fn_sort.hpp"
#include "fn_unique.hpp"
#include "fn_print.hpp"

namespace moris
{
    namespace hmr
    {
//-----------------------------------------------------------------------------

        Mesh::Mesh( HMR & aHMR,
                const uint & aOrder,
                const uint & aActivationPattern ) : mHMR( aHMR )
        {

            // get number of meshes
            uint tNumberOfMeshes = mHMR.get_number_of_lagrange_meshes();


            // find correct block
            for( uint k=0; k<tNumberOfMeshes; ++k )
            {
                auto tMesh = mHMR.get_lagrange_mesh_by_index( k );

                // test if mesh uses active pattern
                if ( tMesh->get_activation_pattern() == aActivationPattern &&
                     tMesh->get_order() == aOrder )
                {
                    mMesh = tMesh;
                    mBlock = new hmr::Block( tMesh, k );
                    break;
                }
            }


        }

//-----------------------------------------------------------------------------

        Mesh::~Mesh()
        {
            delete mBlock;
        }


//-----------------------------------------------------------------------------

        Matrix< IdMat >
        Mesh::get_communication_table() const
        {
            return mHMR.get_communication_table();
        }

//-----------------------------------------------------------------------------

        mtk::Field *
        Mesh::create_field( const std::string & aLabel )
        {
            // fixme: rethink the concept of multiple blocks on HMR

            // create field
            return new Field( aLabel, mBlock );
        }

//-----------------------------------------------------------------------------
// Blocks
//-----------------------------------------------------------------------------

        uint
        Mesh::get_number_of_blocks() const
        {
            return 1;
        }

//-----------------------------------------------------------------------------

        Block *
        Mesh::get_block_by_index( const moris_index & aIndex )
        {
            MORIS_ASSERT( aIndex == 0, "an HMR mesh has only one block");
            return mBlock;

        }

//-----------------------------------------------------------------------------

        const Block *
        Mesh::get_block_by_index( const moris_index & aIndex ) const
        {
            MORIS_ASSERT( aIndex == 0, "an HMR mesh has only one block");
            return mBlock;

        }

//-----------------------------------------------------------------------------
// MTK
//-----------------------------------------------------------------------------

        uint
        Mesh::get_spatial_dim() const
        {
            return mHMR.get_parameters()->get_number_of_dimensions();
        }

//-----------------------------------------------------------------------------

        uint
        Mesh::get_num_entities(
                enum EntityRank aEntityRank ) const
        {
            switch ( aEntityRank )
            {
                case( EntityRank::NODE ) :
                {
                    return this->get_num_nodes();
                    break;
                }
                case( EntityRank::EDGE ) :
                {
                     return this->get_num_edges();
                     break;
                }
                case( EntityRank::FACE ) :
                {
                     return this->get_num_faces();
                     break;
                }
                case( EntityRank::ELEMENT ) :
                {
                    return this->get_num_elems();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown entity rank");
                    return 0;
                    break;
                }
            }
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
            return mMesh->get_number_of_edges();
        }

//-----------------------------------------------------------------------------

        uint
        Mesh::get_num_faces() const
        {
            return mMesh->get_number_of_facets();
        }

//-----------------------------------------------------------------------------

        uint
        Mesh::get_num_elems() const
        {
            return mMesh->get_number_of_elements();
        }

//-----------------------------------------------------------------------------

        Matrix<IndexMat>
        Mesh::get_entity_connected_to_entity_loc_inds(
                           moris_index     aEntityIndex,
                           enum EntityRank aInputEntityRank,
                           enum EntityRank aOutputEntityRank) const
        {
            switch ( aOutputEntityRank )
            {
                case( EntityRank::NODE ) :
                {
                    switch ( aInputEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide node to node connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::EDGE ) :
                        {
                            return this->get_nodes_connected_to_edge_loc_inds( aEntityIndex );
                            break;
                        }
                        case( EntityRank::FACE ) :
                        {
                            return this->get_nodes_connected_to_face_loc_inds( aEntityIndex );
                            break;
                        }
                        case( EntityRank::ELEMENT ) :
                        {
                            return this->get_nodes_connected_to_element_loc_inds( aEntityIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide the requested connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                    }
                    break;
                } // end output rank node
                case( EntityRank::EDGE ) :
                {
                    switch( aInputEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            return get_edges_connected_to_node_loc_inds( aEntityIndex );
                            break;
                        }
                        case( EntityRank::EDGE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide edge to edge connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::FACE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide edge to face connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::ELEMENT ) :
                        {
                            return get_edges_connected_to_element_loc_inds ( aEntityIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide the requested connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                    }
                    break;
                } // end output rank edge
                case( EntityRank::FACE ) :
                {
                    switch( aInputEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            return get_faces_connected_to_node_loc_inds( aEntityIndex );
                            break;
                        }
                        case( EntityRank::EDGE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide face to edge connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::FACE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide face to face connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::ELEMENT ) :
                        {
                            return get_faces_connected_to_element_loc_inds( aEntityIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide the requested connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                    }
                    break;
                } // end output rank face
                case( EntityRank::ELEMENT ) :
                {
                    switch( aInputEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            return get_elements_connected_to_node_loc_inds ( aEntityIndex );
                            break;
                        }
                        case( EntityRank::EDGE ) :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide element to edge connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                        case( EntityRank::FACE ) :
                        {
                            return get_elements_connected_to_face_loc_inds( aEntityIndex );
                            break;
                        }
                        case( EntityRank::ELEMENT ) :
                        {
                            return get_elements_connected_to_element_loc_inds( aEntityIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false,
                                    "HMR does not provide the requested connectivity" );
                            return Matrix<IndexMat>( 0, 0 );
                            break;
                        }
                    }
                    break;
                } // end output rank element
                default :
                {
                    MORIS_ERROR( false,
                            "HMR does not provide the requested connectivity" );
                    return Matrix<IndexMat>( 0, 0 );
                    break;
                }
            }
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
            // get pointer to element
            Element * tElement = mMesh->get_element( aElementIndex );

            // get number of nodes
            uint tNumberOfNodes = tElement->get_number_of_vertices();

            // allocate output
            Matrix< IndexMat > aIndices( tNumberOfNodes, 1 );

            // populate output
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                aIndices( k ) = tElement->get_basis( k )->get_index();
            }

            return aIndices;
        }

//-----------------------------------------------------------------------------

        Matrix < IndexMat >
        Mesh::get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const
        {
            // get pointer to basis
            Basis * tBasis = mMesh->get_node_by_index( aNodeIndex );

            uint tNumberOfEdges = tBasis->get_edge_counter();

            Matrix < IndexMat > tEdgeIndex( tNumberOfEdges , 1 );

            for( uint k=0; k<tNumberOfEdges; ++k )
            {
                tEdgeIndex( k ) = tBasis->get_edge( k )->get_index();
            }

            Matrix < IndexMat > aEdgeIndex;
            sort( tEdgeIndex, aEdgeIndex );
            return aEdgeIndex;
        }

//-----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh::get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const
        {
            // get pointer to element
            Element * tElement = mMesh->get_element( aElementIndex );

            // will be 0 for 2D, 12 for 3D
            uint tNumberOfEdges = tElement->get_background_element()->get_number_of_edges();

            // allocate output
            Matrix< IndexMat > aIndices( tNumberOfEdges, 1 );

            // populate output
            for( uint e=0; e<tNumberOfEdges; ++e )
            {
                aIndices( e ) = tElement->get_hmr_edge( e )->get_index();
            }

            return aIndices;
        }

//-----------------------------------------------------------------------------

        Matrix < IndexMat >
        Mesh::get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const
        {
            // get pointer to basis
            Basis * tBasis = mMesh->get_node_by_index( aNodeIndex );

            uint tNumberOfFacets = tBasis->get_facet_counter();

            Matrix < IndexMat > tFaceIndex( tNumberOfFacets , 1 );

            for( uint k=0; k<tNumberOfFacets; ++k )
            {
                tFaceIndex( k ) = tBasis->get_facet( k )->get_index();
            }

            Matrix < IndexMat > aFaceIndex;
            sort( tFaceIndex, aFaceIndex );
            return aFaceIndex;
        }

//-----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh::get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const
        {
            // get pointer to element
            Element * tElement = mMesh->get_element( aElementIndex );

            uint tNumberOfFaces
            = tElement->get_background_element()->get_number_of_facets();

            Matrix< IndexMat > aIndices( tNumberOfFaces, 1 );

            for( uint f=0; f<tNumberOfFaces; ++f )
            {
                aIndices( f ) = tElement->get_hmr_facet( f )->get_index();
            }

            return aIndices;
        }

//-----------------------------------------------------------------------------

        Matrix < IndexMat >
        Mesh::get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const
        {
            // collect memory indices of active elements
            Matrix< DDLUMat> tMemoryIndices;
            this->collect_memory_indices_of_active_elements_connected_to_node(
                    aNodeIndex, tMemoryIndices );

            Matrix< IndexMat > aIndices;
            this->get_element_indices_from_memory_indices(
                    tMemoryIndices, aIndices );

            return aIndices;
        }

//-----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh::get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const
        {

            // get pointer to facet
            Facet * tFacet = mMesh->get_facet( aFaceIndex );

            Matrix< IndexMat > aIndices;

            if( tFacet->get_slave() == NULL )
            {
                if( tFacet->get_hmr_master()->is_active() )
                {
                    aIndices.set_size( 1, 1 );
                    aIndices( 0 ) = tFacet->get_hmr_master()->get_index();
                }
            }
            else
            {
                Element * tMaster = tFacet->get_hmr_master();
                Element * tSlave  = tFacet->get_hmr_slave();

                if(  tMaster->is_active() && tSlave->is_active() )
                {
                    aIndices.set_size( 2, 1 );
                    aIndices( 0 ) = tMaster->get_index();
                    aIndices( 1 ) = tSlave->get_index();
                }
                else if ( tMaster->is_active() )
                {
                    aIndices.set_size( 1, 1 );
                    aIndices( 0 ) = tMaster->get_index();
                }
                else if ( tSlave->is_active() )
                {
                    aIndices.set_size( 1, 1 );
                    aIndices( 0 ) = tSlave->get_index();
                }
            }

            return aIndices;
        }

//-----------------------------------------------------------------------------

        Matrix< IndexMat >
        Mesh::get_elements_connected_to_element_loc_inds( moris_index aElementIndex ) const
        {

            // collect memory indices of active neighbors
            Matrix< DDLUMat> tMemoryIndices;
            luint tNumberOfNeighbors;

            this->collect_memory_indices_of_active_element_neighbors(
                    aElementIndex, tMemoryIndices, tNumberOfNeighbors );

            // reserve memory for output matrix
            Matrix< IndexMat > tIndices( tNumberOfNeighbors, 1 );

            // copy indices from pointers
            for( luint k=0; k<tNumberOfNeighbors; ++k )
            {
                tIndices( k )
                       = mMesh->get_element_by_memory_index( tMemoryIndices( k ) )->get_index();
            }

            Matrix< IndexMat > aIndices;

            // make result unique
            unique( tIndices, aIndices );
            return aIndices;
        }

//-----------------------------------------------------------------------------
// private
//-----------------------------------------------------------------------------

        void
        Mesh::get_element_indices_from_memory_indices(
                const Matrix< DDLUMat>      & aMemoryIndices,
                      Matrix< IndexMat >    & aIndices ) const
        {
            // get length of vector
            uint tNumberOfElements = aMemoryIndices.length();

            // allocate memory
            aIndices.set_size( tNumberOfElements, 1);

            // copy indices from pointers
            for( uint k=0; k<tNumberOfElements; ++k )
            {
                aIndices( k ) = mMesh->get_element_by_memory_index(
                        aMemoryIndices( k ) )->get_index();
            }

        }

//-----------------------------------------------------------------------------

        void
        Mesh::collect_memory_indices_of_active_element_neighbors(
                            const moris_index  aElementIndex,
                            Matrix< DDLUMat> & aMemoryIndices,
                            luint            & aCounter) const
        {
            // get active index of this mesh
            uint tPattern = mMesh->get_activation_pattern();

            // get pointer to element
            Element * tElement = mMesh->get_element( aElementIndex );

            // get pointer to background element
            Background_Element_Base * tBackElement = tElement->get_background_element();

            // count elements that are not padding
            aCounter = 0;

            uint tNumberOfFacets = tBackElement->get_number_of_facets();

            for( uint k=0; k<tNumberOfFacets; ++k )
            {
                // get pointer to neighbor
                Background_Element_Base * tNeighbor = tBackElement->get_neighbor( k );

                if( ! tNeighbor->is_padding() )
                {
                    if( tNeighbor->is_refined( tPattern ) )
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

            // allocate matrix with memory indices
            aMemoryIndices.set_size( aCounter, 1 );

            // reset counter
            aCounter = 0;

            // loop over all neighbors
            for( uint k=0; k<tNumberOfFacets; ++k )
            {
                // get pointer to neighbor
                Background_Element_Base * tNeighbor = tBackElement->get_neighbor( k );

                if( ! tNeighbor->is_padding() )
                {
                    if( tNeighbor->is_refined( tPattern ) )
                    {
                        tNeighbor->collect_active_descendants_by_memory_index( tPattern, aMemoryIndices, aCounter, k );
                    }
                    else
                    {
                        while ( ! tNeighbor->is_active( tPattern ) )
                        {
                            tNeighbor = tNeighbor->get_parent();
                        }
                        aMemoryIndices( aCounter++ ) = tNeighbor->get_memory_index();
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Mesh::collect_memory_indices_of_active_elements_connected_to_node(
                           const moris_index  aNodeIndex,
                           Matrix< DDLUMat> & aMemoryIndices ) const
        {

            // get active index of this mesh
            uint tPattern = mMesh->get_activation_pattern();

            // get pointer to node
            Basis * tNode = mMesh->get_node_by_index( aNodeIndex );

            // get number of elements
            luint tElementCounter = tNode->get_element_counter();

            // element counter with active elements
            luint tCount = 0;

            for( uint k=0; k<tElementCounter; ++k )
            {

                if( tNode->get_element( k )->get_background_element()->is_active( tPattern ) )
                {
                    // increment counter
                    ++tCount;
                }
            }

            // allocate matrix with memory indices
            aMemoryIndices.set_size( tCount, 1 );

            // reset counter
            tCount = 0;

            for( uint k=0; k<tElementCounter; ++k )
            {
                // get pointer to element
                Background_Element_Base * tBackElement
                    = tNode->get_element( k )->get_background_element();

                if( tBackElement->is_active( tPattern ) )
                {
                    // increment counter
                    aMemoryIndices( tCount++ ) = tBackElement->get_memory_index( );
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
                moris_index     aEntityIndex,
                enum EntityRank aEntityRank) const
        {
            switch ( aEntityRank )
            {
                case( EntityRank::NODE ) :
                {
                    return mMesh->get_node_by_index( aEntityIndex )->get_id();
                    break;
                }
                case( EntityRank::EDGE ) :
                {
                    return mMesh->get_edge( aEntityIndex )->get_id();
                    break;
                }
                case( EntityRank::FACE ) :
                {
                    return mMesh->get_facet( aEntityIndex )->get_id();
                    break;
                }
                case( EntityRank::ELEMENT ) :
                {
                    return mMesh->get_element( aEntityIndex )->get_id();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown entity rank");
                    return 0;
                    break;
                }
            }
        }

//-----------------------------------------------------------------------------

        moris_id
        Mesh::get_entity_owner(  moris_index     aEntityIndex,
                         enum EntityRank aEntityRank ) const
        {
            switch ( aEntityRank )
            {
                case( EntityRank::NODE ) :
                {
                    return mMesh->get_node_by_index( aEntityIndex )->get_owner();
                    break;
                }
                case( EntityRank::EDGE ) :
                {
                    return mMesh->get_edge( aEntityIndex )->get_owner();
                    break;
                }
                case( EntityRank::FACE ) :
                {
                    return mMesh->get_facet( aEntityIndex )->get_owner();
                    break;
                }
                case( EntityRank::ELEMENT ) :
                {
                    return mMesh->get_element( aEntityIndex )->get_owner();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown entity rank");
                    return 0;
                    break;
                }
            }
        }


//-----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
