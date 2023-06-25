/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Mesh_Base.cpp
 *
 */

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

#include <cstdio>
#include <fstream>

#include "cl_HMR_Background_Edge.hpp"
#include "cl_HMR_Background_Facet.hpp"
#include "cl_HMR_Facet.hpp"
#include "cl_HMR_T_Matrix.hpp"
#include "cl_HMR_T_Matrix_Advanced.hpp"
#include "HMR_Tools.hpp"
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_Map.hpp"
#include "cl_HMR_Factory.hpp"

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "HDF5_Tools.hpp"

#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_equal_to.hpp"

namespace moris::hmr
{
    //------------------------------------------------------------------------------
    //   public:
    //------------------------------------------------------------------------------

    Lagrange_Mesh_Base::Lagrange_Mesh_Base (
            Parameters const             * aParameters,
            Background_Mesh_Base         * aBackgroundMesh,
            Cell< BSpline_Mesh_Base *  > & aBSplineMeshes,
            uint                           aOrder,
            uint                           aActivationPattern )
    : Mesh_Base( aParameters,
            aBackgroundMesh,
            aOrder,
            aActivationPattern ),
            mBSplineMeshes( aBSplineMeshes ) // initialize member cell of associated B-Spline meshes by copying aBSplineMeshes into it
    {
        mNumBSplineMeshes = mBSplineMeshes.size();

        this->check_if_bspline_mesh_is_trivial_interpolation();
        this->reset_fields();
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::update_mesh()
    {
        // start timer
        tic tTimer;

        // activate pattern on background mesh
        this->select_activation_pattern();

        // tidy up memory
        this->delete_pointers();
        this->delete_facets();

        // create Lagrange Elements from Background Elements
        this->create_elements();

        // create nodes
        this->create_nodes();

        // update list of used nodes
        this->update_node_list();

        // update element indices
        this->update_element_indices();

        // link elements to B-Spline meshes
        //if(  mBSplineMeshes.size() > 0 )
        //{
        //this->link_twins();
        //}

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // print output
        MORIS_LOG_INFO( "%s Created Lagrange mesh of order %u on pattern %u.",
                proc_string().c_str(),
                ( unsigned int ) mOrder,
                ( unsigned int ) mActivationPattern);

        MORIS_LOG_INFO( "Mesh has %lu active and refined elements and %lu nodes.",
                sum_all( ( long unsigned int ) this->get_number_of_elements() ),
                sum_all( ( long unsigned int ) this->get_number_of_nodes_on_proc() ) );

        MORIS_LOG_INFO( "Creation took %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    // ----------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_to_file( const std::string & aFilePath )
    {
        // get the file extension
        auto tFileExt = aFilePath.substr(aFilePath.find_last_of(".")+1,
                aFilePath.length() );

        // guess routine from extension
        if (tFileExt == "vtk")
        {
            this->save_to_vtk( aFilePath );
        }
        else if(tFileExt == "msh")
        {
            this->save_to_gmsh( aFilePath );
        }
        else
        {
            MORIS_ERROR( false,
                    "Wrong file type passed to Lagrange mesh.\n If you want to save an exodus file, create an MTK object first." );
        }
    }

    //------------------------------------------------------------------------------

    uint Lagrange_Mesh_Base::create_real_scalar_field_data( const std::string     & aLabel,
            const enum EntityRank   aEntityRank )
    {
        // get index for output
        uint aIndex = mRealScalarFieldData.size();

        MORIS_ERROR( mRealScalarFieldData.size() == mRealScalarFieldLabels.size() ,
                "Sizes of Field labels and Data container does not match " );

        // initialize empty matrix. It is populated later
        Matrix< DDRMat > tEmpty;

        mRealScalarFieldLabels       .push_back( aLabel );
        mRealScalarFieldData         .push_back( tEmpty );
        mRealScalarFieldBSplineCoeffs.push_back( tEmpty );
        mRealScalarFieldBSplineOrders.push_back( 0 );
        mRealScalarFieldRanks        .push_back( aEntityRank );

        return aIndex;
    }

    //------------------------------------------------------------------------------
    //   protected:
    // -----------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_nodes_on_higher_levels()
    {
        // get max level of mest
        uint tMaxLevel = mBackgroundMesh->get_max_level();

        // loop over all levels
        for ( uint l=0; l<tMaxLevel; ++l )
        {
            // get all elements from this level
            Cell< Background_Element_Base* > tElements;

            mBackgroundMesh->collect_elements_on_level_including_aura( l, tElements );

            // loop over all elements on this level
            for( auto tElement : tElements )
            {
                // test if this element has children and is not padding and is refined
                if ( tElement->has_children() and ! tElement->is_padding() && tElement->is_refined( mActivationPattern ) )
                {
                    // calculate nodes of children
                    mNumberOfAllBasis += mAllElementsOnProc( tElement->get_memory_index() )->create_basis_for_children( mAllElementsOnProc );
                }
            }
        }
    }

    //------------------------------------------------------------------------------
    //   private:
    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_nodes()
    {
        // nodes on first level are created separately
        this->create_nodes_on_level_zero();

        // it is easier to create all zero level nodes first
        // and remove obsolete ones. Don't worry, there are not that many.
        this->delete_nodes_on_padding_elements();

        // now we do the nodes on higher levels
        this->create_nodes_on_higher_levels();

        // fill cell with all remaining nodes
        this->collect_nodes();

        // calculate system wide IDs (helpful for debugging)
        this->calculate_node_ids();

        // find out node to element connectivity
        this->determine_elements_connected_to_basis();

        // and determine the node ownership to be the smallest proc ID
        // of all elements connected to this node
        this->guess_basis_ownership();

        // Make sure that node ownership is correct. Correct otherwise.
        this->confirm_basis_ownership();

        // calculate node coordinates with respect to user defined offset
        this->calculate_node_coordinates();

        // create node numbers
        this->calculate_node_indices();

        if( mParameters->use_number_aura() and
                mParameters->get_number_of_dimensions() == 3 and
                par_size() >= 4 and
                mParameters->is_output_mesh( mMeshIndex ) )
        {
            this->communicate_missed_node_indices();
        }

        //#ifdef MORIS_HAVE_DEBUG
        //            MORIS_LOG_WARNING("Sanity check for vertex basis Ids and ownership will be performed. This might slow down the execution significantly. \n");
        //            this->sanity_check_for_ids_and_ownership();
        //#endif
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::delete_nodes_on_padding_elements()
    {
        // get number of elements on coarsest level
        auto tNumberOfElements = mAllCoarsestElementsOnProc.size();

        // counter for nodes
        luint tCount = 0;

        // loop over all elements on coarsest level
        for( luint e = 0; e < tNumberOfElements; ++e)
        {
            // get pointer to Lagrange element
            Element* tElement = mAllCoarsestElementsOnProc( e );

            // test if element is not padding
            if ( ! tElement->is_padding() )
            {
                // loop over all nodes of element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to node
                    Basis* tNode = tElement->get_basis( k );

                    if ( ! tNode->is_flagged() )
                    {
                        // set basis as active
                        tNode->flag();

                        // increment counter
                        ++tCount;
                    }
                }
            }
        }

        // ask background mesh for number of elements per direction
        Matrix< DDLUMat > tNumberOfElementsPerDirection = mBackgroundMesh->get_number_of_elements_per_direction_on_proc();

        // assign Cell for nodes to be deleted
        Cell< Basis* > tNodes( mNumberOfAllBasis-tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all elements on coarsest level
        for( luint e = 0; e < tNumberOfElements; ++e)
        {
            // get pointer to Lagrange element
            Element* tElement = mAllCoarsestElementsOnProc( e );

            // test if element is padding
            if ( tElement->is_padding() )
            {
                // loop over all nodes of element
                for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to node
                    Basis* tNode = tElement->get_basis( k );

                    // test if node exists
                    if ( tNode != nullptr )
                    {
                        // test if node is not flagged
                        if ( ! tNode->is_flagged() )
                        {
                            // flag node
                            tNode->flag();

                            // copy node to delete list
                            tNodes( tCount++ ) = tNode;
                        }
                    }
                }
            }
        }

        // delete nodes
        for ( uint k = 0; k < tCount; ++k )
        {
            delete tNodes( k );

            // decrement number of nodes
            --mNumberOfAllBasis;
        }

        // tidy up: unflag all remaining nodes
        // loop over all elements on coarsest level
        for( luint e = 0; e < tNumberOfElements; ++e)
        {
            // get pointer to Lagrange element
            Element* tElement = mAllCoarsestElementsOnProc( e );

            // test if element is not padding
            if ( ! tElement->is_padding() )
            {
                // loop over all nodes of element
                for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to node
                    Basis* tNode = tElement->get_basis( k );

                    // unset flag
                    tNode->unflag();
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::collect_nodes()
    {
        // clear node list
        mAllBasisOnProc.clear();

        // reserve size
        mAllBasisOnProc.resize( mNumberOfAllBasis, nullptr );

        // initialize counter
        luint tCount = 0;

        // get number of active elements on proc
        luint tNumberOfElements = mBackgroundMesh ->get_number_of_active_elements_on_proc_including_aura();

        // get rank
        moris_id tMyRank = par_rank();

        // reset element counter
        mNumberOfElements = 0;
        mNumberOfElementsIncludingAura = 0;

        // loop over all active elements on proc
        for ( luint e = 0; e < tNumberOfElements; ++e )
        {
            // get pointer to background element
            Background_Element_Base* tBackElement = mBackgroundMesh ->get_element_from_proc_domain_including_aura( e );

            // get pointer to Lagrange element
            Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

            if ( ! tBackElement->is_deactive( mActivationPattern )  )
            {
                // flag nodes that are used by this proc
                if ( tBackElement->get_owner() == tMyRank )
                {
                    for ( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                    {
                        tElement->get_basis( k )->use();
                    }

                    // increment element counter
                    ++mNumberOfElements;
                }

                // flag nodes that are owned and shared on this proc. this includes the aura
                for ( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                {
                    tElement->get_basis( k )->use_owned_and_shared();
                }

                // increment element counter including aura
                ++mNumberOfElementsIncludingAura;

                // loop over all nodes of this element
                for ( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to node
                    Basis* tNode = tElement->get_basis( k );

                    // test if node is flagged
                    if ( ! tNode->is_flagged() )
                    {
                        // set index in memory
                        tNode->set_memory_index( tCount );

                        // add node to list
                        mAllBasisOnProc( tCount ++ ) = tNode;

                        // flag node
                        tNode->flag();
                    }
                }
            }
        }

        // make sure that number of nodes is correct
        MORIS_ERROR( tCount == mNumberOfAllBasis, "Number of Nodes does not match." );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::calculate_node_ids()
    {
        switch ( mParameters->get_number_of_dimensions() )
        {
            case( 1 ):
                        {
                for( auto tNode  : mAllBasisOnProc )
                {
                    // get ij position of node
                    const luint * tI = tNode->get_ijk();

                    // calculate ID and write to node
                    tNode->set_domain_id( this->calculate_node_id( tNode->get_level(),
                            tI[0] ) );
                }

                break;
                        }
            case( 2 ):
                        {
                for( auto tNode  : mAllBasisOnProc )
                {
                    // get ij position of node
                    const luint * tIJ = tNode->get_ijk();

                    // calculate ID and write to node
                    tNode->set_domain_id( this->calculate_node_id( tNode->get_level(),
                            tIJ[0],
                            tIJ[1]) );
                }

                break;
                        }
            case( 3 ):
                        {
                // 3D case
                for( auto tNode  : mAllBasisOnProc )
                {
                    // get ij position of node
                    const luint * tIJK = tNode->get_ijk();

                    // calculate ID and write to node
                    tNode->set_domain_id( this->calculate_node_id( tNode->get_level(),
                            tIJK[0],
                            tIJK[1],
                            tIJK[2]) );
                }

                break;
                        }
            default:
            {
                MORIS_ERROR( false, "Lagrange_Mesh: Invalid number of dimensions");
                break;
            }
        }
    }

    //------------------------------------------------------------------------------
    /**
     * calculates system wide unique node indices for MTK
     */
    void Lagrange_Mesh_Base::calculate_node_indices()
    {
        // reset counters
        mNumberOfUsedAndOwnedNodes = 0;
        mNumberOfUsedNodes = 0;

        moris_id tNumberOfProcs = par_size();

        if( tNumberOfProcs == 1 ) // serial mode
        {
            moris_id tCount = 0 ;
            for( auto tNode : mAllBasisOnProc )
            {
                // test if node is used by current setup
                if ( tNode->is_used() )
                {
                    // in serial local index and domain index are identical
                    tNode->set_domain_index( mNumberOfUsedAndOwnedNodes++ );
                    tNode->set_local_index( mNumberOfUsedNodes++ );
                    tCount++;
                }
            }
            mMaxNodeDomainIndex = tCount;
        }
        else // parallel mode
        {
            // STEP 1: label the nodes that I own

            // get my rank
            moris_id tMyRank = par_rank();

            if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
            {
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if node is used by current setup
                    if ( tNode->is_used() )
                    {
                        // test if node is owned
                        if ( tNode->get_owner() == tMyRank )
                        {
                            tNode->set_domain_index( mNumberOfUsedAndOwnedNodes++ );
                        }
                    }
                    if ( tNode->is_use_owned_and_shared() )
                    {
                        // set local index of node including aura
                        tNode->set_local_index( mNumberOfUsedNodes++ );
                    }

                    // make sure that this basis is not flagged
                    tNode->unflag();
                }
            }
            else
            {
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if node is used by current setup
                    if ( tNode->is_used() )
                    {
                        // test if node is owned
                        if ( tNode->get_owner() == tMyRank )
                        {
                            tNode->set_domain_index( mNumberOfUsedAndOwnedNodes++ );
                        }

                        // set local index of node
                        tNode->set_local_index( mNumberOfUsedNodes++ );
                    }

                    // make sure that this basis is not flagged
                    tNode->unflag();
                }
            }

            // STEP 2: Add offset to IDs which I own

            // communicate number of owned nodes with other procs
            Matrix< DDLUMat > tNodesOwnedPerProc;
            comm_gather_and_broadcast( mNumberOfUsedAndOwnedNodes, tNodesOwnedPerProc );

            // get proc neighbors from background mesh
            auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            // calculate node offset table
            Matrix< DDLUMat > tNodeOffset( tNumberOfProcs, 1, 0 );
            for( moris_id p = 1; p < tNumberOfProcs; ++p )
            {
                tNodeOffset( p ) = tNodeOffset( p-1 ) + tNodesOwnedPerProc( p-1 );
            }

            // remember for MTK output
            mMaxNodeDomainIndex = tNodeOffset( tNumberOfProcs-1 ) + tNodesOwnedPerProc( tNumberOfProcs-1 );

            // get my offset
            luint tMyOffset = tNodeOffset( tMyRank );

            // loop over all nodes on proc
            for( auto tNode : mAllBasisOnProc )
            {
                // test if the is used and node belongs to me
                if ( tNode->is_used() )
                {
                    if ( tNode->get_owner() == tMyRank )
                    {
                        // set global node index
                        tNode->set_domain_index( tNode->get_hmr_index() + tMyOffset );
                    }
                }
            }

            // STEP 3: Global node indices of used and owned nodes
            // must be communicated to the other procs

            // get number of proc neighbors
            uint tNumberOfProcNeighbors = mBackgroundMesh->get_number_of_proc_neighbors();

            // create cell of matrices to send
            Matrix< DDLUMat > tEmpty;
            Cell< Matrix< DDLUMat > > tSendID( tNumberOfProcNeighbors, tEmpty );

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                moris_id tNeighbor = tProcNeighbors( p );

                if ( tNeighbor < tNumberOfProcs && tNeighbor != tMyRank )
                {
                    // cell with basis in aura
                    Cell< Basis* > tNodes;

                    // collect nodes within aura
                    this->collect_basis_from_aura( p, 0, tNodes );

                    // count nodes that belong to neighbor
                    uint tCount = 0;

                    for( Basis* tNode : tNodes )
                    {
                        // test if node belongs to neighbor
                        if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_use_owned_and_shared() )
                            {
                                // increment counter
                                ++tCount;
                            }
                        }
                        else
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_used() )
                            {
                                // increment counter
                                ++tCount;
                            }
                        }
                    }

                    // set matrix length
                    tSendID( p ).set_size( tCount, 1 );

                    // reset counter
                    tCount = 0;

                    // fill matrix with IDs
                    for( Basis * tNode : tNodes )
                    {
                        // test if node belongs to neighbor
                        if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_use_owned_and_shared() )
                            {
                                // increment counter
                                tSendID( p )( tCount++ ) = tNode->get_hmr_id();
                            }
                        }
                        else
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_used() )
                            {
                                // increment counter
                                tSendID( p )( tCount++ ) = tNode->get_hmr_id();
                            }
                        }
                    }
                } // end neighbor exists
            } // end loop over all neighbors

            Cell< Matrix< DDLUMat > > tReceiveID;

            // communicate node IDs to neighbors
            communicate_mats( tProcNeighbors,
                    tSendID,
                    tReceiveID );

            // clear memory
            tSendID.clear();

            Cell< Matrix< DDLUMat > > tSendIndex( tNumberOfProcNeighbors, tEmpty );

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                moris_id tNeighbor = tProcNeighbors( p );

                if ( tNeighbor < tNumberOfProcs && tNeighbor != tMyRank )
                {
                    // cell with basis in aura
                    Cell< Basis * > tNodes;

                    // collect nodes within inverse aura
                    this->collect_basis_from_aura( p, 1, tNodes );

                    // create Map
                    map< luint, moris_id > tMap;

                    for( Basis * tNode : tNodes )
                    {
                        if( tNode->get_owner() == tMyRank )
                        {
                            tMap[ tNode->get_hmr_id() ] = tNode->get_hmr_index();
                        }
                    }

                    // get number of nodes
                    uint tNumberOfNodes = tReceiveID( p ).length();

                    // send indices
                    tSendIndex( p ).set_size( tNumberOfNodes, 1 );

                    // fill index with requested IDs
                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        tSendIndex( p )( k ) = tMap.find( tReceiveID( p )( k ) );
                    }
                }
            }

            Cell< Matrix< DDLUMat > > tReceiveIndex;

            // communicate node IDs to neighbors
            communicate_mats( tProcNeighbors,
                    tSendIndex,
                    tReceiveIndex );

            // clear memory
            tSendIndex.clear();

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                moris_id tNeighbor = tProcNeighbors( p );

                if ( tNeighbor < tNumberOfProcs && tNeighbor != tMyRank )
                {
                    // cell with basis in aura
                    Cell< Basis * > tNodes;

                    // collect nodes within aura
                    this->collect_basis_from_aura( p, 0, tNodes );

                    // initialize counter
                    uint tCount = 0;

                    for( Basis* tNode : tNodes )
                    {
                        if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_use_owned_and_shared() )
                            {
                                tNode->set_domain_index( tReceiveIndex( p )( tCount++ ) );
                            }
                        }
                        else
                        {
                            if( tNode->get_owner() == tNeighbor && tNode->is_used() )
                            {
                                tNode->set_domain_index( tReceiveIndex( p )( tCount++ ) );
                            }
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::communicate_missed_node_indices()
    {
        // get number of ranks
        uint tNumberOfProcs = par_size();

        // only needed if in parallel mode
        if ( tNumberOfProcs > 1 )
        {
            // make sure that all basis are unflagged
            this->unflag_all_basis();

            // get number of proc neighbors
            uint tNumberOfProcNeighbors = mBackgroundMesh->get_number_of_proc_neighbors();

            // get proc neighbors from background mesh
            const Matrix< IdMat > & tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            uint tCounter= 0;

            for ( auto tBasis : mAllBasisOnProc )
            {
                if( tBasis->get_hmr_index() == gNoEntityID  )
                {
                    MORIS_ASSERT( tBasis->get_owner() != gNoProcOwner,
                            "Lagrange_Mesh_Base::communicate_missed_node_indices(), Node has no proc owner assigned");

                    tCounter++;
                }
            }

            Cell< Basis* > tBasisWithoutId( tCounter );
            tCounter = 0;

            for ( auto tBasis : mAllBasisOnProc )
            {
                if( tBasis->get_hmr_index() == gNoEntityID  )
                {
                    tBasisWithoutId( tCounter++ ) = tBasis;
                }
            }

            // create cell of matrices to send
            Matrix< DDLUMat > tEmptyLuint;
            Matrix<  DDUMat > tEmptyUint;

            Cell< Matrix< DDLUMat > > tSendAncestor  ( tNumberOfProcNeighbors, tEmptyLuint );
            Cell< Matrix<  DDUMat > > tSendPedigree  ( tNumberOfProcNeighbors, tEmptyUint  );
            Cell< Matrix<  DDUMat > > tSendBasisIndex( tNumberOfProcNeighbors, tEmptyUint  );

            // get my rank
            moris_id tMyRank = par_rank();

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                // get rank of neighbor
                moris_id tNeihgborRank = tProcNeighbors( p );

                if ( tNeihgborRank != tMyRank && tNeihgborRank != gNoProcNeighbor )
                {
                    // calculate addresses of basis to ask for
                    this->encode_foreign_basis_path( tBasisWithoutId,
                            tNeihgborRank,
                            tSendAncestor( p ),
                            tSendPedigree( p ),
                            tSendBasisIndex( p ) );
                }
            }

            // send basis indices
            Cell< Matrix< DDUMat > > tReceiveBasisIndex;

            // communicate basis indices
            communicate_mats(
                    tProcNeighbors,
                    tSendBasisIndex,
                    tReceiveBasisIndex );

            // free memory
            tSendBasisIndex.clear();

            // ancestors to receive
            Cell< Matrix< DDLUMat > > tReceiveAncestor;

            // communicate ancestor list
            communicate_mats(
                    tProcNeighbors,
                    tSendAncestor,
                    tReceiveAncestor );

            // free memory
            tSendAncestor.clear();

            // communicate pedigree list
            Cell< Matrix<  DDUMat > > tReceivePedigree;

            communicate_mats(
                    tProcNeighbors,
                    tSendPedigree,
                    tReceivePedigree );

            // free memory
            tSendPedigree.clear();

            // matrix with owners to send
            Cell< Matrix<  DDUMat > > tSendId( tNumberOfProcNeighbors, tEmptyUint );

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                // get number of basis requested by neighbor
                luint tNumberOfBasis = tReceiveBasisIndex( p ).length();

                // initialize matrix
                tSendId( p ).set_size( tNumberOfBasis, 1 );

                // initialize memory conter
                luint tMemoryCounter = 0;

                // loop over all basis
                for( luint k = 0; k < tNumberOfBasis; ++k )
                {
                    // pick requested element
                    Background_Element_Base*
                    tElement = mBackgroundMesh->decode_pedigree_path(
                            tReceiveAncestor( p )( k ),
                            tReceivePedigree( p ),
                            tMemoryCounter );

                    // pick requested basis
                    Basis * tBasis = mAllElementsOnProc( tElement->get_memory_index() )->
                            get_basis( tReceiveBasisIndex( p )( k ) );

                    // write basis owner into send array
                    tSendId( p )( k ) = tBasis->get_hmr_index();
                }
            }

            // free memory
            tReceiveBasisIndex.clear();
            tReceiveAncestor  .clear();
            tReceivePedigree  .clear();

            // communicate owners
            Cell< Matrix<  DDUMat > > tReceiveId;

            communicate_mats(
                    tProcNeighbors,
                    tSendId,
                    tReceiveId );

            // free memory
            tSendId.clear();

            // loop over all proc neighbors
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                // get rank of neighbor
                moris_id tNeihgborRank = tProcNeighbors( p );

                if ( tNeihgborRank != tMyRank && tNeihgborRank != gNoProcNeighbor )
                {
                    // initialize counter
                    luint tCount = 0;

                    // count number of basis suspected to be owned by neighbor
                    for( auto tBasis : tBasisWithoutId )
                    {
                        if ( tBasis->get_owner() == tNeihgborRank )
                        {
                            // set ownership from received matrix
                            tBasis->set_domain_index( tReceiveId( p )( tCount++ ) );
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::calculate_node_sharing()
    {
        moris_id tNumberOfProcs = par_size();
        //            moris_id tMyProcRank    = par_rank();

        if( tNumberOfProcs == 1 ) // serial mode
        {
            return; // Do nothing (because node sharing is only in parallel)
        }

        if( tNumberOfProcs > 1 )
        {
            // get proc neighbors from background mesh
            const Matrix< IdMat > & tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            // get number of proc neighbors
            uint tNumberOfProcNeighbors = mBackgroundMesh->get_number_of_proc_neighbors();

            // Initialize neighbor proc rank
            uint tNeighborProcRank = 0;

            // Initialize an inverse aura flag
            Cell< Basis* > tAuraNodes;
            Cell< Basis* > tInverseAuraNodes;
            bool           tInverseAuraFlag = true;

            // Counter
            uint tCount = 0;

            // Id
            moris_id tBasisMemoryIndex = 0;

            // Iterate through the processors that we share nodes with
            for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
            {
                // neighbor proc rank
                tNeighborProcRank = tProcNeighbors(p);

                // cell with basis in aura
                tInverseAuraFlag = false;

                // collect nodes within aura
                this->collect_basis_from_aura( p, tInverseAuraFlag, tAuraNodes );

                // flip the aura flag to get inverse flag
                tInverseAuraFlag = true;

                // Collect nodes within inverse aura
                this->collect_basis_from_aura( p, tInverseAuraFlag, tInverseAuraNodes );

                // Get all the ids for the aura and inverse aura nodes
                // The intersection  of the two sets is the nodes on the boundary
                // Cells instead of moris matrix for intersection operation on std::vector
                Cell< moris_id > tAuraNodeMemoryIndex( tAuraNodes.size() );
                Cell< moris_id > tInverseAuraNodeMemoryIndex( tInverseAuraNodes.size() );
                Cell< moris_id > tBoundaryNodeMemoryIndex( std::max( tAuraNodes.size(), tInverseAuraNodes.size() ) );

                // Reset count
                tCount = 0;

                // Iterate through aura nodes and collect ids
                for( auto tBasis : tAuraNodes )
                {
                    tBasisMemoryIndex = tBasis->get_memory_index();

                    //FIXME: Nodes that are in the aura do not all have ids.
                    if(tBasisMemoryIndex != 0)
                    {
                        // Add node id to list
                        tAuraNodeMemoryIndex( tCount++ ) = tBasisMemoryIndex;
                    }
                }

                // Size out extra space
                tAuraNodeMemoryIndex.resize( tCount );

                // Reset count
                tCount = 0;
                // Iterate through inverse aura nodes and collect ids
                for( auto tBasis : tInverseAuraNodes )
                {
                    tBasisMemoryIndex = tBasis->get_memory_index();

                    if( tBasisMemoryIndex != 0 )
                    {
                        // Add node id to list
                        tInverseAuraNodeMemoryIndex( tCount++ ) = tBasisMemoryIndex;
                    }
                }
                // Size out extra space
                tInverseAuraNodeMemoryIndex.resize( tCount );

                // sort inverse aura and aura node ids
                std::sort( tAuraNodeMemoryIndex.data().begin(), tAuraNodeMemoryIndex.data().end() );
                std::sort( tInverseAuraNodeMemoryIndex.data().begin(), tInverseAuraNodeMemoryIndex.data().end() );

                std::vector< moris_id >::iterator it =
                        std::set_intersection( tAuraNodeMemoryIndex.data().begin(),
                                tAuraNodeMemoryIndex.data().end(),
                                tInverseAuraNodeMemoryIndex.data().begin(),
                                tInverseAuraNodeMemoryIndex.data().end(),
                                tBoundaryNodeMemoryIndex.data().begin() );

                // resize boundary node data
                tBoundaryNodeMemoryIndex.data().resize( it - tBoundaryNodeMemoryIndex.data().begin() );

                // Add basis sharing information to the basis itself;
                for( auto tBoundNodeMemoryIndex : tBoundaryNodeMemoryIndex )
                {
                    Basis* tBoundBasis = get_basis_by_memory_index( tBoundNodeMemoryIndex );
                    tBoundBasis->add_node_sharing( tNeighborProcRank );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    Element * Lagrange_Mesh_Base::get_child(       Element * aElement,
            uint aChildIndex )
    {
        // get pointer to background element
        Background_Element_Base* tBackElement = aElement->get_background_element();

        if ( tBackElement->has_children() )
        {
            // get child of background element
            Background_Element_Base* tBackChild = tBackElement->get_child( aChildIndex );

            // grab child from element list
            return mAllElementsOnProc( tBackChild->get_memory_index() );
        }
        else
        {
            // return nothing
            return nullptr;
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::reset_fields()
    {
        mRealScalarFieldLabels       .clear();
        mRealScalarFieldData         .clear();
        mRealScalarFieldBSplineCoeffs.clear();
        mRealScalarFieldRanks        .clear();
        mRealScalarFieldBSplineOrders.clear();

        Matrix< DDRMat > tEmpty;

        // first field is element level
        mRealScalarFieldLabels       .push_back( "Element_Level" );
        mRealScalarFieldRanks        .push_back( EntityRank::ELEMENT );
        mRealScalarFieldData         .push_back( tEmpty );
        mRealScalarFieldBSplineCoeffs.push_back( tEmpty );
        mRealScalarFieldBSplineOrders.push_back( 0 );

        // second field is element owner
        mRealScalarFieldLabels       .push_back( "Element_Owner" );
        mRealScalarFieldRanks        .push_back( EntityRank::ELEMENT );
        mRealScalarFieldData         .push_back( tEmpty );
        mRealScalarFieldBSplineCoeffs.push_back( tEmpty );
        mRealScalarFieldBSplineOrders.push_back( 0 );

        // third field is vertex IDs
        mRealScalarFieldLabels       .push_back( "Node_IDs" );
        mRealScalarFieldRanks        .push_back( EntityRank::NODE );
        mRealScalarFieldData         .push_back( tEmpty );
        mRealScalarFieldBSplineCoeffs.push_back( tEmpty );
        mRealScalarFieldBSplineOrders.push_back( 0 );

        // forth field is vertex IDs
        //            mRealScalarFieldLabels       .push_back( "Element_Indices" );
        //            mRealScalarFieldRanks        .push_back( EntityRank::ELEMENT );
        //            mRealScalarFieldData         .push_back( tEmpty );
        //            mRealScalarFieldBSplineCoeffs.push_back( tEmpty );
        //            mRealScalarFieldBSplineOrders.push_back( 0 );

    }

    //------------------------------------------------------------------------------

    STK * Lagrange_Mesh_Base::create_stk_object( const double aTimeStep )
    {
        MORIS_ERROR( mOrder <= 2 , "Tried to create an STK object for third or higher order. \n This is not supported by Exodus II.");

        // create new MTK object
        STK* aSTK = new STK( this );

        // create data
        aSTK->create_mesh_data( aTimeStep );

        // return MTK object
        return aSTK;
    }

    //------------------------------------------------------------------------------

    bool Lagrange_Mesh_Base::test_for_double_nodes()
    {
        // strategy: fill a matrix with node IDs. Make them unique.
        // each node must appear only once

        this->determine_elements_connected_to_basis();

        // get numnber of nodes
        luint tNumberOfNodes = mAllBasisOnProc.size();

        // matrix which will contain node IDs
        Matrix< DDLUMat > tNodeIDs( tNumberOfNodes, 1 );

        // loop over all nodes
        for( luint k = 0; k < tNumberOfNodes; ++k )
        {
            // get node
            Basis * tNode = mAllBasisOnProc( k );

            // get level of node
            luint tLevel = tNode->get_level();

            if ( tLevel == 0 )
            {
                tNodeIDs( k ) = tNode->get_hmr_id();
            }
            else
            {
                // get ijk of node
                const luint* tNodeIJK = tNode->get_ijk();

                // copy array into writable array
                luint tIJK[ 3 ];

                for( uint i = 0; i < mNumberOfDimensions; ++i )
                {
                    tIJK[ i ] = tNodeIJK[ i ];
                }

                bool tCheck = true;

                // now see if there is any node above
                while ( tLevel > 0 && tCheck )
                {
                    for( uint i = 0; i<mNumberOfDimensions; ++i )
                    {
                        tCheck = tCheck && ( tIJK[ i ] % 2 == 0 );
                    }
                    if ( tCheck )
                    {
                        // go up
                        for( uint i = 0; i < mNumberOfDimensions; ++i )
                        {
                            tIJK[ i ] /= 2;
                        }

                        // decrement level
                        --tLevel;
                    }
                }

                // calculate new id of node
                if ( mNumberOfDimensions == 1 )
                {
                    tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ] );
                }
                else if ( mNumberOfDimensions == 2 )
                {
                    tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ], tIJK[ 1 ] );
                }
                else if ( mNumberOfDimensions == 3 )
                {
                    tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ], tIJK[ 1 ], tIJK[ 2 ] );
                }
            }
        }

        /*        // init a boost bitset
                BoostBitset tBitset( tNodeIDs.max()+1 );

                // loop over all nodes
                for( luint k=0; k<tNumberOfNodes; ++k )
                {
                    // write id of node into array

                    // test if ID was already flagged
                    if ( ! tBitset.test(  tNodeIDs( k ) ) )
                    {
                        tBitset.set( tNodeIDs( k ) );
                    }
                    else
                    {
                        // other coordinate
                        luint j=0;
                        // search for corresponding node
                        for( luint i=0; i<tNumberOfNodes; ++i )
                        {
                            if ( tNodeIDs( i ) == tNodeIDs( k ) )
                            {
                                j = i;
                                break;
                            }

                        }

                        std::fprintf( stdout,
                                "Error: Node %lu has the ID %lu, which is already used bu node %lu.\n\n",
                                ( long unsigned int ) k,
                                ( long unsigned int ) tNodeIDs( k ),
                                ( long unsigned int ) j
                        );

                        // print elements of node 1
                        std::fprintf( stdout, "Elements connected to Node %lu : \n", ( long unsigned int ) k );

                        Matrix< DDLUMat > tElements = mElementsPerNode( mAllBasisOnProc( k )->get_hmr_index() );
                         for( luint i=0; i<tElements.length(); ++i )
                        {
                            // get element
                            Element* tElement = mAllElementsOnProc ( tElements( i ) );

                            std::fprintf( stdout, "    Element %4lu    ID %4lu      Parent:  %4lu\n",
                                    ( long unsigned int ) tElement->get_background_element()->get_subdomain_index(),
                                    ( long unsigned int ) tElement->get_background_element()->get_hmr_id(),
                                    ( long unsigned int ) tElement->get_background_element()->get_parent()->get_hmr_id() );
                        }
                        std::fprintf( stdout, "\n" );

                        // print elements of node2
                        std::fprintf( stdout, "Elements connected to Node %lu : \n", ( long unsigned int ) j );

                        tElements = mElementsPerNode( mAllBasisOnProc( j )->get_hmr_index() );

                        for( luint i=0; i<tElements.length(); ++i )
                        {
                            // get element
                            Element* tElement = mAllElementsOnProc ( tElements( i ) );

                            std::fprintf( stdout, "    Element %4lu    ID %4lu     Parent:  %4lu\n",
                                    ( long unsigned int ) tElement->get_background_element()->get_subdomain_index(),
                                    ( long unsigned int ) tElement->get_background_element()->get_hmr_id(),
                                    ( long unsigned int ) tElement->get_background_element()->get_parent()->get_hmr_id()  );
                        }
                        std::fprintf( stdout, "\n" );

                        exit(-1);
                    }
                } */

        // make matrix unique
        Matrix< DDLUMat > tNodeUniqueIDs;
        unique( tNodeIDs, tNodeUniqueIDs );

        // make sure that number of nodes is the same
        return tNodeUniqueIDs.length() == tNumberOfNodes;
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_to_gmsh( const std::string & aFilePath )
    {
        // start timer
        tic tTimer;

        // get my rank
        moris_id tMyRank = par_rank();

        // modify filename
        std::string tFilePath;
        if ( moris::par_size() > 1 )
        {
            tFilePath = aFilePath + "." +  std::to_string( par_size() ) + "." +  std::to_string( tMyRank );
        }
        else
        {
            tFilePath = aFilePath;
        }

        // create output file
        std::FILE * tFile = std::fopen( tFilePath.c_str(), "w+");

        // write header
        std::fprintf( tFile, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

        // write node coordinates
        std::fprintf( tFile, "$Nodes\n%lu\n", ( long unsigned int ) mNumberOfUsedNodes );

        // get mesh scale factor
        real tScale = mParameters->get_gmsh_scale();

        switch ( mParameters->get_number_of_dimensions() )
        {
            case( 2 ) :
                        {
                // loop over all nodes
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if this node is relevant
                    if ( tNode->is_used() )
                    {
                        // get coordinates of node
                        const real* tXY = tNode->get_xyz();

                        // write coordinates to ASCII file
                        std::fprintf( tFile,
                                "%lu %.17f %.17f 0\n",
                                ( long unsigned int ) tNode->get_hmr_index()+1,
                                ( double ) tXY[ 0 ]*tScale,
                                ( double ) tXY[ 1 ]*tScale );
                    }
                }
                break;
                        }
            case( 3 ) :
                        {
                // loop over all nodes
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if this node is relevant
                    if ( tNode->is_used() )
                    {
                        // get coordinates of node
                        const real* tXYZ = tNode->get_xyz();

                        // write coordinates to ASCII file
                        std::fprintf( tFile,
                                "%lu %.17f %.17f %.17f\n",
                                ( long unsigned int ) tNode->get_hmr_index()+1,
                                ( double ) tXYZ[ 0 ]*tScale,
                                ( double ) tXYZ[ 1 ]*tScale,
                                ( double ) tXYZ[ 2 ]*tScale );
                    }
                }
                break;
                        }
            default :
            {
                MORIS_ERROR( false, "wrong number of dimensions\n");
                break;
            }
        }

        // end node tag
        std::fprintf( tFile, "$EndNodes\n" );

        // write element tag
        std::fprintf( tFile, "$Elements\n%lu\n",
                ( long unsigned int ) mNumberOfElements );

        // loop over all elements
        for( auto tElement: mAllElementsOnProc )
        {
            // test if this element is relevant
            if ( tElement->get_owner() == tMyRank && tElement->is_active() )
            {
                // print element line
                std::fprintf( tFile, "%lu %s\n",
                        ( long unsigned int ) tElement->get_hmr_index()+1,
                        tElement->get_gmsh_string().c_str() );
            }
        }

        // finish element tag
        std::fprintf( tFile, "$EndElements\n" );

        // close file
        std::fclose( tFile );

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // print output
        MORIS_LOG_INFO( "%s Created GMSH File: %s",
                proc_string().c_str(),
                tFilePath.c_str() );
        MORIS_LOG_INFO( "Writing took %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );

    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::link_twins( )
    {
        // get number of elements of interest
        auto tNumberOfElements = this->get_number_of_elements();

        // get number of meshes
        uint tNumberOfTwins = mBSplineMeshes.size();

        // allocate twin container
        for( uint e=0; e<tNumberOfElements; ++e )
        {
            this->get_element( e )->allocate_twin_container( tNumberOfTwins );
        }

        // loop over list of entries of lagrange_to_bspline map set by user
        for( uint k = 1; k < tNumberOfTwins; ++k )
        {
            // check if a B-Spline mesh has been assigned to current Lagrange mesh ...
            if( mBSplineMeshes( k ) != nullptr )
            {
                // ... if so, match twins for all elements in B-Sp. and Lag. meshes
                for( uint e=0; e<tNumberOfElements; ++e )
                {
                    this->get_element( e )->set_twin( k, mBSplineMeshes( k )->get_element( e ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_to_vtk( const std::string & aFilePath )
    {
        // start timer
        tic tTimer;

        // modify filename
        std::string tFilePath = parallelize_path( aFilePath );

        // open the file
        std::ofstream tFile(tFilePath, std::ios::binary);

        // containers
        //float tFValue = 0;
        //int   tIValue = 0;
        float tFChar = 0;
        int   tIChar = 0;

        tFile << "# vtk DataFile Version 3.0" << std::endl;
        tFile << "GO BUFFS!" << std::endl;
        tFile << "BINARY" << std::endl;
        //tFile << "ASCII" << std::endl;
        luint tNumberOfNodes = mAllBasisOnProc.size();

        // write node data
        tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;

        // ask settings for numner of dimensions
        auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

        if ( tNumberOfDimensions == 2 )
        {
            // loop over all nodes
            for ( luint k = 0; k < tNumberOfNodes; ++k )
            {
                // get coordinate from node
                const real* tXY = mAllBasisOnProc( k )->get_xyz();

                // write coordinates to mesh
                tFChar = swap_byte_endian( (float) tXY[ 0 ] );
                tFile.write( (char*) &tFChar, sizeof(float));
                tFChar = swap_byte_endian( (float) tXY[ 1 ] );
                tFile.write( (char*) &tFChar, sizeof(float));
                tFChar = swap_byte_endian( (float) 0 );
                tFile.write( (char*) &tFChar, sizeof(float));
            }
        }
        else if ( tNumberOfDimensions == 3 )
        {
            // loop over all nodes
            for ( luint k = 0; k < tNumberOfNodes; ++k )
            {
                // get coordinate from node
                const real* tXYZ = mAllBasisOnProc( k )->get_xyz();

                // write coordinates to mesh
                tFChar = swap_byte_endian( (float) tXYZ[ 0 ] );
                tFile.write( (char*) &tFChar, sizeof(float));
                tFChar = swap_byte_endian( (float) tXYZ[ 1 ] );
                tFile.write( (char*) &tFChar, sizeof(float));
                tFChar = swap_byte_endian( (float) tXYZ[ 2 ] );
                tFile.write( (char*) &tFChar, sizeof(float));
            }
        }

        tFile << std::endl;

        // write element topology
        int tCellType = mAllElementsOnProc( 0 )->get_vtk_type();

        // count number of non padding elements
        luint tNumberOfElements = 0;

        // can only write element data if vtk map exists
        if ( tCellType != 0 )
        {
            int tNumberOfNodesPerElement = swap_byte_endian( (int) mNumberOfBasisPerElement );

            luint tNumberOfAllElementsOnProc = mAllElementsOnProc.size();

            for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    // increment element counter
                    ++tNumberOfElements;
                }
            }

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                    << ( mNumberOfBasisPerElement + 1 )*tNumberOfElements  << std::endl;

            // matrix containing node indices
            Matrix< DDLUMat > tNodes( mNumberOfBasisPerElement, 1 );

            // loop over all elements
            for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    tFile.write( (char*) &tNumberOfNodesPerElement, sizeof(int)) ;

                    // ask element for nodes
                    mAllElementsOnProc( k )->get_basis_indices_for_vtk( tNodes );

                    for( uint i=0; i <mNumberOfBasisPerElement; ++i )
                    {
                        tIChar = swap_byte_endian( (int) tNodes( i ) );
                        tFile.write((char *) &tIChar, sizeof(int));
                        //tFile << " " << (int) mAllElementsOnProc( k )->get_basis( i )->get_memory_index();
                    }
                    //tFile << std::endl;
                }
            }

            tFile << std::endl;

            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
            tIChar = swap_byte_endian( tCellType );
            for ( luint k = 0; k < tNumberOfElements; ++k)
            {
                tFile.write( (char*) &tIChar, sizeof(int));
                //tFile << tCellType << std::endl;
            }

            // write element data
            tFile << "CELL_DATA " << tNumberOfElements << std::endl;

            // write element ID
            tFile << "SCALARS ELEMENT_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    //tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_hmr_id() );
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_id() );
                    tFile.write( (char*) &tIChar, sizeof(int));
                }
            }
            tFile << std::endl;

            // write element ID
            tFile << "SCALARS ELEMENT_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    //tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_hmr_id() );
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_index() );
                    tFile.write( (char*) &tIChar, sizeof(int));
                }
            }
            tFile << std::endl;

            // write proc owner
            tFile << "SCALARS ELEMENT_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_owner() );
                    tFile.write( (char*) &tIChar, sizeof(float));
                }
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS ELEMENT_LEVEL int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_level() );
                    tFile.write( (char*) &tIChar, sizeof(float));
                }
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS ELEMENT_CHILD_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_child_index() );
                    tFile.write( (char*) &tIChar, sizeof(float));
                }
            }
            tFile << std::endl;

            // write memory index
            tFile << "SCALARS ELEMENT_MEMORY_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
            {
                if ( mAllElementsOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_memory_index() );
                    tFile.write( (char*) &tIChar, sizeof(float));
                }
            }
            tFile << std::endl;
        }

        // write node data
        tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

        tFile << "SCALARS NODE_ID int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
        {

            tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_id() );
            tFile.write( (char*) &tIChar, sizeof(float));
        }
        tFile << std::endl;

        tFile << "SCALARS NODE_INDEX int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
        {

            tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_index() );
            tFile.write( (char*) &tIChar, sizeof(float));
        }
        tFile << std::endl;

        tFile << "SCALARS NODE_OWNER int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
        {

            tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_owner() );
            tFile.write( (char*) &tIChar, sizeof(float));
        }
        tFile << std::endl;

        tFile << "SCALARS DOMAIN_ID int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( moris::uint k = 0; k <  tNumberOfNodes; ++k)
        {

            tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_hmr_id() );
            tFile.write( (char*) &tIChar, sizeof(float));
        }
        tFile << std::endl;

        tFile << "SCALARS DOMAIN_INDEX int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( moris::uint k = 0; k <  tNumberOfNodes; ++k)
        {

            tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_hmr_index() );
            tFile.write( (char*) &tIChar, sizeof(float));
        }
        tFile << std::endl;

        // close the output file
        tFile.close();

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // print output
        MORIS_LOG_INFO( "%s Created VTK debug file.",
                proc_string().c_str());

        MORIS_LOG_INFO( "Mesh has %lu active and refined Elements and %lu Nodes.",
                ( long unsigned int ) tNumberOfElements,
                ( long unsigned int ) tNumberOfNodes );

        MORIS_LOG_INFO("Creation took %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_nodes_on_level_zero()
    {
        // ask background mesh for number of dimensions
        Matrix< DDLUMat > tNumberOfElements = mBackgroundMesh->get_number_of_elements_per_direction_on_proc();

        if( mNumberOfDimensions == 2 )
        {
            // get maximum numbers in i direction
            luint tImax = tNumberOfElements( 0, 0 );

            // get maximum numbers in j direction
            luint tJmax = tNumberOfElements( 1, 0 );

            // initialize element counter
            luint tCount = 0;

            // loop over all elements
            for( luint j=0; j<tJmax; ++j )
            {
                for ( luint i=0; i<tImax; ++i )
                {
                    mNumberOfAllBasis += mAllCoarsestElementsOnProc( tCount++ )->create_basis_on_level_zero( mAllElementsOnProc );
                }
            }
        }
        else if( mNumberOfDimensions == 3 )
        {
            // get maximum numbers in i direction
            luint tImax = tNumberOfElements( 0, 0 );

            // get maximum numbers in j direction
            luint tJmax = tNumberOfElements( 1, 0 );

            // get maximum numbers in k direction
            luint tKmax = tNumberOfElements( 2, 0 );

            // initialize element counter
            luint tCount = 0;

            // loop over all elements
            for( luint k=0; k<tKmax; ++k )
            {
                for( luint j=0; j<tJmax; ++j )
                {
                    for ( luint i=0; i<tImax; ++i )
                    {
                        mNumberOfAllBasis += mAllCoarsestElementsOnProc( tCount++ )->create_basis_on_level_zero( mAllElementsOnProc );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::update_node_list()
    {
        // tidy up memory
        mNodes.clear();

        // assign memory
        mNodes.resize( mNumberOfUsedNodes, nullptr );

        // initialize counter
        luint tCount = 0;

        // reset max level
        mMaxLevel = 0;

        if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
        {
            for( auto tNode : mAllBasisOnProc )
            {
                if ( tNode->is_use_owned_and_shared() )
                {
                    mNodes( tCount++ ) = tNode;

                    mMaxLevel = std::max( tNode->get_level(), mMaxLevel );
                }
            }
        }
        else
        {
            for( auto tNode : mAllBasisOnProc )
            {
                if ( tNode->is_used() )
                {
                    mNodes( tCount++ ) = tNode;

                    mMaxLevel = std::max( tNode->get_level(), mMaxLevel );
                }
            }
        }

        MORIS_ERROR( tCount == mNumberOfUsedNodes, "Number Of used Nodes does not match" );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_facets()
    {
        tic tTimer;

        // get my rank
        moris_id tMyRank = par_rank();

        // delete existing lagrange facets
        this->delete_facets();

        // step 1: unflag all facets

        // loop over all elements
        for( Element * tElement :  mAllElementsOnProc )
        {
            // make sure that all background faces are unflagged
            tElement->get_background_element()->reset_flags_of_facets();
        }

        // step 2 : determine number of facets per element
        uint tNumberOfFacetsPerElement = 0;
        if ( mParameters->get_number_of_dimensions() == 2 )
        {
            tNumberOfFacetsPerElement = 4;
        }
        else if( mParameters->get_number_of_dimensions() == 3 )
        {
            tNumberOfFacetsPerElement = 6;
        }

        // step 2: count number of active or refined facets on mesh

        // initialize counter
        uint tCount = 0;

        // loop over all active elements
        for( Element * tElement : mAllElementsOnProc )
        {
            // test if element is not padding
            if( ! tElement->is_deactive() && ! tElement->is_padding() )
            {
                // get pointer to Element
                Background_Element_Base * tBackElement = tElement->get_background_element();

                for( uint f=0; f<tNumberOfFacetsPerElement; ++f )
                {
                    // get pointer to face
                    Background_Facet * tBackFacet = tBackElement->get_facet( f );

                    MORIS_ASSERT(tBackFacet != nullptr, " background facet is nullptr");

                    // test if background facet is not flagged and element
                    if( ! tBackFacet->is_flagged() )
                    {
                        // flag facet
                        tBackFacet->flag();

                        // increment counter
                        ++tCount;
                    }
                }
            }
        }

        // step 2: create lagrange facets
        mFacets.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // counter for owned facets
        uint tOwnedCount = 0;

        // loop over all active elements
        for( Element * tElement : mAllElementsOnProc )
        {
            // pick pointer to element
            if( ! tElement->is_deactive() && ! tElement->is_padding() )
            {
                Background_Element_Base * tBackElement = tElement->get_background_element();

                for( uint f=0; f<tNumberOfFacetsPerElement; ++f )
                {
                    // get pointer to facet
                    Background_Facet * tBackFacet = tBackElement->get_facet( f );

                    // test if facet is flagged
                    if( tBackFacet->is_flagged() )
                    {
                        // create facet
                        Facet * tFacet = this->create_facet( tBackFacet );

                        MORIS_ASSERT(tFacet != nullptr, " Facet facet is nullptr");

                        // test owner of facet
                        if( tFacet->get_owner() == tMyRank )
                        {
                            tFacet->set_id( tOwnedCount++ );
                        }

                        // set index for this facet
                        tFacet->set_index( tCount );

                        // copy facet into array
                        mFacets( tCount++ ) = tFacet;

                        // unflag facet
                        tBackFacet->unflag();
                    }
                }
            }
        }

        // step 5: write facets into cells
        for( Facet * tFacet : mFacets )
        {
            // get leader
            Element * tLeader = tFacet->get_hmr_leader();

            // get follower
            Element * tFollower = tFacet->get_hmr_follower();

            // leader is always active
            tLeader->set_hmr_facet( tFacet,
                    tFacet->get_index_on_leader() );

            if( tFollower != nullptr )
            {
                // insert element into follower
                tFollower->set_hmr_facet( tFacet,
                        tFacet->get_index_on_follower() );
            }
        }

        // step 6: synchronize proc IDs if parallel
        if( par_size() > 1 )
        {
            this->synchronize_facet_ids( tOwnedCount );
        }

        // step 7 : link facets to basis

        // reset facet containers
        for( Basis * tBasis : mAllBasisOnProc )
        {
            tBasis->delete_facet_container();
        }

        // count facets and increment each ID by 1, because IDs are supposed to
        // be 1-based
        for( Facet * tFacet : mFacets )
        {
            // only connect active facets
            if ( tFacet->is_active() )
            {
                // get number of connected basis
                uint tNumberOfBasis = tFacet->get_number_of_vertices();

                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tFacet->get_basis( k )->increment_facet_counter();
                }
            }

            // increment facet ID
            tFacet->increment_id();
        }

        // insert facet containers
        for( Basis * tBasis : mAllBasisOnProc )
        {
            tBasis->init_facet_container();
        }

        for( Facet * tFacet : mFacets )
        {
            // only connect active facets
            if ( tFacet->is_active() )
            {
                // get number of connected basis
                uint tNumberOfBasis = tFacet->get_number_of_vertices();

                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tFacet->get_basis( k )->insert_facet( tFacet );
                }
            }
        }

        /*std::cout << par_rank() << " flag 1" << std::endl;
        // step 7 : link facets with children
        if( mParameters->get_number_of_dimensions() == 2 )
        {
            this->link_facet_children_2d();
        }
        else if ( mParameters->get_number_of_dimensions() == 3 )
        {
            this->link_facet_children_3d();
        }
        std::cout << par_rank() << " flag 2" << std::endl; */

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        MORIS_LOG_INFO( "%s Created Faces for Lagrange Mesh.",
                proc_string().c_str());

        MORIS_LOG_INFO( "Creation %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }
    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_facet_clusters()
    {
    }
    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::create_edges()
    {
        tic tTimer;

        // get my rank
        moris_id tMyRank = par_rank();

        // delete existing edges
        this->delete_edges();

        // step 1: unflag all edges on background mesh
        for( Element * tElement :  mAllElementsOnProc )
        {
            // make sure that all background faces are unflagged
            tElement->get_background_element()->reset_flags_of_edges();
        }

        // step 2: count number of active or refined edges on mesh

        // initialize counter
        uint tCount = 0;

        // loop over all active elements
        for( Element * tElement : mAllElementsOnProc )
        {
            // test if element is not deactive
            if( ! tElement->is_deactive() && ! tElement->is_padding() )
            {
                // get pointer to Element
                Background_Element_Base * tBackElement = tElement->get_background_element();

                // loop over all edges
                for( uint e=0; e<12; ++e )
                {
                    // get pointer to edge
                    Background_Edge * tBackEdge = tBackElement->get_edge( e );

                    // test if edge is not flagged
                    if( ! tBackEdge->is_flagged() )
                    {
                        // flag edge
                        tBackEdge->flag();

                        // increment counter
                        ++tCount;
                    }
                }
            }
        }

        // step 3: create Lagrange edges

        mEdges.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all active elements
        for( Element * tElement : mAllElementsOnProc )
        {
            // pick pointer to element
            if( ! tElement->is_deactive() && ! tElement->is_padding() )
            {
                Background_Element_Base * tBackElement = tElement->get_background_element();

                for( uint e=0; e<12; ++e )
                {
                    // get pointer to facet
                    Background_Edge* tBackEdge = tBackElement->get_edge( e );

                    // test if facet is flagged
                    if( tBackEdge->is_flagged() )
                    {

                        // create edge
                        Edge * tEdge = this->create_edge( tBackEdge );

                        // set index for this facet
                        tEdge->set_index( tCount );

                        // copy facet into array
                        mEdges( tCount++ ) = tEdge;

                        // unflag edge
                        tBackEdge->unflag();
                    }
                }
            }
        }

        // step 5: write edges into cells
        for( Edge * tEdge: mEdges )
        {
            // get number of elements
            uint tNumberOfElements = tEdge->get_number_of_elements();

            // loop over all elements of this edge
            for( uint e = 0; e<tNumberOfElements; ++e )
            {
                // insert edge into element
                tEdge->get_element( e )->set_hmr_edge( tEdge,
                        tEdge->get_index_on_element( e ) );
            }
        }

        // fix ownership of uncertain edges
        // count owned edges

        uint tOwnedCount = 0;

        if( par_size() > 1 )
        {
            this->negotiate_edge_ownership();

            for( Edge * tEdge: mEdges )
            {
                if( tEdge->get_owner() == tMyRank )
                {
                    tEdge->set_id( tOwnedCount++ );
                }
            }
            this->synchronize_edge_ids( tOwnedCount );
        }
        else
        {
            for( Edge * tEdge: mEdges )
                if( tEdge->get_owner() == tMyRank )
                {
                    tEdge->set_id( tOwnedCount++ );
                }
        }

        // step 7 : link edges to basis
        // reset edge containers
        for( Basis * tBasis : mAllBasisOnProc )
        {
            tBasis->delete_edge_container();
        }

        // count edges
        for( Edge * tEdge : mEdges )
        {
            // only connect active edges
            if ( tEdge->is_active() )
            {
                // get number of connected basis
                uint tNumberOfBasis = tEdge->get_number_of_vertices();

                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tEdge->get_basis( k )->increment_edge_counter();
                }
            }

            // increment edge ID because edge IDs are 1-based
            tEdge->set_id( tEdge->get_id() + 1 );
        }

        // insert edge containers
        for( Basis * tBasis : mAllBasisOnProc )
        {
            tBasis->init_edge_container();
        }

        for( Edge * tEdge : mEdges )
        {
            // only connect active edges
            if ( tEdge->is_active() )
            {
                // get number of connected basis
                uint tNumberOfBasis = tEdge->get_number_of_vertices();

                for( uint k=0; k<tNumberOfBasis; ++k )
                {
                    tEdge->get_basis( k )->insert_edge( tEdge );
                }
            }
        }

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        MORIS_LOG_INFO( "%s Created Edges for Lagrange Mesh.",
                proc_string().c_str());

        MORIS_LOG_INFO( "Creation %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );

        MORIS_LOG_INFO( " " );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::negotiate_edge_ownership()
    {
        moris_id tParSize = par_size();

        if( tParSize > 1 )
        {
            // step 1: create map for neighbors
            Matrix< IdMat > tProcMap( tParSize, 1, tParSize );

            // get proc neighbors from background mesh
            auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

            // number of proc neighbors
            uint tNumberOfNeighbors = tProcNeighbors.length();

            // loop over all neighbors
            for( uint p=0; p<tNumberOfNeighbors; ++p )
            {
                if( tProcNeighbors( p ) < tParSize )
                {
                    tProcMap( tProcNeighbors( p ) ) = p;
                }
            }

            // step 2: determine memory for matrices to send

            Matrix< DDUMat > tElementCount( tNumberOfNeighbors, 1, 0 );
            Matrix< DDLUMat > tMemoryCount( tNumberOfNeighbors, 1, 0 );

            for( Edge * tEdge : mEdges )
            {
                // number of proc
                uint p = tProcMap( tEdge->get_owner() );

                // increment element counter
                ++tElementCount( p );

                // incrememet memory counter
                tMemoryCount( p ) += tEdge->get_hmr_leader()
                                                ->get_background_element()->get_length_of_pedigree_path();
            }

            // create cell of matrices to send
            Matrix< DDLUMat > tEmptyLuint;
            Cell< Matrix< DDLUMat > > tAncestorListSend;
            tAncestorListSend.resize( tNumberOfNeighbors, { tEmptyLuint } );

            Matrix< DDUMat > tEmptyUint;
            Cell< Matrix< DDUMat > > tPedigreeListSend;
            tPedigreeListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

            Cell< Matrix< DDUMat > > tEdgeIndexListSend;
            tEdgeIndexListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

            // step 3: create matrices to send

            // allocate matrices
            for( uint p=0; p< tNumberOfNeighbors; ++p )
            {
                tEdgeIndexListSend( p ).set_size( tElementCount( p ), 1 );
                tAncestorListSend( p ).set_size( tElementCount( p ), 1 );
                tPedigreeListSend( p ).set_size( tMemoryCount( p ), 1 );
            }

            // reset counters
            tElementCount.fill( 0 );
            tMemoryCount.fill( 0 );

            for( Edge * tEdge : mEdges )
            {
                // index of proc
                uint p = tProcMap( tEdge->get_owner() );

                // save index on leader
                tEdgeIndexListSend( p )( tElementCount( p ) )
                                            = tEdge->get_index_on_leader();

                // calculate path o
                tEdge->get_hmr_leader()->get_background_element()
                                                               ->endcode_pedigree_path( tAncestorListSend( p )( tElementCount( p )++ ),
                                                                       tPedigreeListSend( p ),
                                                                       tMemoryCount( p ) );
            }

            // step 4: communicate matrices

            // communicate edge Indices to neighbors

            Cell< Matrix< DDUMat > > tEdgeIndexListReceive;
            communicate_mats( tProcNeighbors,
                    tEdgeIndexListSend,
                    tEdgeIndexListReceive );

            // clear memory
            tEdgeIndexListSend.clear();

            // communicate ancestors to neighbors
            Cell< Matrix< DDLUMat > > tAncestorListReceive;
            communicate_mats( tProcNeighbors,
                    tAncestorListSend,
                    tAncestorListReceive );

            // clear memory
            tAncestorListSend.clear();

            // communicate path to neighbors
            Cell< Matrix< DDUMat > > tPedigreeListReceive;
            communicate_mats( tProcNeighbors,
                    tPedigreeListSend,
                    tPedigreeListReceive );

            // clear memory
            tPedigreeListSend.clear();

            Cell< Matrix< DDUMat > > tOwnerListSend;
            tOwnerListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

            // loop over all received lists
            for ( uint p=0; p<tNumberOfNeighbors; ++p )
            {
                // get number of elements on refinement list
                luint tNumberOfElements = tAncestorListReceive( p ).length();

                // reset memory counter
                luint tMemoryCounter = 0;

                // resize  sending list
                tOwnerListSend( p ).set_size( tNumberOfElements, 1 );

                // loop over all received elements
                for ( uint k=0; k<tNumberOfElements; ++k )
                {
                    // decode path and get pointer to element
                    Background_Element_Base*
                    tBackElement = mBackgroundMesh->decode_pedigree_path( tAncestorListReceive( p )( k ),
                            tPedigreeListReceive( p ),
                            tMemoryCounter );

                    // get pointer to leader
                    Element * tLeader = this->get_element_by_memory_index( tBackElement->get_memory_index() );

                    // get pointer to facet
                    Edge * tEdge = tLeader->get_hmr_edge( tEdgeIndexListReceive( p )( k ) );

                    // copy owner into matrix to send
                    tOwnerListSend( p )( k ) = tEdge->get_owner();
                }
            }  /* end loop over all procs */

            // reset receive lists
            tAncestorListReceive .clear();
            tPedigreeListReceive .clear();
            tEdgeIndexListReceive.clear();

            Cell< Matrix< DDUMat > > tOwnerListReceive;

            // communicate mats
            // note: this is a lot of communication. A more elegant way
            //       could be to just end the edge owners that have changed
            communicate_mats( tProcNeighbors,
                    tOwnerListSend,
                    tOwnerListReceive );

            // clear memory
            tOwnerListSend.clear();

            // reset element counter
            tElementCount.fill( 0 );

            moris_id tMyRank = par_rank();

            // loop over all edges
            for( Edge * tEdge : mEdges )
            {
                if( tEdge->get_owner() != tMyRank )
                {
                    uint p = tProcMap( tEdge->get_owner() );

                    // fix ownership of this edge
                    tEdge->set_owner( tOwnerListReceive( p )( tElementCount( p )++) );
                }
            }
        } // end if parallel
    }

    //------------------------------------------------------------------------------

    Facet * Lagrange_Mesh_Base::create_facet( Background_Facet * aFacet )
    {
        MORIS_ERROR( false, "create_facet() must not be called from base class" );
        return nullptr;
    }

    //------------------------------------------------------------------------------

    Edge * Lagrange_Mesh_Base::create_edge( Background_Edge * aEdge )
    {
        MORIS_ERROR( false, "create_edge() must not be called from base class" );
        return nullptr;
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::delete_facets()
    {
        for( auto tFacet : mFacets )
        {
            delete tFacet;
        }

        mFacets.clear();
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::delete_edges()
    {
        for( auto tEdge : mEdges )
        {
            delete tEdge;
        }

        mEdges.clear();
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::synchronize_facet_ids( uint aOwnedCount )
    {

        // get number of procs
        moris_id tNumberOfProcs = par_size();

        // get my rank
        moris_id tMyRank = par_rank();

        // communicate number of owned nodes with other procs
        Matrix< DDUMat > tFacetsOwnedPerProc;
        comm_gather_and_broadcast( aOwnedCount, tFacetsOwnedPerProc );

        // calculate node offset table
        Matrix< DDUMat > tFacetOffset( tNumberOfProcs, 1, 0 );
        for( moris_id p=1; p<tNumberOfProcs; ++p )
        {
            tFacetOffset( p ) =   tFacetOffset( p-1 ) + tFacetsOwnedPerProc( p-1 );
        }

        // remember max index for vtk
        mMaxFacetDomainIndex = tFacetOffset( tNumberOfProcs -1 ) + tFacetsOwnedPerProc( tNumberOfProcs-1 );

        moris_id tMyOffset = tFacetOffset( tMyRank );

        // update owned nodes
        for( Facet * tFacet : mFacets )
        {
            if( tFacet->get_owner() == tMyRank )
            {
                tFacet->set_id( tFacet->get_id() + tMyOffset );
            }
        }

        // step 4b: synchronize IDs for non owned facets

        // get proc neighbors from background mesh
        auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

        // get number of proc neighbors
        uint tNumberOfNeighbors = mBackgroundMesh->get_number_of_proc_neighbors();

        // create cell of matrices to send
        Matrix< DDLUMat > tEmptyLuint;
        Cell< Matrix< DDLUMat > > tAncestorListSend;
        tAncestorListSend.resize( tNumberOfNeighbors, { tEmptyLuint } );

        Matrix< DDUMat > tEmptyUint;
        Cell< Matrix< DDUMat > > tPedigreeListSend;
        tPedigreeListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

        Cell< Matrix< DDUMat > > tFacetIndexListSend;
        tFacetIndexListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

        // loop over all proc neighbors
        for ( uint p = 0; p<tNumberOfNeighbors; ++p )
        {
            auto tNeighbor = tProcNeighbors( p );

            if ( tNeighbor < tNumberOfProcs && tNeighbor != tMyRank )
            {

                // count facets that belong to neighbor
                luint tElementCounter = 0;

                // initialize counter for memory needed for pedigree tree
                luint tMemoryCounter = 0;

                // loop over all faces on this mesh
                for( Facet * tFacet : mFacets )
                {
                    if( tFacet->get_owner() == tNeighbor )
                    {
                        // increment counter
                        ++tElementCounter;

                        // get memory needed for pedigree path
                        tMemoryCounter += tFacet->get_hmr_leader()->get_background_element()
                                                                                          ->get_length_of_pedigree_path();
                    }
                }

                if ( tElementCounter > 0 )
                {
                    // prepare matrix containing ancestors
                    tAncestorListSend( p ).set_size( tElementCounter, 1 );

                    // prepare matrix containing pedigree list
                    tPedigreeListSend( p ).set_size( tMemoryCounter, 1 );

                    // prepare matrix containing face indices
                    tFacetIndexListSend( p ).set_size( tElementCounter, 1 );

                    // reset counter for elements
                    tElementCounter = 0;

                    // reset pedigree memory counter
                    tMemoryCounter = 0;

                    // loop over all faces on this mesh
                    for( Facet * tFacet : mFacets )
                    {
                        if( tFacet->get_owner() == tNeighbor )
                        {
                            // save index on leader
                            tFacetIndexListSend( p )( tElementCounter ) = tFacet->get_index_on_leader();

                            // calculate path of facet
                            tFacet->get_hmr_leader()->get_background_element()
                                                                            ->endcode_pedigree_path( tAncestorListSend( p )( tElementCounter++ ),
                                                                                    tPedigreeListSend( p ),
                                                                                    tMemoryCounter );
                        }
                    }
                }
            }
        } /* end loop over all procs */

        // initialize matrices for receiving
        Cell< Matrix< DDLUMat > > tAncestorListReceive;
        Cell< Matrix< DDUMat > >  tPedigreeListReceive;
        Cell< Matrix< DDUMat > >  tFacetIndexListReceive;

        // communicate ancestor IDs
        communicate_mats( tProcNeighbors,
                tAncestorListSend,
                tAncestorListReceive );

        // clear memory
        tAncestorListSend.clear();

        // communicate pedigree list
        communicate_mats( tProcNeighbors,
                tPedigreeListSend,
                tPedigreeListReceive );

        // clear memory
        tPedigreeListSend.clear();

        // communicate indices
        communicate_mats( tProcNeighbors,
                tFacetIndexListSend,
                tFacetIndexListReceive );

        // loop over all received lists
        for ( uint p=0; p<tNumberOfNeighbors; ++p )
        {
            // get number of elements on refinement list
            luint tNumberOfElements = tAncestorListReceive( p ).length();

            // reset memory counter
            luint tMemoryCounter = 0;

            // resize  sending list
            tFacetIndexListSend( p ).set_size( tNumberOfElements, 1 );

            // loop over all received elements
            for ( uint k=0; k<tNumberOfElements; ++k )
            {
                // decode path and get pointer to element
                Background_Element_Base * tBackElement = mBackgroundMesh->decode_pedigree_path(
                        tAncestorListReceive( p )( k ),
                        tPedigreeListReceive( p ),
                        tMemoryCounter );

                // get pointer to leader
                Element * tLeader = this->get_element_by_memory_index(
                        tBackElement->get_memory_index() );

                // get pointer to facet
                Facet * tFacet = tLeader->get_hmr_facet(
                        tFacetIndexListReceive( p )( k ) );

                // copy ID into send index
                tFacetIndexListSend( p )( k ) = tFacet->get_id();
            }
        }  /* end loop over all procs */

        // reset receive list
        tFacetIndexListReceive.clear();

        // communicate ids
        communicate_mats( tProcNeighbors,
                tFacetIndexListSend,
                tFacetIndexListReceive );

        // reset send list
        tFacetIndexListSend.clear();

        // loop over all received lists
        for ( uint p=0; p<tNumberOfNeighbors; ++p )
        {
            if( tFacetIndexListReceive( p ).length() > 0 )
            {
                // get neighbor id
                auto tNeighbor = tProcNeighbors( p );

                // reset counter
                uint tCount = 0;

                // loop over all faces on this mesh
                for( Facet * tFacet : mFacets )
                {
                    if( tFacet->get_owner() == tNeighbor )
                    {
                        // set index of facet
                        tFacet->set_id( tFacetIndexListReceive( p )( tCount++ ) );
                    }
                }

            }
        } /* end loop over all procs */
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::synchronize_edge_ids( uint aOwnedCount )
    {
        // get number of procs
        moris_id tNumberOfProcs = par_size();

        // get my rank
        moris_id tMyRank = par_rank();

        // communicate number of owned nodes with other procs
        Matrix< DDUMat > tEdgesOwnedPerProc;
        comm_gather_and_broadcast( aOwnedCount, tEdgesOwnedPerProc );

        // calculate node offset table
        Matrix< DDUMat > tEdgeOffset( tNumberOfProcs, 1, 0 );
        for( moris_id p=1; p<tNumberOfProcs; ++p )
        {
            tEdgeOffset( p ) =   tEdgeOffset( p-1 ) + tEdgesOwnedPerProc( p-1 );
        }

        mMaxEdgeDomainIndex = tEdgeOffset( tNumberOfProcs-1 )
                                                    + tEdgesOwnedPerProc( tNumberOfProcs-1 );

        moris_id tMyOffset = tEdgeOffset( tMyRank );

        // update owned nodes
        for( Edge * tEdge : mEdges )
        {
            if( tEdge->get_owner() == tMyRank )
            {
                tEdge->set_id( tEdge->get_id() + tMyOffset );
            }
        }

        // step 4b: synchronize IDs for non owned facets

        // get proc neighbors from background mesh
        auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

        // get number of proc neighbors
        uint tNumberOfNeighbors = mBackgroundMesh->get_number_of_proc_neighbors();

        // create cell of matrices to send
        Matrix< DDLUMat > tEmptyLuint;
        Cell< Matrix< DDLUMat > > tAncestorListSend;
        tAncestorListSend.resize( tNumberOfNeighbors, { tEmptyLuint } );

        Matrix< DDUMat > tEmptyUint;
        Cell< Matrix< DDUMat > > tPedigreeListSend;
        tPedigreeListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

        Cell< Matrix< DDUMat > > tEdgeIndexListSend;
        tEdgeIndexListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

        // loop over all proc neighbors
        for ( uint p = 0; p<tNumberOfNeighbors; ++p )
        {
            auto tNeighbor = tProcNeighbors( p );

            if ( tNeighbor < tNumberOfProcs && tNeighbor != tMyRank )
            {

                // count facets that belong to neighbor
                luint tElementCounter = 0;

                // initialize counter for memory needed for pedigree tree
                luint tMemoryCounter = 0;

                // loop over all faces on this mesh
                for( Edge * tEdge : mEdges )
                {
                    if( tEdge->get_owner() == tNeighbor )
                    {
                        // increment counter
                        ++tElementCounter;

                        // get memory needed for pedigree path

                        tMemoryCounter += tEdge->get_hmr_leader()
                                                        ->get_background_element()
                                                        ->get_length_of_pedigree_path();
                    }
                }

                if ( tElementCounter > 0 )
                {
                    // prepare matrix containing ancestors
                    tAncestorListSend( p ).set_size( tElementCounter, 1 );

                    // prepare matrix containing pedigree list
                    tPedigreeListSend( p ).set_size( tMemoryCounter, 1 );

                    // prepare matrix containing face indices
                    tEdgeIndexListSend( p ).set_size( tElementCounter, 1 );

                    // reset counter for elements
                    tElementCounter = 0;

                    // reset pedigree memory counter
                    tMemoryCounter = 0;

                    // loop over all edges of this mesh
                    for( Edge * tEdge : mEdges )
                    {
                        if( tEdge->get_owner() == tNeighbor )
                        {
                            // save index on leader
                            tEdgeIndexListSend( p )( tElementCounter ) = tEdge->get_index_on_leader();

                            // calculate path of facet
                            tEdge->get_hmr_leader()->get_background_element()
                                                                           ->endcode_pedigree_path(
                                                                                   tAncestorListSend( p )( tElementCounter++ ),
                                                                                   tPedigreeListSend( p ),
                                                                                   tMemoryCounter );
                        }
                    }
                }
            }
        } /* end loop over all procs */

        // initialize matrices for receiving
        Cell< Matrix< DDLUMat > > tAncestorListReceive;
        Cell< Matrix< DDUMat > >  tPedigreeListReceive;
        Cell< Matrix< DDUMat > >  tEdgeIndexListReceive;

        // communicate ancestor IDs
        communicate_mats( tProcNeighbors,
                tAncestorListSend,
                tAncestorListReceive );

        // clear memory
        tAncestorListSend.clear();

        // communicate pedigree list
        communicate_mats( tProcNeighbors,
                tPedigreeListSend,
                tPedigreeListReceive );

        // clear memory
        tPedigreeListSend.clear();

        // communicate indices
        communicate_mats( tProcNeighbors,
                tEdgeIndexListSend,
                tEdgeIndexListReceive );

        // loop over all received lists
        for ( uint p=0; p<tNumberOfNeighbors; ++p )
        {
            // get number of elements on refinement list
            luint tNumberOfElements = tAncestorListReceive( p ).length();

            // reset memory counter
            luint tMemoryCounter = 0;

            // resize  sending list
            tEdgeIndexListSend( p ).set_size( tNumberOfElements, 1 );

            // loop over all received elements
            for ( uint k=0; k<tNumberOfElements; ++k )
            {
                // decode path and get pointer to element
                Background_Element_Base* tBackElement = mBackgroundMesh->decode_pedigree_path(
                        tAncestorListReceive( p )( k ),
                        tPedigreeListReceive( p ),
                        tMemoryCounter );

                // get pointer to leader
                Element * tLeader = this->get_element_by_memory_index(
                        tBackElement->get_memory_index() );

                // get pointer to facet
                Edge * tEdge = tLeader->get_hmr_edge( tEdgeIndexListReceive( p )( k ) );

                // copy ID into send index
                tEdgeIndexListSend( p )( k ) = tEdge->get_id();
            }
        }  /* end loop over all procs */

        // reset receive list
        tEdgeIndexListReceive.clear();

        // communicate ids
        communicate_mats( tProcNeighbors,
                tEdgeIndexListSend,
                tEdgeIndexListReceive );

        // reset send list
        tEdgeIndexListSend.clear();

        // loop over all received lists
        for ( uint p=0; p<tNumberOfNeighbors; ++p )
        {
            if( tEdgeIndexListReceive( p ).length() > 0 )
            {
                // get neighbor id
                auto tNeighbor = tProcNeighbors( p );

                // reset counter
                uint tCount = 0;

                // loop over all faces on this mesh
                for( Edge * tEdge : mEdges )
                {
                    if( tEdge->get_owner() == tNeighbor )
                    {
                        // set index of facet
                        tEdge->set_id( tEdgeIndexListReceive( p )( tCount++ ) );
                    }
                }
            }
        } /* end loop over all procs */
    }

    //------------------------------------------------------------------------------

    /*       void
    Lagrange_Mesh_Base::link_facet_children_2d()
    {
        for( Facet * tFacet : mFacets )
        {
            // get leader
            Element * tLeader = tFacet->get_hmr_leader();

            if( tLeader->is_refined() )
            {
                // reserve memory for children
                tFacet->allocate_child_container( 2 );
                if( par_rank() == 0 )
                {
                    std::cout << "Leader " << tLeader->get_hmr_id() << std::endl;
                    for( uint k = 0; k<4; ++k )
                    {
                        std::cout << k << " " << ( tLeader->get_child( mAllElementsOnProc, k) != nullptr ) << std::endl;
                    }
                }
                // get index on leader
                switch( tFacet->get_index_on_leader() )
                {
                    case( 0 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 0 )->get_hmr_facet( 0 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 1 )->get_hmr_facet( 0 ), 1 );
                        break;
                    }
                    case( 1 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 1 )->get_hmr_facet( 1 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 3 )->get_hmr_facet( 1 ), 1 );
                        break;
                    }
                    case( 2 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 3 )->get_hmr_facet( 2 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 2 )->get_hmr_facet( 2 ), 1 );
                        break;
                    }
                    case( 3 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 2 )->get_hmr_facet( 3 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 0 )->get_hmr_facet( 3 ), 1 );
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "invalid child index" );
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Lagrange_Mesh_Base::link_facet_children_3d()
    {
        for( Facet * tFacet : mFacets )
        {
            // get leader
            Element * tLeader = tFacet->get_hmr_leader();

            if( tLeader->is_refined() )
            {
                // reserve memory for children
                tFacet->allocate_child_container( 4 );

                // get index on leader
                switch( tFacet->get_index_on_leader() )
                {
                    case( 0 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 0 )->get_hmr_facet( 0 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 1 )->get_hmr_facet( 0 ), 1 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 4 )->get_hmr_facet( 0 ), 2 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 5 )->get_hmr_facet( 0 ), 3 );
                        break;
                    }
                    case( 1 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 1 )->get_hmr_facet( 1 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 3 )->get_hmr_facet( 1 ), 1 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 5 )->get_hmr_facet( 1 ), 2 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 7 )->get_hmr_facet( 1 ), 3 );
                        break;
                    }
                    case( 2 ) :
                    {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 3 )->get_hmr_facet( 2 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 2 )->get_hmr_facet( 2 ), 1 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 7 )->get_hmr_facet( 2 ), 2 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 6 )->get_hmr_facet( 2 ), 3 );
                        break;
                   }
                   case( 3 ) :
                   {
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 2 )->get_hmr_facet( 3 ), 0 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 0 )->get_hmr_facet( 3 ), 1 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 6 )->get_hmr_facet( 3 ), 2 );
                        tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 4 )->get_hmr_facet( 3 ), 3 );
                        break;
                   }
                   case( 4 ) :
                   {
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 2 )->get_hmr_facet( 4 ), 0 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 3 )->get_hmr_facet( 4 ), 1 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 0 )->get_hmr_facet( 4 ), 2 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 1 )->get_hmr_facet( 4 ), 3 );
                       break;
                   }
                   case( 5 ) :
                   {
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 4 )->get_hmr_facet( 5 ), 0 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 5 )->get_hmr_facet( 5 ), 1 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 6 )->get_hmr_facet( 5 ), 2 );
                       tFacet->insert_child( tLeader->get_child( mAllElementsOnProc, 7 )->get_hmr_facet( 5 ), 3 );
                       break;
                   }
                   default :
                   {
                       MORIS_ERROR( false, "invalid child index" );
                   }
                }
            }
        }
    } */

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_faces_to_vtk( const std::string & aPath )
    {
        if( mFacets.size() > 0 )
        {
            std::string tFilePath = parallelize_path( aPath );

            // open the file
            std::ofstream tFile(tFilePath, std::ios::binary);

            // containers
            float tFChar = 0;
            int   tIChar = 0;

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;
            int tNumberOfNodes = mAllBasisOnProc.size();

            // write node data
            tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;

            // ask settings for numner of dimensions
            auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

            if ( tNumberOfDimensions == 2 )
            {
                // loop over all nodes
                for ( int k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXY = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXY[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXY[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) 0 );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                // loop over all nodes
                for ( int k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXYZ = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXYZ[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 2 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }

            tFile << std::endl;

            // get vtk index for edge
            int tCellType = 0;

            if ( mParameters->get_number_of_dimensions() == 2 )
            {
                switch( this->get_order() )
                {
                    case( 1 ) :
                                {
                        tCellType = 3;
                        break;
                                }
                    case( 2 ) :
                                {
                        tCellType = 21;
                        break;
                                }
                    case( 3 ) :
                                {
                        tCellType = 35;

                        break;
                                }
                    default :
                    {
                        tCellType = 68;
                        break;
                    }
                }
            }
            else if( mParameters->get_number_of_dimensions() == 3 )
            {
                switch( this->get_order() )
                {
                    case( 1 ) :
                                {
                        tCellType = 9;
                        break;
                                }
                    case( 2 ) :
                                {
                        tCellType = 28;
                        break;
                                }
                    default :
                    {
                        tCellType = 70;
                        break;
                    }
                }
            }

            // get number of nodes per cell
            int tNumberOfNodesPerElement = mFacets( 0 )->get_vertex_ids().length();

            // value to write in VTK file
            int tNumberOfNodesVTK = swap_byte_endian( (int) tNumberOfNodesPerElement );

            // get number of faces
            int tNumberOfElements = mFacets.size();

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                    << ( tNumberOfNodesPerElement + 1 )*tNumberOfElements  << std::endl;

            // loop over all faces
            for( Facet * tFacet : mFacets )
            {
                tFile.write( (char*) &tNumberOfNodesVTK, sizeof(int) );

                // loop over all nodes of this element
                for( int k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    // write node to mesh file
                    tIChar = swap_byte_endian( ( int ) tFacet->get_basis( k )->get_memory_index() );
                    tFile.write((char *) &tIChar, sizeof(int));
                }
            }

            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
            tIChar = swap_byte_endian( tCellType );
            for ( int k = 0; k < tNumberOfElements; ++k)
            {
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write element data
            tFile << "CELL_DATA " << tNumberOfElements << std::endl;

            // write face ID
            tFile << "SCALARS FACET_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Facet * tFacet : mFacets )
            {
                tIChar = swap_byte_endian( (int) tFacet->get_id() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write face index
            tFile << "SCALARS FACET_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Facet * tFacet : mFacets )
            {
                tIChar = swap_byte_endian( (int) tFacet->get_index() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write owner
            tFile << "SCALARS FACET_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Facet * tFacet : mFacets )
            {
                tIChar = swap_byte_endian( (int) tFacet->get_owner() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS FACET_LEVEL int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Facet * tFacet : mFacets )
            {
                tIChar = swap_byte_endian( (int) tFacet->get_hmr_leader()
                        ->get_background_element()->get_level() );

                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS FACET_STATE int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Facet * tFacet : mFacets )
            {
                if( tFacet->is_active() )
                {
                    tIChar = swap_byte_endian( (int) 1 );
                }
                else
                {
                    tIChar = swap_byte_endian( (int) 0 );
                }
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write node data
            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS NODE_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( int k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_id() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            tFile << "SCALARS NODE_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( int k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_index() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            // close the output file
            tFile.close();
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_edges_to_vtk( const std::string & aPath )
    {
        if( mEdges.size() > 0 )
        {
            std::string tFilePath = parallelize_path( aPath );

            // open the file
            std::ofstream tFile(tFilePath, std::ios::binary);

            // containers
            float tFChar = 0;
            int   tIChar = 0;

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;
            int tNumberOfNodes = mAllBasisOnProc.size();

            // write node data
            tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;

            // ask settings for numner of dimensions
            auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

            if ( tNumberOfDimensions == 2 )
            {
                // loop over all nodes
                for ( int k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXY = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXY[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXY[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) 0 );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                // loop over all nodes
                for ( int k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXYZ = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXYZ[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 2 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }

            tFile << std::endl;

            // get vtk index for edge
            int tCellType = 0;

            switch( this->get_order() )
            {
                case( 1 ) :
                        {
                    tCellType = 3;
                    break;
                        }
                case( 2 ) :
                        {
                    tCellType = 21;
                    break;
                        }
                case( 3 ) :
                        {
                    tCellType = 35;

                    break;
                        }
                default :
                {
                    tCellType = 68;
                    break;
                }
            }

            // get number of nodes per cell
            int tNumberOfNodesPerElement = mEdges( 0 )->get_vertex_ids().length();

            // value to write in VTK file
            int tNumberOfNodesVTK = swap_byte_endian( (int) tNumberOfNodesPerElement );

            // get number of faces
            int tNumberOfElements = mEdges.size();

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                    << ( tNumberOfNodesPerElement + 1 )*tNumberOfElements  << std::endl;

            // loop over all faces
            for( Edge * tEdge : mEdges )
            {
                tFile.write( (char*) &tNumberOfNodesVTK, sizeof(int) );

                // loop over all nodes of this element
                for( int k=0; k<tNumberOfNodesPerElement; ++k )
                {
                    // write node to mesh file
                    tIChar = swap_byte_endian( ( int ) tEdge->get_basis( k )->get_memory_index() );
                    tFile.write((char *) &tIChar, sizeof(int));
                }
            }

            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
            tIChar = swap_byte_endian( tCellType );
            for ( int k = 0; k < tNumberOfElements; ++k)
            {
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write element data
            tFile << "CELL_DATA " << tNumberOfElements << std::endl;

            // write element ID
            tFile << "SCALARS EDGE_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_id() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            tFile << "SCALARS EDGE_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_index() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write owner
            tFile << "SCALARS EDGE_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_owner() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS EDGE_LEVEL int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_hmr_leader()
                        ->get_background_element()->get_level() );

                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write leader
            tFile << "SCALARS EDGE_LEADER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_hmr_leader()->get_hmr_id() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write index
            tFile << "SCALARS EDGE_INDEX_ON_LEADER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for( Edge * tEdge : mEdges )
            {
                tIChar = swap_byte_endian( (int) tEdge->get_index_on_leader() );
                tFile.write( (char*) &tIChar, sizeof(int));
            }
            tFile << std::endl;

            // write level
            /* tFile << "SCALARS FACET_STATE int" << std::endl;
                    tFile << "LOOKUP_TABLE default" << std::endl;
                    for( Facet * tFacet : mFacets )
                    {
                        if( tFacet->is_active() )
                        {
                            tIChar = swap_byte_endian( (int) 1 );
                        }
                        else
                        {
                            tIChar = swap_byte_endian( (int) 0 );
                        }
                        tFile.write( (char*) &tIChar, sizeof(int));
                    }
                    tFile << std::endl; */

            // write node data
            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS NODE_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( int k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_id() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            tFile << "SCALARS NODE_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( int k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_index() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            // close the output file
            tFile.close();
        }
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::save_coeffs_to_binary_file(
            const uint          aOrder,
            const std::string & aFilePath )
    {
        // start timer
        tic tTimer;

        // make path parallel
        std::string tFilePath = parallelize_path( aFilePath );

        // Step 1: determine size of output matrix
        uint tCount = 0;

        // get number of nodes from mesh
        uint tNumberOfNodes = this->get_number_of_nodes_on_proc();

        // increment counter: first entry is number of nodes
        ++tCount;

        // loop over all nodes
        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            // get pointer to node
            mtk::Vertex * tNode = this->get_node_by_index( k );

            // increment counter for node ID, node index and number of coeffs
            // + 2*number of coefficients
            tCount += 3 + 2*tNode->get_interpolation( aOrder )->get_number_of_coefficients();
        }

        // Step 2: allocate output matrix and populate it with data

        // allocate output matrix
        Matrix< DDRMat > tOutput( tCount, 1 );

        // reset counter
        tCount = 0;

        // write number of nodes
        tOutput( tCount++ ) = tNumberOfNodes;

        // loop over all nodes
        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            // get pointer to node
            mtk::Vertex * tNode = this->get_node_by_index( k );

            // write node Index to matrix
            tOutput( tCount++ ) = tNode->get_index();

            // write node ID to matrix
            tOutput( tCount++ ) = tNode->get_id();

            // get number of coeffs
            uint tNumberOfCoeffs = tNode ->get_interpolation( aOrder )->get_number_of_coefficients();

            // write number of coeffs to matrix
            tOutput( tCount++ ) = tNumberOfCoeffs;

            // get IDs
            Matrix< IdMat > tIDs = tNode ->get_interpolation( aOrder )->get_ids();

            // get weights
            const Matrix< DDRMat > & tWeights = * tNode->get_interpolation( aOrder )->get_weights();

            // loop over all coeffs and write dof ids
            for( uint i=0; i<tNumberOfCoeffs; ++i )
            {
                tOutput( tCount++ ) =  tIDs( i );
            }

            // loop over all coeffs and write weights
            for( uint i=0; i<tNumberOfCoeffs; ++i )
            {
                tOutput( tCount++ ) =  tWeights( i );
            }
        }

        MORIS_ASSERT( tCount = tOutput.length(), "Something went wrong while writing coeffs to file." );

        // step 3: store output matrix into file
        save_matrix_to_binary_file( tOutput, tFilePath );

        // stop timer
        real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

        // Log output
        MORIS_LOG_INFO( "%s Saved coefficients to binary file:  %s.",
                proc_string().c_str(),
                tFilePath.c_str());
        MORIS_LOG_INFO( "Saving took %5.3f seconds.",
                ( double ) tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::delete_t_matrix_lagrange_mesh()
    {
        for( auto tLagrangeMesh : mLagrangeMeshForTMatrix )
        {
             delete tLagrangeMesh;
        }

        mLagrangeMeshForTMatrix.clear();
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::get_elements_in_bspline_element(
            moris_index const aBspElementIndex,
            moris_index const aDiscretizationMeshIndex,
            moris::Cell< mtk::Cell * > & aCells )
    {
        // get pointer to b-spline and background elements
        Element * tBsplineElement = mBSplineMeshes( aDiscretizationMeshIndex )->get_element( aBspElementIndex );
        Background_Element_Base * tBackgroundElement  = tBsplineElement->get_background_element();

        // check that current element is actually active on the current activation pattern
        MORIS_ERROR( tBackgroundElement->is_active( mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern() ),
            "Lagrange_Mesh_Base::get_elements_in_bspline_element() - trying to get non-active B-Spline element." );

        // get the number of active Lagrange elmements within active B-spline element
        luint tNumActiveLagrangeElements = 0;
        uint tLagrangePattern = this->get_activation_pattern();
        tBackgroundElement->get_number_of_active_descendants( tLagrangePattern, tNumActiveLagrangeElements );

        // collect the active Lagrange elements
        tNumActiveLagrangeElements = 0;
        moris::Cell< Background_Element_Base * > tActiveElements( tNumActiveLagrangeElements, nullptr );
        tBackgroundElement->collect_active_descendants( tLagrangePattern, tActiveElements, tNumActiveLagrangeElements );

        // initialize output cell with correct size
        aCells.resize( tNumActiveLagrangeElements, nullptr );

        //
        for( uint iLagElem = 0; iLagElem < tNumActiveLagrangeElements; iLagElem ++ )
        {
            luint tMemoryIndex = tActiveElements( iLagElem )->get_memory_index();
            aCells( iLagElem ) = this->get_element_by_memory_index( tMemoryIndex );
        }
    }

    //------------------------------------------------------------------------------

    void
    Lagrange_Mesh_Base::get_lagrange_elements_in_bspline_elements(
            moris_index const                          aDiscretizationMeshIndex,
            moris::Cell< moris::Cell< mtk::Cell* > >&  aCells,
            moris::Cell< moris::Cell< moris_index > >& aCellIndices,
            moris::Cell< moris_index >&                aLagToBspCellIndices,
            moris::Cell< uint >&                       aBspCellRefineLevels )
    {
        // get B-Spline pattern of this mesh
        uint tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

        // get Lagrange pattern of this mesh
        uint tLagrangePattern = this->get_activation_pattern();

        // get number of Lagrange elements on current mesh
        luint tNumLagElems = this->get_background_mesh()->get_number_of_active_elements_on_proc_including_aura();
        aLagToBspCellIndices.resize( tNumLagElems, -1 );

        // set the activation pattern to that of the B-spline mesh
        MORIS_LOG_INFO( "Set activation pattern to B-spline mesh to retrieve B-spline information" );
        mBackgroundMesh->set_activation_pattern( tBSplinePattern );

        // get number of B-spline elements on bg mesh
        luint tNumBsplineElems = this->get_num_active_bg_elements_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex );

        // get indices of active B-spline elements on bg mesh
        Matrix< DDLUMat > tElementIndices;
        this->get_active_bg_element_indices_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex, tElementIndices );

        // check for debug
        MORIS_ASSERT( tElementIndices.numel() == tNumBsplineElems,
            "Lagrange_Mesh_Base::get_lagrange_elements_in_bspline_elements() - number of Bsp-elems and size of Bsp-elem index list don't match." );

        // initialize the output list sizes
        aCells.resize( tNumBsplineElems );
        aCellIndices.resize( tNumBsplineElems );
        aBspCellRefineLevels.resize( tNumBsplineElems );

        // for each B-spline element find and store the Lagrange elements within it
        for ( luint iBspElem = 0; iBspElem < tNumBsplineElems; iBspElem++ )
        {
            // get pointer to b-spline and background elements
            Element* tBsplineElement = mBSplineMeshes( aDiscretizationMeshIndex )->get_element_including_aura( iBspElem );
            Background_Element_Base* tBackgroundElement  = tBsplineElement->get_background_element();

            // get and store the refinement level of the current B-spline element
            aBspCellRefineLevels( iBspElem ) = tBsplineElement->get_level();

            // check that current element is actually active on the current activation pattern
            MORIS_ASSERT( tBackgroundElement->is_active( mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern() ),
                "Lagrange_Mesh_Base::get_elements_in_bspline_element() - trying to get non-active B-Spline element." );

            // get the number of active Lagrange elmements within active B-spline element
            luint tNumActiveLagrangeElements = 0;
            tBackgroundElement->get_number_of_active_descendants( tLagrangePattern, tNumActiveLagrangeElements );

            // collect the active Lagrange elements
            luint tNumActiveLagrangeElementsCheck = 0;
            moris::Cell< Background_Element_Base * > tActiveElements( tNumActiveLagrangeElements, nullptr );
            tBackgroundElement->collect_active_descendants( tLagrangePattern, tActiveElements, tNumActiveLagrangeElementsCheck );

            // sanity check that the number of descendants in list matches number reported
            MORIS_ASSERT( tNumActiveLagrangeElementsCheck == tNumActiveLagrangeElements,
                "Lagrange_Mesh_Base::get_lagrange_elements_in_bspline_elements() - Number of Active descendants doesn't match list of descendants." );

            // initialize output cell with correct size
            aCells( iBspElem ).resize( tNumActiveLagrangeElements, nullptr );
            aCellIndices( iBspElem ).resize( tNumActiveLagrangeElements, -1 );

            // get the Lagrange elements and their indices
            for( uint iLagElem = 0; iLagElem < tNumActiveLagrangeElements; iLagElem ++ )
            {
                // find Lagrange elments via memory index and store them in output lists
                luint tMemoryIndex = tActiveElements( iLagElem )->get_memory_index();
                aCells( iBspElem )( iLagElem ) = this->get_element_by_memory_index( tMemoryIndex );
                moris_index tLagElemIndex = aCells( iBspElem )( iLagElem )->get_index();
                aCellIndices( iBspElem )( iLagElem ) = tLagElemIndex;

                // store which B-spline element current Larange element belongs to
                aLagToBspCellIndices( (uint) tLagElemIndex ) = iBspElem;
            }
        }

        // set the activation pattern back to the pattern of the Lagrange mesh
        MORIS_LOG_INFO( "Reset activation pattern back to Lagrange mesh" );
        mBackgroundMesh->set_activation_pattern( tLagrangePattern );
    }

    //------------------------------------------------------------------------------

    const luint*
    Lagrange_Mesh_Base::get_bspline_element_ijk_level(
            moris_index         aDiscretizationMeshIndex,
            moris_index         aBsplineElementIndex,
            uint                aLevel )
    {
        // get B-Spline pattern of this mesh
        uint tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

        // get Lagrange pattern of this mesh
        uint tLagrangePattern = this->get_activation_pattern();

        // set the activation pattern to that of the B-spline mesh
        mBackgroundMesh->set_activation_pattern( tBSplinePattern );

        // use the b-spline mesh index to get the correct b-spline element
        Element*                 tBsplineElement    = mBSplineMeshes( aDiscretizationMeshIndex )->get_element_including_aura( aBsplineElementIndex );
        aLevel = tBsplineElement->get_level(); // FIXME what is going on here?

        // set the activation pattern back to the pattern of the Lagrange mesh
        mBackgroundMesh->set_activation_pattern( tLagrangePattern );

        // get the i-j-k indices of the b-spline element
        return tBsplineElement->get_ijk();
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::get_elements_in_interpolation_cluster(
            moris_index const aElementIndex,
            moris_index const aDiscretizationMeshIndex,
            moris::Cell< mtk::Cell * > & aCells)
    {
        // get B-Spline pattern of this mesh
        auto tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

        // get Lagrange pattern of this mesh
        auto tLagrangePattern = this->get_activation_pattern();

        // select pattern
        //this->select_activation_pattern();

        // get pointer to element
        auto tLagrangeElement = this->get_element_including_aura( aElementIndex );

        // get pointer to background element
        auto tBackgroundElement = tLagrangeElement->get_background_element();

        // keep jumping to parent until we have the active element
        // Note: underlying assumption that the Lagrange element is at least as refined as the active bspline element
        while( ! tBackgroundElement->is_active( tBSplinePattern ) )
        {
            // jump to parent
            tBackgroundElement = tBackgroundElement->get_parent();
        }

        // initialize counter
        luint tCount = 0;

        tBackgroundElement->get_number_of_active_descendants( tLagrangePattern, tCount );

        moris::Cell< Background_Element_Base * > tActiveElements( tCount, nullptr );

        // reset counter
        tCount = 0;

        tBackgroundElement->collect_active_descendants( tLagrangePattern, tActiveElements, tCount );

        aCells.resize( tCount, nullptr );

        for( uint Ik = 0; Ik < tCount; Ik ++ )
        {
            luint tMemoryIndex = tActiveElements( Ik )->get_memory_index();

            aCells( Ik ) = this->get_element_by_memory_index( tMemoryIndex );
        }
    }

    //------------------------------------------------------------------------------

    void
    Lagrange_Mesh_Base::get_elements_in_bspline_element_and_side_ordinal(
            moris_index const          aBsplineElementIndex,
            moris_index const          aDiscretizationMeshIndex,
            moris_index const          aSideOrdinal,
            moris::Cell< mtk::Cell* >& aCells )
    {
        // get pattern of the current B-spline mesh
        uint tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

        // get Lagrange pattern of this mesh
        uint tLagrangePattern = this->get_activation_pattern();

        // set activation pattern to B-splines
        mBackgroundMesh->set_activation_pattern( tBSplinePattern );

        // get pointer to B-spline and background element
        Element * tBsplineElement = mBSplineMeshes( aDiscretizationMeshIndex )->get_element_including_aura( aBsplineElementIndex );

        // get pointer to background element
        Background_Element_Base* tBackgroundElement = tBsplineElement->get_background_element();

        // // since the currently active pattern may be finer, jump back to parent elements until level associated with B-spline mesh is found
        // while( ! tBackgroundElement->is_active( tBSplinePattern ) )
        // {
        //     // jump to parent
        //     tBackgroundElement = tBackgroundElement->get_parent();
        // }

        // initialize variable storing number of Lagrange elements on side ordinal
        luint tCount = 0;

        // retrieve number of lagrange elements on side ordinal
        switch( aSideOrdinal )
        {
            case( 1 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_1( tLagrangePattern, tCount );
                break;
            }
            case( 2 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_2( tLagrangePattern, tCount );
                break;
            }
            case( 3 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_3( tLagrangePattern, tCount );
                break;
            }
            case( 4 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_4( tLagrangePattern, tCount );
                break;
            }
            case( 5 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_5( tLagrangePattern, tCount );
                break;
            }
            case( 6 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_6( tLagrangePattern, tCount );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid Side set ordinal.");
            }
        }

        // initialize list of background elements on side ordinal
        moris::Cell< Background_Element_Base * > tActiveElements( tCount, nullptr );

        // reset counter
        tCount = 0;

        // retrieve lagrange elements on side ordinal
        switch( aSideOrdinal )
        {
            case( 1 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_1( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 2 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_2( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 3 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_3( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 4 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_4( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 5 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_5( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 6 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_6( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid Side set ordinal.");
            }
        }

        // initialize list of mtk::cells on side ordinal
        aCells.resize( tCount, nullptr );

        // retrieve corresponding IP elements from corresponding background elements
        for( uint Ik = 0; Ik < tCount; Ik ++ )
        {
            luint tMemoryIndex = tActiveElements( Ik )->get_memory_index();
            aCells( Ik ) = this->get_element_by_memory_index( tMemoryIndex );
        }

        // set activation pattern back to the Lagrange
        mBackgroundMesh->set_activation_pattern( tLagrangePattern );
    }

    //------------------------------------------------------------------------------

    void Lagrange_Mesh_Base::get_elements_in_interpolation_cluster_and_side_ordinal(
            moris_index const            aElementIndex,
            moris_index const            aDiscretizationMeshIndex,
            moris_index const            aSideOrdinal,
            moris::Cell< mtk::Cell * > & aCells )
    {
        // get B-Spline pattern of this mesh
        auto tBSplinePattern = mBSplineMeshes( aDiscretizationMeshIndex )->get_activation_pattern();

        // get Lagrange pattern of this mesh
        auto tLagrangePattern = this->get_activation_pattern();

        // select pattern
        //this->select_activation_pattern();

        // get pointer to element
        auto tLagrangeElement = this->get_element_including_aura( aElementIndex );

        // get pointer to background element
        auto tBackgroundElement = tLagrangeElement->get_background_element();

        while( ! tBackgroundElement->is_active( tBSplinePattern ) )
        {
            // jump to parent
            tBackgroundElement = tBackgroundElement->get_parent();
        }

        // initialize counter
        luint tCount = 0;

        switch( aSideOrdinal )
        {
            case( 1 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_1( tLagrangePattern, tCount );
                break;
            }
            case( 2 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_2( tLagrangePattern, tCount );
                break;
            }
            case( 3 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_3( tLagrangePattern, tCount );
                break;
            }
            case( 4 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_4( tLagrangePattern, tCount );
                break;
            }
            case( 5 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_5( tLagrangePattern, tCount );
                break;
            }
            case( 6 ) :
            {
                tBackgroundElement->get_number_of_active_descendants_on_side_6( tLagrangePattern, tCount );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid Side set ordinal.");
            }
        }

        moris::Cell< Background_Element_Base * > tActiveElements( tCount, nullptr );

        // reset counter
        tCount = 0;

        switch( aSideOrdinal )
        {
            case( 1 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_1( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 2 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_2( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 3 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_3( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 4 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_4( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 5 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_5( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            case( 6 ) :
            {
                tBackgroundElement->collect_active_descendants_on_side_6( tLagrangePattern, tActiveElements, tCount );
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Invalid Side set ordinal.");
            }
        }

        aCells.resize( tCount, nullptr );

        for( uint Ik = 0; Ik < tCount; Ik ++ )
        {
            luint tMemoryIndex =tActiveElements( Ik )->get_memory_index();

            aCells( Ik ) = this->get_element_by_memory_index( tMemoryIndex );
        }
    }

    //------------------------------------------------------------------------------

    // BIG HACK for femdoc with explicit consent of Kurt. Only tested in serial and linear meshes.
    void Lagrange_Mesh_Base::nodes_renumbering_hack_for_femdoc()
    {
        MORIS_ERROR( par_size() <= 1, "Lagrange_Mesh_Base::nodes_renumbering_hack_for_femdoc(), this function is intended to work only in serial");

        moris::uint tCounter  = 0;
        moris::uint tCounter2 = 0;

        moris::sint tMaxID = 0;

        uint tNumberOfNodes = this->get_number_of_nodes_on_proc();

        for( uint Ik = 0; Ik<tNumberOfNodes; Ik ++ )
        {
            Basis * tBasis = this->get_node_by_index( Ik );

            MORIS_ERROR( tBasis->has_interpolation (1), "Lagrange_Mesh_Base:: node has no first order basis");

            mtk::Vertex_Interpolation * tInterp = tBasis->get_interpolation( 1 );

            Matrix< IdMat > tLocalIDs = tInterp->get_ids();

            tMaxID= std::max( tMaxID, tLocalIDs.max() );
        }

        Matrix< DDSMat > tReverseIndexMap( tMaxID+1, 1, -1 );
        Matrix< DDSMat > tReverseIDMap( tMaxID+1, 1, -1 );

        moris::Cell< Basis * >tNonBSplineBasis( mAllBasisOnProc.size(), nullptr );

        this->calculate_t_matrices( false );

        for( uint Ik = 0; Ik<tNumberOfNodes; Ik ++ )
        {
            Basis * tBasis = this->get_node_by_index( Ik );

            MORIS_ERROR( tBasis->has_interpolation (1), "Lagrange_Mesh_Base:: node has no first order basis");

            mtk::Vertex_Interpolation * tInterp = tBasis->get_interpolation( 1 );

            Matrix< IdMat > tLocalIDs = tInterp->get_ids();

            const Matrix< DDRMat > & tLocalWeights = *tInterp->get_weights();

            if ( tLocalIDs.numel()==1 || equal_to(tLocalWeights.max(),1.0) )
            {
                uint tIndex1 = 0;

                // value 1 should be always on the first entry in the vector
                // FIXME: should be replace by proper moris::mat function
                for( uint Ia = 0; Ia<tLocalWeights.numel(); Ia ++ )
                {
                    if( equal_to(tLocalWeights(Ia),1.0) )
                    {
                        tIndex1 = Ia;
                        break;
                    }
                }

                // check whether the same basis is used twice for being the only basis interpolating at a node
                MORIS_ASSERT( tReverseIDMap( tLocalIDs( tIndex1, 0 )-1 ) == -1, "Node Id %-5i appears twice", tLocalIDs( tIndex1, 0 )-1 );

                tReverseIndexMap( tLocalIDs( tIndex1, 0 ) ) = tBasis->get_index();
                tReverseIDMap( tLocalIDs( tIndex1, 0 )-1 )  = tBasis->get_hmr_index();

                tBasis->set_local_index( tLocalIDs( tIndex1, 0 ) );
                tBasis->set_domain_index( tLocalIDs( tIndex1, 0 ) - 1 );

                tCounter++;
            }

            else
            {
                tNonBSplineBasis( tCounter2++ ) = tBasis;
            }

        }
        tNonBSplineBasis.resize( tCounter2 );

        //            uint tMaxID = tReverseIDMap.max();

        tReverseIndexMap.resize( tMaxID + tNonBSplineBasis.size()+1, 1 );
        tReverseIDMap.resize( tMaxID + tNonBSplineBasis.size()+1, 1 );

        for( Basis * tBasis : tNonBSplineBasis )
        {
            tReverseIndexMap( tMaxID ) = tBasis->get_index();
            tReverseIDMap( tMaxID )    = tBasis->get_hmr_index();

            tBasis->set_local_index( tMaxID );
            tBasis->set_domain_index( tMaxID );

            tCounter++;
            tMaxID++;
        }

        // write reverse map for linear mesh
        std::string tFilePath = "Reverse_Map_1.hdf5";

        // make path parallel
        tFilePath = parallelize_path( tFilePath );

        // Create a new file using default properties
        hid_t tFileID = H5Fcreate( tFilePath.c_str(),
                H5F_ACC_TRUNC,
                H5P_DEFAULT,
                H5P_DEFAULT);

        // error handler
        herr_t tStatus;

        // save reverse indices to file
        save_matrix_to_hdf5_file( tFileID,
                "Index",
                tReverseIndexMap,
                tStatus );

        // save reverse ids to file
        save_matrix_to_hdf5_file( tFileID,
                "Id",
                tReverseIDMap,
                tStatus );
    }

} /* namespace moris */
