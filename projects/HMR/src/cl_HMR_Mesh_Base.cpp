/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh_Base.cpp
 *
 */

#include "cl_HMR_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
        //------------------------------------------------------------------------------

        Mesh_Base::Mesh_Base (
                const Parameters           * aParameters,
                Background_Mesh_Base       * aBackgroundMesh,
                uint aOrder,
                uint aActivationPattern )
        : mParameters( aParameters ),
          mBackgroundMesh( aBackgroundMesh ),
          mOrder( aOrder ),
          mNumberOfDimensions( mParameters->get_number_of_dimensions() ),
          mNumberOfBasisPerElement( pow( aOrder+1, mParameters->get_number_of_dimensions() ) ),
          mNumberOfNeighborsPerElement( pow( 3, mParameters->get_number_of_dimensions() ) - 1 ),
          mActivationPattern( aActivationPattern )
        {
        }
        //------------------------------------------------------------------------------
        // protected:
        // -----------------------------------------------------------------------------

        void Mesh_Base::delete_pointers()
        {
            if ( mAllBasisOnProc.size() > 0 )
            {
                // delete all nodes
                for ( auto tBasis: mAllBasisOnProc )
                {
                    delete tBasis;
                }

                // clear cell
                mAllBasisOnProc.clear();
            }

            if ( mAllElementsOnProc.size() > 0 )
            {
                // delete all elements
                for ( auto tElement: mAllElementsOnProc )
                {
                    delete tElement;
                }

                // clear cell
                mAllElementsOnProc.clear();
            }

            // reset node counter
            mNumberOfAllBasis = 0;

            // reset element counter
            mNumberOfElements = 0;
        }

        // -----------------------------------------------------------------------------

        void Mesh_Base::create_elements()
        {
            // cell for background elements
            Cell< Background_Element_Base* > tAllBackgroundElements;

            // collect all background elements on proc
            mBackgroundMesh->collect_all_elements( tAllBackgroundElements );

            // set size of all elements
            mNumberOfAllElementsOnProc = tAllBackgroundElements.size();

            // reserve memory
            mAllElementsOnProc.resize( mNumberOfAllElementsOnProc, nullptr );

            // loop over all background elements
            for ( luint k = 0; k < mNumberOfAllElementsOnProc; ++k )
            {
                mAllElementsOnProc( k ) = this->create_element( tAllBackgroundElements( k ) );
            }

            this->collect_coarsest_elements();
        }

        // -----------------------------------------------------------------------------

        /**
         * collects all coarsest elements on proc including aura
         * Memory index of background elements and all B-spline elements is the same.
         *
         * @return void
         */
        void Mesh_Base::collect_coarsest_elements()
        {
            // count number of coarsest elements
            luint tNumberOfElements = mBackgroundMesh ->get_number_of_coarsest_elements_on_proc_including_aura();

            // reset cell
            mAllCoarsestElementsOnProc.clear();

            // assign memory for cell
            mAllCoarsestElementsOnProc.resize( tNumberOfElements, nullptr );

            // loop over all coarsest elements
            for( luint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to background mesh element
                Background_Element_Base* tBackElement = mBackgroundMesh->get_coarsest_element_by_subdomain_id( e );

                // copy pointer into cell of coarsest elements
                mAllCoarsestElementsOnProc( e ) = mAllElementsOnProc( tBackElement->get_memory_index() );
            }
        }

        // -----------------------------------------------------------------------------

        void Mesh_Base::unflag_all_basis()
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                tBasis->unflag();
            }
        }

        // -----------------------------------------------------------------------------

        void Mesh_Base::unuse_all_basis()
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                tBasis->unuse();
            }
        }

        // -----------------------------------------------------------------------------

        /**
         * sets the flag of all nodes on proc
         */
        void Mesh_Base::flag_all_basis()
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                tBasis->flag();
            }
        }

        // -----------------------------------------------------------------------------

        /**
         * Determines which nodes are connected to which element.
         * Writes value in mesh data struct
         * @return void
         */
        void Mesh_Base::determine_elements_connected_to_basis()
        {
            // ask background mesh about active elements in proc
            luint tNumberOfElements = mBackgroundMesh ->get_number_of_active_elements_on_proc_including_aura();

            // loop over all background elements
            for ( luint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to background element
                Background_Element_Base* tBackElement = mBackgroundMesh->get_element_from_proc_domain_including_aura( e );

                // get pointer to Lagrange element
                Element * tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                // loop over all nodes connected to element
                for( uint k = 0; k < mNumberOfBasisPerElement; ++k  )
                {
                    // get node pointer
                    Basis* tNode = tElement->get_basis( k );

                    // increment element counter
                    tNode->increment_element_counter();
                }
            }

            // loop over all nodes
            for ( auto tBasis : mAllBasisOnProc )
            {
                // reserve memory for node container and reset element counter
                tBasis->init_element_container() ;
            }

            // loop over all background elements
            for ( luint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to background element
                Background_Element_Base* tBackElement = mBackgroundMesh->get_element_from_proc_domain_including_aura( e );

                // get pointer to Lagrange element
                Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                // loop over all nodes connected to element
                for( uint k = 0; k < mNumberOfBasisPerElement; ++k  )
                {
                    // get node pointer
                    Basis* tNode = tElement->get_basis( k );

                    // link this element with the node
                    tNode->insert_element( tElement );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Mesh_Base::guess_basis_ownership()
        {
            // only need to do something in parallel mode
            if( par_size() > 1 )
            {
                // loop over all basis
                for ( auto tBasis : mAllBasisOnProc )
                {
                    // get number of connected elements
                    uint tNumberOfElements = tBasis->get_element_counter();

                    // loop over all connected elements
                    for( uint k = 0; k < tNumberOfElements; ++k )
                    {
                        // get pointer to element
                        Element * tElement = tBasis->get_element( k );

                        // get owner of element
                        auto tElementOwner = tElement->get_owner();

                        // basis belongs to proc with smallest ID
                        if ( tBasis->get_owner() > tElementOwner )
                        {
                            tBasis->set_owner( tElementOwner );
                        }
                    }
                }
            }
            else
            {
                // in serial, claim ownership of all basis
                moris_id tMyRank = par_rank();
                for ( auto tBasis : mAllBasisOnProc )
                {
                    tBasis->set_owner( tMyRank );
                }
            }
        }

        //------------------------------------------------------------------------------

        void Mesh_Base::confirm_basis_ownership()
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
                        // cell containing basis pointers
                        Cell< Basis* > tBasisInAura;

                        // collect basis within inverse aura
                        this->collect_basis_from_aura( p, 2, tBasisInAura );

                        // calculate addresses of basis to ask for
                        this->encode_foreign_basis_path( tBasisInAura,
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
                Cell< Matrix<  DDUMat > > tSendOwner( tNumberOfProcNeighbors, tEmptyUint );

                // loop over all proc neighbors
                for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
                {
                    // get number of basis requested by neighbor
                    luint tNumberOfBasis = tReceiveBasisIndex( p ).length();

                    // initialize matrix
                    tSendOwner( p ).set_size( tNumberOfBasis, 1 );

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
                        tSendOwner( p )( k ) = tBasis->get_owner();
                    }
                }

                // free memory
                tReceiveBasisIndex.clear();
                tReceiveAncestor  .clear();
                tReceivePedigree  .clear();

                // communicate owners
                Cell< Matrix<  DDUMat > > tReceiveOwner;

                communicate_mats(
                        tProcNeighbors,
                        tSendOwner,
                        tReceiveOwner );

                // free memory
                tSendOwner.clear();

                // loop over all proc neighbors
                for ( uint p = 0; p < tNumberOfProcNeighbors; ++p )
                {
                    // get rank of neighbor
                    moris_id tNeihgborRank = tProcNeighbors( p );

                    if ( tNeihgborRank != tMyRank && tNeihgborRank != gNoProcNeighbor )
                    {
                        // cell containing basis pointers
                        Cell< Basis* > tBasisInAura;

                        // collect basis within inverse aura
                        this->collect_basis_from_aura( p, 2, tBasisInAura );

                        // initialize counter
                        luint tCount = 0;

                        // count number of basis suspected to be owned by neighbor
                        for( auto tBasis : tBasisInAura )
                        {
                            if ( tBasis->get_owner() == tNeihgborRank )
                            {
                                // set ownership from received matrix
                                tBasis->set_owner( tReceiveOwner( p )( tCount++ ) );
                            }
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------

        void Mesh_Base::collect_basis_from_aura(
                uint            aProcNeighborIndex,
                uint            aMode,
                Cell< Basis* >& aBasisList )
        {
            // clear basis list
            aBasisList.clear();

            // only do something in parallel mode
            if ( par_size() > 1 )
            {
                // cell of elements on background mesh
                Cell< Background_Element_Base* > tBackElements;

                // get element list from background mesh
                mBackgroundMesh->collect_active_elements_from_aura( aProcNeighborIndex,
                        aMode,
                        tBackElements );

                // initialize basis counter
                luint tBasisCount = 0;

                // loop over all background elements
                for ( auto tBackElement : tBackElements )
                {
                    // get pointer to Lagrange element
                    Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                    // loop over all basis
                    for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

                        // test if basis exists
                        if ( tBasis != NULL )
                        {
                            // test if basis is not flagged
                            if ( ! tBasis->is_flagged() )
                            {
                                // increment counter
                                ++tBasisCount;

                                // flag this basis
                                tBasis->flag();
                            }
                        }
                    }
                }

                // reserve memory for output cell
                aBasisList.resize( tBasisCount, nullptr );

                // reset basis counter
                tBasisCount = 0;

                // loop over all background elements
                for ( auto tBackElement : tBackElements )
                {
                    // get pointer to Lagrange element
                    Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                    // loop over all basis
                    for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

                        // test if basis exists
                        if ( tBasis != NULL )
                        {
                            // test if basis is flagged
                            if ( tBasis->is_flagged() )
                            {
                                // increment counter
                                aBasisList( tBasisCount++ ) = tBasis;

                                // unflag this basis
                                tBasis->unflag();
                            }
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------------------

        void
        Mesh_Base::encode_foreign_basis_path(
                Cell< Basis* >    & aBasis,
                const moris_id    & aOwner,
                Matrix< DDLUMat > & aElementAncestors,
                Matrix< DDUMat >  & aElementPedigree,
                Matrix< DDUMat >  & aElementLocalIndex )
        {
            // initialize counter
            luint tCount = 0;

            // count basis belonging to owner
            for( auto tBasis: aBasis )
            {
                // test if basis belongs to defined owner
                if ( tBasis->get_owner() == aOwner )
                {
                    // increment counter
                    ++tCount;
                }
            }

            // allocate matrix for local index
            aElementLocalIndex.set_size( tCount, 1 );

            // allocate matrices for element ancestor
            // contains element memory indes first and is overwritten later
            aElementAncestors.set_size( tCount, 1 );

            // reset counter
            tCount = 0;

            // loop over all basis and determine reference elements
            for( auto tBasis: aBasis )
            {
                // test if basis belongs to defined owner
                if ( tBasis->get_owner() == aOwner )
                {
                    this->get_reference_element_of_basis( tBasis,
                            aElementAncestors( tCount ),
                            aElementLocalIndex( tCount ) );

                    // increment counter
                    ++tCount;
                }
            }

            // initialize memory counter
            luint tMemoryCount = 0;

            // determine memory needs for pedigree list
            for( luint k=0; k<tCount; ++k )
            {
                tMemoryCount +=
                        mAllElementsOnProc( aElementAncestors( k ) )->
                        get_background_element()->get_length_of_pedigree_path();
            }

            // allocate pedigree path
            aElementPedigree.set_size( tMemoryCount, 1 );

            // reset counter
            tMemoryCount = 0;

            // encode pedigree path
            for( luint k = 0; k < tCount; ++k )
            {
                // get pointer to element
                Background_Element_Base* tElement =
                        mAllElementsOnProc( aElementAncestors( k ) )->get_background_element();

                // encode path and overwrite aElementAncestor with Ancestor Index
                tElement->endcode_pedigree_path(
                        aElementAncestors( k ),
                        aElementPedigree,
                        tMemoryCount );
            }
        }

        //------------------------------------------------------------------------------

        void
        Mesh_Base::get_reference_element_of_basis(
                Basis * aBasis,
                luint & aElementMemoryIndex,
                uint  & aElementLocalBasisIndex )
        {
            // copy owner of element into temporary variable
            moris_id tOwner = aBasis->get_owner();

            // find out how many elements are connected to this basis
            uint tNumberOfElements = aBasis->get_element_counter();

            // pointer to background element
            Element* tElement = nullptr;

            // loop over all elements of this basis
            for( uint e = 0; e < tNumberOfElements; ++ e )
            {
                // get pointer to element of this basis
                tElement = aBasis->get_element( e );

                // test if element is active and has same owner
                if ( tElement->is_active() && tElement->get_owner() == tOwner )
                {
                    // exit the loop
                    break;
                }
            }

            // get my memory index
            uint tID = aBasis->get_hmr_id();

            // set basis to invalid value
            aElementLocalBasisIndex = mNumberOfBasisPerElement;

            // find out local node id on element
            for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
            {
                // get pointer to basis
                Basis* tBasis = tElement->get_basis( k );

                // test if basis exists
                if( tBasis != NULL )
                {
                    if ( tBasis->get_hmr_id() == tID )
                    {
                        // copy local ID
                        aElementLocalBasisIndex = k;

                        // break loop
                        break;
                    }
                }
            }

            MORIS_ASSERT( aElementLocalBasisIndex != mNumberOfBasisPerElement,
                    "something went wrong in Mesh_Base::get_reference_element_of_basis()" );

            // get pointer to background element
            aElementMemoryIndex = tElement->get_background_element()->get_memory_index();
        }

        // -----------------------------------------------------------------------------

        void
        Mesh_Base::get_basis_coords_of_element(
                Matrix< DDRMat >  & aBasisCoords,
                luint               aElementIndex )
        {
            // set number of basis per element
            uint tNumberOfBasisPerElement;

            // test if serendipity flag is set
            /*if ( aSerendipityFlag )
            {
               // copy only serendipity basis
                tNumberOfBasisPerElement =
                  this->get_number_of_serendipity_basis();
            }
            else
            { */
            // copy all basis
            tNumberOfBasisPerElement = mNumberOfBasisPerElement;
            //}

            // set size of output matrix
            aBasisCoords.set_size( tNumberOfBasisPerElement,
                    mNumberOfDimensions );

            // get pointer to element
            Element* tElement = this->get_element( aElementIndex );

            // loop over all nodes of this element
            for( uint k = 0; k < tNumberOfBasisPerElement; ++k )
            {
                // get pointer to basis
                Basis* tBasis = tElement->get_basis( k );

                // get array of node coordinates
                const real * tXYZ = tBasis->get_xyz();

                // write coordinates into output matrix
                for( uint i = 0; i < mNumberOfDimensions; ++i )
                {
                    aBasisCoords( k, i ) = tXYZ[ i ];
                }
            }
        }

        // -----------------------------------------------------------------------------

        /*uint
        Mesh_Base::get_number_of_serendipity_basis()
        {
            switch ( mNumberOfBasisPerElement )
            {

                case( 9 ) : // QUAD9 -> QUAD8
                {
                    return 8;
                    break;
                }
                case( 27 ) :  // HEX27 -> HEX20
                {
                    return 20;
                    break;
                }
                //  the following two would work in theory,
                //  however, the interpolation functions do not exist in MORIS
                //
                //  case( 16 ) : // QUAD16 -> QUAD12
                //{
                //    return 12;
                //    break;
                //}
                //case( 64 ) : // HEX64 -> HEX32
                //{
                //   return 32;
                //    break;
                //}
                default :
                {

                    MORIS_ERROR( false, "unknown number of nodes for serendipity");
                    return 0;
                    break;
                }
            }

        } */

        // -----------------------------------------------------------------------------

        void
        Mesh_Base::update_element_indices()
        {
            // select activation pattern
            this->select_activation_pattern();

            // get number of active elements
            if( mParameters->use_number_aura() and mParameters->is_output_mesh( mMeshIndex ) )
            {
                moris_index tNumberOfElements = this->get_number_of_elements_including_aura();

                // loop over all active elements
                for ( moris_index e = 0; e < tNumberOfElements; ++e )
                {
                    // write index into active element
                    this->get_element_including_aura( e )->set_index( e );
                }
            }
            else
            {
                moris_index tNumberOfElements = this->get_number_of_elements();

                // loop over all active elements
                for ( moris_index e = 0; e < tNumberOfElements; ++e )
                {
                    // write index into active element
                    this->get_element( e )->set_index( e );
                }
            }
        }

        // -----------------------------------------------------------------------------

        moris::luint Mesh_Base::get_max_basis_hmr_id()
        {
            luint tHMRID = 0;
            for( auto tBasis : mAllBasisOnProc )
            {
                tHMRID = std::max(  tHMRID, (luint)tBasis->get_hmr_id() );
            }

            luint tHMRID_return = max_all( tHMRID );

            return tHMRID_return;
        }

        // -----------------------------------------------------------------------------

        void Mesh_Base::sanity_check_for_ids_and_ownership()
        {
            luint tMax_HMR_ID = this->get_max_basis_hmr_id();

            // create Map
            std::map< luint, luint> tIdMap;
            std::map< luint, moris_id > tOwnerMap;
            for( auto tBasis : mAllBasisOnProc )
            {
                tIdMap[ tBasis->get_hmr_id() ] = tBasis->get_hmr_index();
                tOwnerMap[ tBasis->get_hmr_id() ] = tBasis->get_owner();
            }

            for( int Ik = 0; Ik < (int)tMax_HMR_ID; Ik++ )
            {
                luint tId = gNoEntityID;
                moris_id tOwner = gNoProcID;

                auto tIter = tIdMap.find(Ik);

                if ( tIter!=tIdMap.end()  )
                {
                    tId = tIter->second;
                }

                auto tIterOwner = tOwnerMap.find(Ik);

                if ( tIterOwner!=tOwnerMap.end()  )
                {
                    tOwner = tIterOwner->second;
                }

                sint tNumProcs = par_size();
                moris::Cell< luint > tReciveBuffer(tNumProcs);
                moris::Cell< moris_id > tReciveBufferOwner(tNumProcs);

                moris::Cell< sint > tReciveBufferNumber( tNumProcs, 1 );
                moris::Cell< sint > tReciveBufferOFFSET( tNumProcs, 0 );

                for( int Ii = 0; Ii < tNumProcs; Ii++ )
                {
                    tReciveBufferOFFSET( Ii) = Ii;
                }

                MPI_Gatherv( &tId,
                        1,
                        MPI_UNSIGNED_LONG,
                        (tReciveBuffer.data()).data(),
                        (tReciveBufferNumber.data()).data(),
                        (tReciveBufferOFFSET.data()).data(),
                        MPI_UNSIGNED_LONG,
                        0,
                        MPI_COMM_WORLD );

                MPI_Gatherv( &tOwner,
                        1,
                        MPI_INT,
                        (tReciveBufferOwner.data()).data(),
                        (tReciveBufferNumber.data()).data(),
                        (tReciveBufferOFFSET.data()).data(),
                        MPI_INT,
                        0,
                        MPI_COMM_WORLD );

                if ( par_rank() == 0 )
                {
                    // Sort this created list
                    std::sort( ( tReciveBuffer.data() ).data(), ( tReciveBuffer.data() ).data() + tReciveBuffer.size() );

                    // use std::unique and std::distance to create list containing all used dof types. This list is unique
                    auto last = std::unique( ( tReciveBuffer.data() ).data(), ( tReciveBuffer.data() ).data() + tReciveBuffer.size() );
                    auto pos  = std::distance( ( tReciveBuffer.data() ).data(), last );

                    tReciveBuffer.resize( pos );

                    MORIS_ERROR( tReciveBuffer.size() <= 2, " Basis has more than 1 uniqe Id" );

                    if( tReciveBuffer.size() == 2 )
                    {
                        MORIS_ERROR( tReciveBuffer( 1 ) == gNoEntityID, " Basis Id are wrong" );
                        MORIS_ERROR( tReciveBuffer( 0 ) != gNoEntityID, " Basis Id is MORIS_MAX_LONGUINT" );
                    }
                }

                if ( par_rank() == 0 )
                {
                    // Sort this created list
                    std::sort( ( tReciveBufferOwner.data() ).data(), ( tReciveBufferOwner.data() ).data() + tReciveBufferOwner.size() );

                    // use std::unique and std::distance to create list containing all used dof types. This list is unique
                    auto last = std::unique( ( tReciveBufferOwner.data() ).data(), ( tReciveBufferOwner.data() ).data() + tReciveBufferOwner.size() );
                    auto pos  = std::distance( ( tReciveBufferOwner.data() ).data(), last );

                    tReciveBufferOwner.resize( pos );

                    MORIS_ERROR( tReciveBufferOwner.size() <= 2, " Basis has more than one owning processor" );

                    if( tReciveBufferOwner.size() == 2 )
                    {
                        MORIS_ERROR( tReciveBufferOwner( 1 ) == gNoProcID, " Basis owners are wrong" );
                        MORIS_ERROR( tReciveBufferOwner( 0 ) != gNoProcID, " Basis Id is MORIS_SINT_MAX" );
                    }
                }
            }

            MORIS_LOG_ERROR("Basis Id and ownership sanity check passed.");
        }

        // -----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

