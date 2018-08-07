/*
 * cl_HMR_Mesh_Base.cpp
 *
 *  Created on: Jun 12, 2018
 *      Author: messe
 */

#include "cl_HMR_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Mesh_Base::Mesh_Base (
                const Parameters       * aParameters,
                Background_Mesh_Base * aBackgroundMesh,
                const uint           & aOrder ) :
                        mParameters( aParameters ),
                        mBackgroundMesh( aBackgroundMesh ),
                        mOrder( aOrder ),
                        mNumberOfDimensions( mParameters->get_number_of_dimensions() ),
                        mNumberOfBasisPerElement(
                                pow( aOrder+1, mParameters->get_number_of_dimensions() ) ),
                        mNumberOfNeighborsPerElement(
                            pow( 3, mParameters->get_number_of_dimensions() ) - 1 )

    {

    }
//------------------------------------------------------------------------------
// protected:
// -----------------------------------------------------------------------------

        void
        Mesh_Base::delete_pointers()
        {

            if ( mAllBasisOnProc.size() > 0 )
            {
                // delete all nodes
                for ( auto tBasis:  mAllBasisOnProc )
                {
                    delete tBasis;
                }

                // clear cell
                mAllBasisOnProc.clear();
            }

            if ( mAllElementsOnProc.size() > 0 )
            {
                // delete all elements
                for ( auto tElement:  mAllElementsOnProc )
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

        void
        Mesh_Base::create_elements()
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
            for ( luint k=0; k<mNumberOfAllElementsOnProc; ++k )
            {
                mAllElementsOnProc( k )
                    = this->create_element( tAllBackgroundElements( k ) );
            }

            this->collect_coarsest_elements();
        }

// -----------------------------------------------------------------------------

        /**
         * collects all coarsest elements on proc including aura
         *
         * @return void
         */
        void
        Mesh_Base::collect_coarsest_elements()
        {
            // count number of coarsest elements
            luint tNumberOfElements = mBackgroundMesh
                    ->get_number_of_coarsest_elements_on_proc_including_aura();

            // reset cell
            mAllCoarsestElementsOnProc.clear();

            // assign memory for cell
            mAllCoarsestElementsOnProc.resize( tNumberOfElements, nullptr );

            // loop over all coarsest elements
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to background mesh element
                Background_Element_Base* tBackElement
                = mBackgroundMesh->get_coarsest_element_by_subdomain_id( e );

                // copy pointer into cell of coarsest elements
                mAllCoarsestElementsOnProc( e ) =
                        mAllElementsOnProc( tBackElement->get_memory_index() );

            }
        }

// -----------------------------------------------------------------------------

        /**
         * unsets the flag of all nodes on proc
         */
        void
        Mesh_Base::unflag_all_basis()
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                tBasis->unflag();
            }
        }

// -----------------------------------------------------------------------------

        /**
         * sets the flag of all nodes on proc
         */
        void
        Mesh_Base::flag_all_basis()
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
        void
        Mesh_Base::determine_elements_connected_to_basis()
        {
            // ask background mesh about active elements in proc
            luint tNumberOfElements = mBackgroundMesh
                    ->get_number_of_active_elements_on_proc_including_aura();

            // loop over all background elements
            for ( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to background element
                Background_Element_Base* tBackElement
                = mBackgroundMesh->get_element_from_proc_domain_including_aura( e );

                // get pointer to Lagrange element
                Element* tElement =
                        mAllElementsOnProc( tBackElement->get_memory_index() );
                // loop over all nodes connected to element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k  )
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
            for ( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to background element
                Background_Element_Base* tBackElement
                = mBackgroundMesh->get_element_from_proc_domain_including_aura( e );

                // get pointer to Lagrange element
                Element* tElement =
                        mAllElementsOnProc( tBackElement->get_memory_index() );
                // loop over all nodes connected to element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k  )
                {
                    // get node pointer
                    Basis* tNode = tElement->get_basis( k );

                    // link this element with the node
                    tNode->insert_element( tElement );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mesh_Base::guess_basis_ownership()
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
                    for( uint k=0; k<tNumberOfElements; ++k )
                    {
                        // get pointer to element
                        Element* tElement = tBasis->get_element( k );

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
                uint tMyRank = par_rank();
                for ( auto tBasis : mAllBasisOnProc )
                {
                    tBasis->set_owner( tMyRank );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Mesh_Base::confirm_basis_ownership()
        {
            // get number of ranks
            uint tNumberOfProcs = par_size();

            // only needed if in parallel mode
            if ( tNumberOfProcs > 1 )
            {
                // make sure that all basis are unflagged
                this->unflag_all_basis();

                // get number of proc neighbors
                uint tNumberOfProcNeighbors
                    = mBackgroundMesh->get_number_of_proc_neighbors();

                // get proc neighbors from background mesh
                auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

                // create cell of matrices to send
                Mat< luint > tEmptyLuint;
                Mat<  uint > tEmptyUint;

                Cell< Mat< luint > > tSendAncestor( tNumberOfProcNeighbors, tEmptyLuint );
                Cell< Mat<  uint > > tSendPedigree( tNumberOfProcNeighbors, tEmptyUint );
                Cell< Mat<  uint > > tSendBasisIndex(  tNumberOfProcNeighbors, tEmptyUint );

                // get my rank
                uint tMyRank = par_rank();

                // loop over all proc neighbors
                for ( uint p = 0; p<tNumberOfProcNeighbors; ++p )
                {
                    // get rank of neighbor
                    uint tNeihgborRank = tProcNeighbors( p );

                    if ( tNeihgborRank != tMyRank && tNeihgborRank != gNoProcNeighbor )
                    {
                        // cell containing basis pointers
                        Cell< Basis* > tBasisInAura;

                        // collect basis within inverse aura
                        this->collect_basis_from_aura( p, 1, tBasisInAura );

                        // calculate addresses of basis to ask for
                        this->encode_foreign_basis_path(
                                tBasisInAura,
                                tNeihgborRank,
                                tSendAncestor( p ),
                                tSendPedigree( p ),
                                tSendBasisIndex( p ) );
                    }
                }

                // send basis indices
                Cell< Mat< uint > > tReceiveBasisIndex;

                // communicate basis indices
                communicate_mats(
                        tProcNeighbors,
                        tSendBasisIndex,
                        tReceiveBasisIndex );

                // free memory
                tSendBasisIndex.clear();

                // ancestors to receive
                Cell< Mat< luint > > tReceiveAncestor;

                // communicate ancestor list
                communicate_mats(
                        tProcNeighbors,
                        tSendAncestor,
                        tReceiveAncestor );

                // free memory
                tSendAncestor.clear();

                // communicate pedigree list
                Cell< Mat<  uint > > tReceivePedigree;
                communicate_mats(
                        tProcNeighbors,
                        tSendPedigree,
                        tReceivePedigree );

                // free memory
                tSendPedigree.clear();

                // matrix with owners to send
                Cell< Mat<  uint > > tSendOwner( tNumberOfProcNeighbors, tEmptyUint );

                // loop over all proc neighbors
                for ( uint p = 0; p<tNumberOfProcNeighbors; ++p )
                {
                    // get number of basis requested by neighbor
                    luint tNumberOfBasis = tReceiveBasisIndex( p ).length();

                    // initialize matrix
                    tSendOwner( p ).set_size( tNumberOfBasis, 1 );

                    // initialize memory conter
                    luint tMemoryCounter = 0;

                    // loop over all basis
                    for( luint k=0; k<tNumberOfBasis; ++k )
                    {
                        // pick requested element
                        Background_Element_Base*
                        tElement = mBackgroundMesh->decode_pedigree_path(
                                tReceiveAncestor( p )( k ),
                                tReceivePedigree( p ),
                                tMemoryCounter );


                        // pick requested basis
                        Basis* tBasis
                            = mAllElementsOnProc( tElement->get_memory_index() )
                                ->get_basis( tReceiveBasisIndex( p )( k ) );

                        // write basis owner into send array
                        tSendOwner( p )( k ) = tBasis->get_owner();
                    }
                }

                // free memory
                tReceiveBasisIndex.clear();
                tReceiveAncestor.clear();
                tReceivePedigree.clear();

                // communicate owners
                Cell< Mat<  uint > > tReceiveOwner;
                communicate_mats(
                        tProcNeighbors,
                        tSendOwner,
                        tReceiveOwner );

                // free memory
                tSendOwner.clear();

                // loop over all proc neighbors
                for ( uint p = 0; p<tNumberOfProcNeighbors; ++p )
                {
                    // get rank of neighbor
                    uint tNeihgborRank = tProcNeighbors( p );

                    if ( tNeihgborRank != tMyRank && tNeihgborRank != gNoProcNeighbor )
                    {
                        // cell containing basis pointers
                        Cell< Basis* > tBasisInAura;

                        // collect basis within inverse aura
                        this->collect_basis_from_aura( p, 1, tBasisInAura );

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

        void
        Mesh_Base::collect_basis_from_aura(
                const uint                & aProcNeighbor,
                const bool                & aUseInverseAura,
                Cell< Basis* >            & aBasisList )
        {

            // clear basis list
            aBasisList.clear();

            // only do something in parallel mode
            if ( par_size() > 1 )
            {
                // cell of elements on background mesh
                Cell< Background_Element_Base* > tBackElements;

                // get element list from background mesh
                mBackgroundMesh->collect_active_elements_from_aura(
                        aProcNeighbor,
                        aUseInverseAura,
                        tBackElements);

                // initialize basis counter
                luint tBasisCount = 0;

                // loop over all background elements
                for ( auto tBackElement : tBackElements )
                {
                    // get pointer to Lagrange element
                    Element* tElement =
                            mAllElementsOnProc( tBackElement->get_memory_index() );

                    // loop over all basis
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

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

                // reserve memory for output cell
                aBasisList.resize( tBasisCount, nullptr );

                // reset basis counter
                tBasisCount = 0;

                // loop over all background elements
                for ( auto tBackElement : tBackElements )
                {
                    // get pointer to Lagrange element
                    Element* tElement =
                            mAllElementsOnProc(
                                    tBackElement->get_memory_index() );

                    // loop over all basis
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

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

// -----------------------------------------------------------------------------

        void
        Mesh_Base::encode_foreign_basis_path(
                           Cell< Basis* > & aBasis,
                           const uint     & aOwner,
                           Mat< luint >   & aElementAncestors,
                           Mat<  uint >   & aElementPedigree,
                           Mat<  uint >   & aElementLocalIndex )
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
                    this->get_reference_element_of_basis(
                            tBasis,
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
                tMemoryCount += mAllElementsOnProc( aElementAncestors( k ) )
                                ->get_background_element()->
                                 get_length_of_pedigree_path();

            }

            // allocate pedigree path
            aElementPedigree.set_size( tMemoryCount, 1 );

            // reset counter
            tMemoryCount = 0;

            // encode pedigree path
            for( luint k=0; k<tCount; ++k )
            {
                // get pointer to element
                Background_Element_Base* tElement
                    =  mAllElementsOnProc( aElementAncestors( k )  )
                        ->get_background_element();
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
                Basis                   * aBasis,
                luint                   & aElementMemoryIndex,
                uint                    & aElementLocalBasisIndex )

        {
            // copy owner of element into temporary variable
            uint tOwner = aBasis->get_owner();

            // find out how many elements are connected to this basis
            uint tNumberOfElements = aBasis->get_element_counter();

            // pointer to background element
            Element* tElement = nullptr;

            // loop over all elements of this basis
            for( uint e=0; e<tNumberOfElements; ++ e )
            {
                // get pointer to element of this basis
                tElement = aBasis->get_element( e );

                // test if element is active and has same owner
                if (    tElement->is_active() && tElement->get_owner() == tOwner )
                {
                    // exit the loop
                    break;
                }
            }

            // get my memory index
            uint tID = aBasis->get_domain_id();

            // set basis to invalid value
            aElementLocalBasisIndex = mNumberOfBasisPerElement;

            // find out local node id on element
            for( uint k=0; k<mNumberOfBasisPerElement; ++k )
            {
                // get pointer to basis
                Basis* tBasis = tElement->get_basis( k );

                // test if basis exists
                if( tBasis != NULL )
                {
                    if ( tBasis->get_domain_id() == tID )
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
            aElementMemoryIndex
                = tElement->get_background_element()->get_memory_index();
        }

// -----------------------------------------------------------------------------

        void
        Mesh_Base::get_basis_coords_of_element(
                      Mat<real>   & aBasisCoords,
                      const luint & aElementIndex,
                      const bool    aSerendipityFlag = false )
        {

            // set number of basis per element
            uint tNumberOfBasisPerElement;

            // test if serendipity flag is set
            if ( aSerendipityFlag )
            {
               // copy only serendipity basis
                tNumberOfBasisPerElement =
                  this->get_number_of_serendipity_basis();
            }
            else
            {
                // copy all basis
                tNumberOfBasisPerElement = mNumberOfBasisPerElement;
            }

            // set size of output matrix
            aBasisCoords.set_size( tNumberOfBasisPerElement,
                                   mNumberOfDimensions );

            // get pointer to element
            Element* tElement = this->get_element( aElementIndex );

            // loop over all nodes of this element
            for( uint k=0; k<tNumberOfBasisPerElement; ++k )
            {
                // get pointer to basis
                Basis* tBasis = tElement->get_basis( k );

                // get array of node coordinates
                const real * tXYZ = tBasis->get_xyz();

                // write coordinates into output matrix
                for( uint i=0; i<mNumberOfDimensions; ++i )
                {
                    aBasisCoords( k, i ) = tXYZ[ i ];
                }
            }
        }

// -----------------------------------------------------------------------------

        uint
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
                /*
                 * the following two would work in theory,
                 * however, the interpolation functions do not exist in MORIS
                 *
                 * case( 16 ) : // QUAD16 -> QUAD12
                {
                    return 12;
                    break;
                }
                case( 64 ) : // HEX64 -> HEX32
                {
                    return 32;
                    break;
                } */
                default :
                {

                    MORIS_ERROR( false, "unknown number of nodes for serendipity");
                    return 0;
                    break;
                }
            }

        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
