/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Mesh_Base.cpp
 *
 */

#include "cl_HMR_Background_Mesh_Base.hpp"    //HMR/src

#include <fstream>

#include "cl_HMR_Background_Facet.hpp"
#include "HMR_Globals.hpp"     //HMR/src
#include "HMR_Tools.hpp"       //HMR/src
#include "cl_Stopwatch.hpp"    //CHR/src
#include "cl_Vector.hpp"         //CNT/src
#include "cl_Bitset.hpp"       //CNT/src
#include "cl_Tracer.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace moris::hmr
{

    //-------------------------------------------------------------------------------

    Background_Mesh_Base::Background_Mesh_Base( const Parameters* aParameters )
            : mParameters( aParameters )
            , mNumberOfDimensions( aParameters->get_number_of_dimensions() )
            , mMaxPolynomial( aParameters->get_max_polynomial() )
            , mPaddingRefinement( ceil( 0.5 * (real)aParameters->get_max_polynomial() ) )
            , mPaddingSize( aParameters->get_padding_size() )
            , mBufferSize( aParameters->get_staircase_buffer() )
            , mNumberOfChildrenPerElement( pow( 2, aParameters->get_number_of_dimensions() ) )
            , mMyRank( par_rank() )
    {
        // make sure that settings are OK
        aParameters->check_sanity();

        // initialize size of Aura Cells ( note this relates to processors not elements in aura )
        uint tSize = pow( 3, aParameters->get_number_of_dimensions() );

        // create empty matrix to initialize fixed size cell
        Matrix< DDLUMat > tEmpty;

        // resize aura cell and use empty matrix as default entry
        mCoarsestAura.resize( tSize, tEmpty );

        // resize inverse aura cell and use empty matrix as default entry
        mCoarsestInverseAura.resize( tSize, tEmpty );

        // reset elements per level table
        for ( uint l = 0; l < gMaxNumberOfLevels; ++l )
        {
            mNumberOfElementsPerLevel[ l ] = 0;
        }

        // set number of neighbors per element
        mNumberOfNeighborsPerElement = std::pow( 3, mNumberOfDimensions ) - 1;
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::print_active_elements()
    {
        luint tNumberOfActiveElements = mActiveElements.size();
        std::fprintf( stdout, "  active elements: %lu\n", (long unsigned int)tNumberOfActiveElements );

        for ( luint k = 0; k < tNumberOfActiveElements; ++k )
        {
            Background_Element_Base* tElement = mActiveElements( k );

            uint tLevel = tElement->get_level();

            if ( tLevel != 0 )
            {
                std::fprintf( stdout, "  el: %4lu   id: %4lu l: %u   o: %u \n", k, (long unsigned int)tElement->get_hmr_id(), (unsigned int)tLevel, (unsigned int)tElement->get_owner() );
            }
            else
            {
                std::fprintf( stdout, "  el: %4lu   id: %4lu l: 0   o: %u \n", k, (long unsigned int)tElement->get_hmr_id(), (unsigned int)tElement->get_owner() );
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::print_active_elements_including_aura()
    {

        this->collect_active_elements_including_aura();

        luint tNumberOfActiveElements = mActiveElementsIncludingAura.size();
        std::fprintf( stdout, "  active elements including aura: %lu\n", (long unsigned int)tNumberOfActiveElements );

        for ( luint k = 0; k < tNumberOfActiveElements; ++k )
        {
            Background_Element_Base*
                    tElement = mActiveElementsIncludingAura( k );

            std::fprintf( stdout, "  el: %lu id: %lu a: %d  r: %d  o: %u\n", k, (long unsigned int)tElement->get_hmr_id(), (int)tElement->is_active( mActivePattern ), (int)tElement->is_refined( mActivePattern ), (unsigned int)tElement->get_owner() );
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::synchronize_coarsest_aura()
    {
        // only synchronize when in parallel mode
        if ( par_size() > 1 )
        {
            uint tNumberOfNeighbors = mMyProcNeighbors.length();

            // matrix to send
            Matrix< IdMat >         tEmpty;
            Vector< Matrix< IdMat > > tSendOwners( tNumberOfNeighbors, tEmpty );

            // loop over all neighbor procs
            for ( uint k = 0; k < tNumberOfNeighbors; ++k )
            {
                // get number of aura elements
                uint tNumberOfAuraElements = mCoarsestInverseAura( k ).length();

                // test if inverse aura exists
                if ( tNumberOfAuraElements > 0 )
                {
                    // insert state switch
                    tSendOwners( k ).set_size( tNumberOfAuraElements, 1 );

                    // loop over all aura elements
                    for ( uint e = 0; e < tNumberOfAuraElements; ++e )
                    {
                        // copy owner of element into state
                        tSendOwners( k )( e ) = mCoarsestElementsIncludingAura( mCoarsestInverseAura( k )( e ) )->get_owner();
                    }
                }
            }

            // matrix to receive
            Vector< Matrix< IdMat > > tReceiveOwners;

            // communicate ownership to neighbors
            communicate_mats( mMyProcNeighbors,
                    tSendOwners,
                    tReceiveOwners );

            // now test which elements are padding and active

            // loop over all neighbors
            for ( uint k = 0; k < tNumberOfNeighbors; ++k )
            {
                uint tNumberOfAuraElements = tReceiveOwners( k ).length();

                if ( tNumberOfAuraElements > 0 )
                {
                    // loop over all aura elements
                    for ( uint e = 0; e < tNumberOfAuraElements; ++e )
                    {
                        // test if owner is MORIS_UINT_MAX. Then this is a padding element
                        if ( tReceiveOwners( k )( e ) == gNoProcOwner )
                        {
                            mCoarsestElementsIncludingAura( mCoarsestAura( k )( e ) )->set_padding_flag();
                        }
                        else
                        {
                            mCoarsestElementsIncludingAura( mCoarsestAura( k )( e ) )->set_active_flag( mActivePattern );
                        }
                    }
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::get_active_elements_on_proc( Matrix< DDLUMat >& aElementIDs )
    {
        // update element table
        this->collect_active_elements();

        // get number of active elements
        luint tNumberOfElements = mActiveElements.size();

        // assign memory
        aElementIDs.set_size( tNumberOfElements, 1 );

        // loop over all active elements
        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // get ID of element
            aElementIDs( k ) = mActiveElements( k )->get_hmr_id();
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::get_active_elements_on_proc_including_aura(
            Matrix< DDLUMat >& aElementIDs )
    {
        // update element table
        this->collect_active_elements_including_aura();

        // get number of active elements
        luint tNumberOfElements = mActiveElementsIncludingAura.size();

        // assign memory
        aElementIDs.set_size( tNumberOfElements, 1 );

        // loop over all active elements
        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // get ID of element
            aElementIDs( k ) = mActiveElementsIncludingAura( k )->get_hmr_id();
        }
    }

    //-------------------------------------------------------------------------------

    luint
    Background_Mesh_Base::count_active_elements() const
    {
        // get number of elements in proc domain
        luint tNumberOfElementsInFrame = mCoarsestElements.size();

        // initialize counter
        luint aCount = 0;

        // loop over frame and count active descendants
        for ( luint k = 0; k < tNumberOfElementsInFrame; ++k )
        {
            mCoarsestElements( k )->get_number_of_active_descendants( mActivePattern, aCount );
        }

        return aCount;
    }

    //-------------------------------------------------------------------------------

    luint
    Background_Mesh_Base::count_active_elements_including_aura() const
    {
        // get number of coarsest elements on proc
        luint tNumberOfCoarsestElements = mCoarsestElementsIncludingAura.size();

        // initialize counter
        luint aCount = 0;

        // loop over frame and count active descendants
        for ( luint iElem = 0; iElem < tNumberOfCoarsestElements; iElem++ )
        {
            mCoarsestElementsIncludingAura( iElem )->get_number_of_active_descendants( mActivePattern, aCount );
        }

        return aCount;
    }

    //-------------------------------------------------------------------------------

    luint
    Background_Mesh_Base::count_all_elements_including_aura() const
    {
        // get number of elements in proc domain
        luint tNumberOfCoarsestElements = mCoarsestElementsIncludingAura.size();

        // initialize counter
        luint aCount = 0;

        // loop over frame and count active descendants
        for ( luint k = 0; k < tNumberOfCoarsestElements; ++k )
        {
            mCoarsestElementsIncludingAura( k )->get_number_of_descendants( aCount );
        }

        return aCount;
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_active_elements()
    {
        // get number of active elements
        luint tCount = this->count_active_elements();

        // clear active element list
        mActiveElements.clear();

        // assign memory
        mActiveElements.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // get number of elements in proc domain
        // frame = domain + aura + padding visible to proc
        luint tNumberOfElementsInFrame = mCoarsestElements.size();

        // loop over all elements in proc domain
        for ( luint iElem = 0; iElem < tNumberOfElementsInFrame; ++iElem )
        {
            // add children to array  mActiveElements
            mCoarsestElements( iElem )->collect_active_descendants( mActivePattern, mActiveElements, tCount );
        }

        // print to console that this operation has been performed
        // MORIS_LOG_INFO( "Updated list of active elements on HMR background mesh." );
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_active_elements_including_aura()
    {
        // get number of active elements
        luint tCount = this->count_active_elements_including_aura();

        // clear active element list
        mActiveElementsIncludingAura.clear();

        // assign memory
        mActiveElementsIncludingAura.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // get number of elements in proc domain (including aura/padding)
        luint tNumberOfElementsInFrame = mCoarsestElementsIncludingAura.size();

        for ( luint iElem = 0; iElem < tNumberOfElementsInFrame; ++iElem )
        {
            mCoarsestElementsIncludingAura( iElem )->collect_active_descendants( mActivePattern, mActiveElementsIncludingAura, tCount );
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::synchronize_refinement_queue()
    {
        // only do something in parallel mode
        if ( par_size() > 1 )
        {
            uint tNumberOfNeighbors = mMyProcNeighbors.length();

            // initialize matrices for sending

            // create empty matrix n
            Matrix< DDLUMat >         tEmptyLuint;
            Vector< Matrix< DDLUMat > > tAncestorListSend;
            tAncestorListSend.resize( tNumberOfNeighbors, { tEmptyLuint } );

            Matrix< DDUMat >         tEmptyUint;
            Vector< Matrix< DDUMat > > tPedigreeListSend;
            tPedigreeListSend.resize( tNumberOfNeighbors, { tEmptyUint } );

            // loop over all procs
            for ( uint p = 0; p < tNumberOfNeighbors; ++p )
            {
                // only do this if there is a neighbor
                if ( mMyProcNeighbors( p ) != gNoProcNeighbor && mMyProcNeighbors( p ) != par_rank() )
                {
                    // initialize element counter
                    luint tNumberOfActiveElements = 0;

                    // get number of elements in inverse aura
                    luint tNumberOfCoarsestElementsOnInverseAura = mCoarsestInverseAura( p ).length();

                    // count active elements on inverse aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnInverseAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( p )( e ) )
                                ->get_number_of_active_descendants( mActivePattern, tNumberOfActiveElements );
                    }

                    // get number of elements on aura
                    luint tNumberOfCoarsestElementsOnAura = mCoarsestAura( p ).length();

                    // count active elements on aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestAura( p )( e ) )
                                ->get_number_of_active_descendants( mActivePattern, tNumberOfActiveElements );
                    }

                    // assign memory
                    Vector< Background_Element_Base* > tActiveElements( tNumberOfActiveElements, nullptr );

                    // reset counter
                    luint tCount = 0;

                    // get active elements on inverse aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnInverseAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( p )( e ) )
                                ->collect_active_descendants( mActivePattern, tActiveElements, tCount );
                    }

                    // get active elements on aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestAura( p )( e ) )
                                ->collect_active_descendants( mActivePattern, tActiveElements, tCount );
                    }

                    // initialize counter for elements
                    luint tElementCounter = 0;

                    // initialize counter for memory needed for pedigree tree
                    luint tMemoryCounter = 0;

                    // create matrix for communication
                    for ( luint k = 0; k < tNumberOfActiveElements; ++k )
                    {
                        // test if element is queued
                        if ( tActiveElements( k )->is_queued_for_refinement() )
                        {
                            // increment counter
                            ++tElementCounter;

                            // get memory needed for pedigree path
                            tMemoryCounter += tActiveElements( k )->get_length_of_pedigree_path();
                        }
                    }

                    if ( tElementCounter > 0 )
                    {
                        // prepare matrix containing ancestors
                        tAncestorListSend( p ).set_size( tElementCounter, 1 );

                        // prepare matrix containing pedigree list
                        tPedigreeListSend( p ).set_size( tMemoryCounter, 1 );

                        // assign matrix for pedigree path

                        // reset counter for elements
                        tElementCounter = 0;

                        // reset pedigree memory counter
                        tMemoryCounter = 0;

                        // create matrices for communication
                        for ( auto tElement : tActiveElements )
                        {
                            if ( tElement->is_queued_for_refinement() )
                            {
                                tElement->encode_pedigree_path( tAncestorListSend( p )( tElementCounter++ ),
                                        tPedigreeListSend( p ),
                                        tMemoryCounter );
                            }
                        }
                    }
                }
            } /* end loop over all procs */

            // initialize matrices for receiving
            Vector< Matrix< DDLUMat > > tAncestorListReceive;
            Vector< Matrix< DDUMat > >  tPedigreeListReceive;

            // communicate ancestor IDs
            communicate_mats( mMyProcNeighbors,
                    tAncestorListSend,
                    tAncestorListReceive );

            // communicate pedigree list
            communicate_mats( mMyProcNeighbors,
                    tPedigreeListSend,
                    tPedigreeListReceive );

            // loop over all received lists
            for ( uint p = 0; p < tNumberOfNeighbors; ++p )
            {
                // get number of elements on refinement list
                luint tNumberOfElements = tAncestorListReceive( p ).length();

                // reset memory counter
                luint tMemoryCounter = 0;

                // loop over all received elements
                for ( uint k = 0; k < tNumberOfElements; ++k )
                {
                    // decode path and get pointer to element
                    Background_Element_Base* tElement = this->decode_pedigree_path( tAncestorListReceive( p )( k ),
                            tPedigreeListReceive( p ),
                            tMemoryCounter );

                    // flag this element for refinement
                    tElement->put_on_refinement_queue();
                }
            }
        }
    }

    //-------------------------------------------------------------------------------

    bool
    Background_Mesh_Base::collect_refinement_queue()
    {
        // synchronize with other procs ( this must be called twice )
        this->synchronize_refinement_queue();
        this->synchronize_refinement_queue();

        // update element active cell
        this->collect_active_elements_including_aura();

        // step 1:  count flagged elements from active list
        // count active elements on proc including aura
        uint tNumActiveElementsOnProc = mActiveElementsIncludingAura.size();

        // reset counter
        luint tCount = 0;

        for ( luint iActiveElem = 0; iActiveElem < tNumActiveElementsOnProc; ++iActiveElem )
        {
            // check if element is flagged
            if ( mActiveElementsIncludingAura( iActiveElem )->is_queued_for_refinement() )
            {
                // increment counter
                ++tCount;

                // perform padding test
                this->check_queued_element_for_padding( mActiveElementsIncludingAura( iActiveElem ) );
            }
        }

        // step 2:  count flagged elements from padding list
        // get number of padding elements
        luint tNumberOfCoarsestPaddingElements = mCoarsestPaddingElements.size();

        // initialize padding counter
        luint tPaddingCount = 0;

        // loop over all coarsest padding elements
        for ( luint iPaddingElem = 0; iPaddingElem < tNumberOfCoarsestPaddingElements; ++iPaddingElem )
        {
            // count descendants
            mCoarsestPaddingElements( iPaddingElem )->get_number_of_descendants( tPaddingCount );
        }

        // array for padding elements
        Vector< Background_Element_Base* > tAllPaddingElements( tPaddingCount, nullptr );

        // reset counter
        tPaddingCount = 0;

        // collect array
        for ( luint iPaddingElem = 0; iPaddingElem < tNumberOfCoarsestPaddingElements; ++iPaddingElem )
        {
            // count descendants
            mCoarsestPaddingElements( iPaddingElem )->collect_descendants( tAllPaddingElements, tPaddingCount );
        }

        // loop over all padding elements
        for ( luint iPaddingElem = 0; iPaddingElem < tPaddingCount; ++iPaddingElem )
        {
            // test if element is flagged for refinement
            if ( tAllPaddingElements( iPaddingElem )->is_queued_for_refinement() )
            {
                // increment counter
                ++tCount;
            }
        }

        // assign memory for queue
        mRefinementQueue.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // step 3: add flagged elements from active list
        for ( luint iActiveElem = 0; iActiveElem < tNumActiveElementsOnProc; ++iActiveElem )
        {
            if ( mActiveElementsIncludingAura( iActiveElem )->is_queued_for_refinement() )
            {
                mRefinementQueue( tCount++ ) = mActiveElementsIncludingAura( iActiveElem );
            }
        }

        // step 4: add flagged elements from padding list

        // loop over all padding elements
        for ( luint iPaddingElem = 0; iPaddingElem < tPaddingCount; ++iPaddingElem )
        {
            // test if element is flagged for refinement
            if ( tAllPaddingElements( iPaddingElem )->is_queued_for_refinement() )
            {
                // copy pointer to queue
                mRefinementQueue( tCount++ ) = tAllPaddingElements( iPaddingElem );
            }
        }

        if ( par_size() == 1 )
        {
            // serial mode, only this counter is important
            return tCount > 0;
        }
        else
        {
            // check number of flagged elements per proc
            Matrix< DDLUMat > tCountersOfAllProcs;
            comm_gather_and_broadcast( tCount, tCountersOfAllProcs );

            // return true if at least one proc wants to refine
            return tCountersOfAllProcs.max() > 0;
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::apply_refinement_queue_to_pattern( const uint aPattern )
    {
        // synchronize with other procs ( this must be called twice )
        this->synchronize_refinement_queue();
        this->synchronize_refinement_queue();

        // update element active cell
        this->collect_active_elements_including_aura();

        // step 1:  count flagged elements from active list
        // count active elements on proc including aura
        uint tActiveElementsOnProc = mActiveElementsIncludingAura.size();

        for ( luint k = 0; k < tActiveElementsOnProc; ++k )
        {
            Background_Element_Base* tBackGroundElement = mActiveElementsIncludingAura( k );

            // check if element is flagged
            if ( tBackGroundElement->is_queued_for_refinement() )
            {
                tBackGroundElement->check_refinement_queue_for_pattern( aPattern );
            }
        }
    }

    //-------------------------------------------------------------------------------

    bool
    Background_Mesh_Base::perform_refinement( const uint aPattern )
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "Background Mesh", "Perform queued refinement on pattern #" + std::to_string( aPattern ) );

        // update buffer size
        mBufferSize = mParameters->get_staircase_buffer();

        // get number of procs
        uint tNumberOfProcs = par_size();

        // update number of elements on queue
        Matrix< DDUMat > tElementCount( tNumberOfProcs, 1, 0 );
        Matrix< DDUMat > tElementCountOld( tNumberOfProcs, 1 );

        uint tNumberOfElementsToBeRefined = 0;

        bool tPatternChange = false;
        uint tOldPattern    = mActivePattern;

        if ( mActivePattern != aPattern )
        {
            tPatternChange = true;

            this->apply_refinement_queue_to_pattern( aPattern );

            // change activation pattern
            mActivePattern = aPattern;
        }

        bool aFlag = false;

        if ( mBufferSize > 0 )
        {
            // update refinement queue with the elements added in the buffer
            this->collect_refinement_queue();
        }

        while ( true )
        {
            // create buffer,if set
            this->create_staircase_buffer();

            // update refinement queue. Refinement queue changes after creating the staircase buffer. so it has to be called again
            aFlag = aFlag || this->collect_refinement_queue();

            // update number of elements on queue
            tNumberOfElementsToBeRefined = mRefinementQueue.size();

            // shift counters
            tElementCountOld = tElementCount;

            // broadcast element count
            comm_gather_and_broadcast( tNumberOfElementsToBeRefined, tElementCount );

            // make sure that element count is the same
            if ( all_true( tElementCount == tElementCountOld ) )
            {
                // exit the loop
                break;
            }
        }

        // empty list of active elements including aura
        mActiveElementsIncludingAura.clear();

        // empty list of active elements
        mActiveElements.clear();

        // perform refinement
        for ( luint iElemToBeRefined = 0; iElemToBeRefined < tNumberOfElementsToBeRefined; ++iElemToBeRefined )
        {
            this->refine_element( mRefinementQueue( iElemToBeRefined ), false );
        }

        // empty queue
        mRefinementQueue.clear();

        if ( tPatternChange )
        {
            mActivePattern = tOldPattern;
        }

        // update database
        this->update_database();

        // return flag
        return aFlag;
    }

    //-------------------------------------------------------------------------------

    uint
    Background_Mesh_Base::calc_child_index( luint aI )
    {
        return aI % 2;
    }

    //-------------------------------------------------------------------------------

    uint
    Background_Mesh_Base::calc_child_index(
            luint aI,
            luint aJ )
    {
        uint tModI = aI % 2;
        uint tModJ = aJ % 2;

        if ( tModJ == 0 )
        {
            if ( tModI == 0 )
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }
        else
        {
            if ( tModI == 0 )
            {
                return 2;
            }
            else
            {
                return 3;
            }
        }
    }

    //-------------------------------------------------------------------------------

    uint
    Background_Mesh_Base::calc_child_index(
            luint aI,
            luint aJ,
            luint aK )
    {
        uint tModI = aI % 2;
        uint tModJ = aJ % 2;
        uint tModK = aK % 2;

        if ( tModK == 0 )
        {
            if ( tModJ == 0 )
            {
                if ( tModI == 0 )
                {
                    return 0;
                }
                else
                {
                    return 1;
                }
            }
            else
            {
                if ( tModI == 0 )
                {
                    return 2;
                }
                else
                {
                    return 3;
                }
            }
        }
        else
        {
            if ( tModJ == 0 )
            {
                if ( tModI == 0 )
                {
                    return 4;
                }
                else
                {
                    return 5;
                }
            }
            else
            {
                if ( tModI == 0 )
                {
                    return 6;
                }
                else
                {
                    return 7;
                }
            }
        }
    }

    //--------------------------------------------------------------------------------

    luint
    Background_Mesh_Base::count_elements_on_level( uint aLevel )
    {
        // get number of coarsest elements
        luint tNumberOfCoarsestElements = mCoarsestElements.size();

        // reset counter
        luint aCount = 0;

        if ( aLevel == 0 )
        {
            return tNumberOfCoarsestElements;
        }
        else
        {
            // loop over all elements on coarsest level
            for ( uint k = 0; k < tNumberOfCoarsestElements; ++k )
            {
                // call recursive function to count elements
                mCoarsestElements( k )->count_elements_on_level( aLevel,
                        aCount );
            }
        }

        // return counter
        return aCount;
    }

    //--------------------------------------------------------------------------------

    luint
    Background_Mesh_Base::count_elements_on_level_including_aura( uint aLevel )
    {
        // get number of coarsest elements
        luint tNumberOfCoarsestElements = mCoarsestElementsIncludingAura.size();

        // reset counter
        luint aCount = 0;

        if ( aLevel == 0 )
        {
            return tNumberOfCoarsestElements;
        }
        else
        {
            // loop over all elements on coarsest level
            for ( uint k = 0; k < tNumberOfCoarsestElements; ++k )
            {
                // call recursive function to count elements
                mCoarsestElementsIncludingAura( k )->count_elements_on_level( aLevel,
                        aCount );
            }
        }

        // return counter
        return aCount;
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::count_elements()
    {
        // reset lookup table
        mNumberOfElementsPerLevel[ 0 ] = mCoarsestElements.size();
        for ( uint l = 1; l < gMaxNumberOfLevels; ++l )
        {
            mNumberOfElementsPerLevel[ l ] = 0;
        }

        // loop over all levels
        for ( uint l = 1; l < gMaxNumberOfLevels; ++l )
        {
            // count elements on this level
            mNumberOfElementsPerLevel[ l ] = this->count_elements_on_level( l );

            // if the current level has no elements,
            // the we cancel the loop
            if ( mNumberOfElementsPerLevel[ l ] == 0 )
            {
                break;
            }
        }

        // part of the counting is also knowing the max level of
        // the subdomain including the aura. This is done here.

        // reset max level
        mMaxLevel = 0;

        for ( uint l = 1; l < gMaxNumberOfLevels; ++l )
        {
            // update max level
            mMaxLevel = l - 1;

            // get number of elements on this level
            luint tNumberOfElementsPerLevel = this->count_elements_on_level_including_aura( l );

            if ( tNumberOfElementsPerLevel == 0 )
            {
                break;
            }
        }

        // Communicate mMaxLevel in case a processor has no flagged elements
        moris::uint tMaxLevelResult = max_all( mMaxLevel );

        mMaxLevel = tMaxLevelResult;
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_elements_on_level(
            uint                              aLevel,
            Vector< Background_Element_Base* >& aElementList )
    {
        // get number of elements from lookup table
        auto tNumberOfElements = mNumberOfElementsPerLevel[ aLevel ];

        // empty cell
        aElementList.clear();

        // initialize counter
        luint tCount = 0;

        if ( tNumberOfElements > 0 )
        {
            // assign memory for cell
            aElementList.resize( tNumberOfElements, nullptr );
            // get number of coarsest elements
            luint tNumberOfCoarsestElements = mNumberOfElementsPerLevel[ 0 ];

            // loop over all elements on coarsest level
            for ( uint k = 0; k < tNumberOfCoarsestElements; ++k )
            {
                // call recursive function to collect elements
                mCoarsestElements( k )->collect_elements_on_level( aLevel,
                        aElementList,
                        tCount );
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_elements_on_level_including_aura(
            uint                              aLevel,
            Vector< Background_Element_Base* >& aElementList )
    {
        // get number of elements from lookup table
        auto tNumberOfElements = this->count_elements_on_level_including_aura( aLevel );

        // empty cell
        aElementList.clear();

        // initialize counter
        luint tCount = 0;

        if ( tNumberOfElements > 0 )
        {
            // assign memory for cell
            aElementList.resize( tNumberOfElements, nullptr );
            // get number of coarsest elements
            luint tNumberOfCoarsestElements = mCoarsestElementsIncludingAura.size();

            // loop over all elements on coarsest level
            for ( uint k = 0; k < tNumberOfCoarsestElements; ++k )
            {
                // call recursive function to collect elements
                mCoarsestElementsIncludingAura( k )->collect_elements_on_level( aLevel,
                        aElementList,
                        tCount );
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_elements_on_level_within_proc_domain(
            uint                              aLevel,
            Vector< Background_Element_Base* >& aElementList )
    {
        // get number of elements from lookup table
        auto tNumberOfElements = this->count_elements_on_level( aLevel );

        // empty cell
        aElementList.clear();

        // initialize counter
        luint tCount = 0;

        if ( tNumberOfElements > 0 )
        {
            // assign memory for cell
            aElementList.resize( tNumberOfElements, nullptr );
            // get number of coarsest elements
            luint tNumberOfCoarsestElements = mCoarsestElements.size();

            // loop over all elements on coarsest level
            for ( uint k = 0; k < tNumberOfCoarsestElements; ++k )
            {
                // call recursive function to collect elements
                mCoarsestElements( k )->collect_elements_on_level( aLevel,
                        aElementList,
                        tCount );
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_all_elements( Vector< Background_Element_Base* >& aElementList )
    {
        // clear element list
        aElementList.clear();

        // initialize counter
        luint tCount = this->count_all_elements_including_aura();

        // reserve memory on element list
        aElementList.resize( tCount, nullptr );

        // get number of elements on coarsest level including aura and padding elements
        luint tNumberOfCoarsestElements = mCoarsestElementsIncludingAura.size();

        // reset counter
        tCount = 0;

        // loop over all elements on coarsest level
        for ( luint iElem = 0; iElem < tNumberOfCoarsestElements; ++iElem )
        {
            mCoarsestElementsIncludingAura( iElem )->collect_descendants( aElementList, tCount );
        }

        // sets memory index for all elements
        for ( luint iActiveElem = 0; iActiveElem < tCount; ++iActiveElem )
        {
            aElementList( iActiveElem )->set_memory_index( iActiveElem );
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_neighbors()
    {
        // update counting table
        this->count_elements();

        this->collect_neighbors_on_level_zero();

        // element list on level
        Vector< Background_Element_Base* > tElementList;

        // collect higher level neighbors if any refinement exists
        if ( mMaxLevel > 0 )
        {
            // loop over all levels
            for ( uint l = 1; l <= mMaxLevel; ++l )
            {
                // pick all elements from this level
                this->collect_elements_on_level_including_aura( l, tElementList );

                // get number of elements
                luint tNumberOfElements = tElementList.size();

                // loop over all elements on this level
                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    // collect neighbors
                    tElementList( k )->collect_neighbors( mActivePattern );
                }
            }
        }

        // output that this operation was performed
        // MORIS_LOG_INFO( "Updated neighbors for active or refined elements." );
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::update_element_indices()
    {
        // get number of elements on current proc
        luint tNumberOfElements = mActiveElements.size();

        // FIXME Loop over all patterns

        // update local element indices
        for ( luint iElem = 0; iElem < tNumberOfElements; ++iElem )
        {
            mActiveElements( iElem )->set_domain_index( mActivePattern, iElem );
        }

        // communicate indices to other proc
        this->synchronize_local_indices();

        // communicate number of elements with other procs
        Matrix< DDLUMat > tElementsPerProc;
        comm_gather_and_broadcast( tNumberOfElements, tElementsPerProc );

        // get number of procs
        luint tNumberOfProcs = par_size();

        // create offset per proc
        Matrix< DDLUMat > tProcOffset( tNumberOfProcs, 1 );
        tProcOffset( 0 ) = 0;

        for ( uint iProc = 1; iProc < tNumberOfProcs; ++iProc )
        {
            tProcOffset( iProc ) = tProcOffset( iProc - 1 ) + tElementsPerProc( iProc - 1 );
        }

        mMaxElementDomainIndex = tProcOffset( tNumberOfProcs - 1 ) + tElementsPerProc( tNumberOfProcs - 1 );

        // get number of all elements on proc
        tNumberOfElements = mActiveElementsIncludingAura.size();

        // create global indices
        for ( luint iBgElem = 0; iBgElem < tNumberOfElements; ++iBgElem )
        {
            // get pointer to element
            Background_Element_Base* tElement = mActiveElementsIncludingAura( iBgElem );

            // get owner of element
            auto tOwner = tElement->get_owner();

            // get local index of element
            auto tIndex = tElement->get_hmr_index( mActivePattern );

            // calculate global index
            tElement->set_domain_index( mActivePattern,
                    tIndex + tProcOffset( tOwner ) );
        }

        // print to console that this operation has been performed
        // MORIS_LOG_INFO( "Updated element indices on background mesh." );
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::synchronize_local_indices()
    {
        // only do something in parallel mode
        if ( par_size() > 1 )
        {
            uint tNumberOfNeighbors = mMyProcNeighbors.length();

            // create empty matrix n
            Matrix< DDLUMat >         tEmpty;
            Vector< Matrix< DDLUMat > > tIndexListSend;
            tIndexListSend.resize( tNumberOfNeighbors, { tEmpty } );

            // initialize matrices for receiving
            Vector< Matrix< DDLUMat > > tIndexListReceive;

            // loop over all procs
            for ( uint p = 0; p < tNumberOfNeighbors; ++p )
            {
                // only do this if there is a neighbor
                if ( mMyProcNeighbors( p ) != gNoProcNeighbor && mMyProcNeighbors( p ) != par_rank() )
                {
                    // number of elements on inverse aura
                    luint tNumberOfElementsOnInverseAura = 0;

                    // get number of elements in inverse aura
                    luint tNumberOfCoarsestElementsOnInverseAura = mCoarsestInverseAura( p ).length();

                    // count active elements on inverse aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnInverseAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( p )( e ) )
                                ->get_number_of_active_descendants( mActivePattern, tNumberOfElementsOnInverseAura );
                    }

                    // assign memory for cell of elements
                    Vector< Background_Element_Base* > tElements( tNumberOfElementsOnInverseAura, nullptr );

                    // reset counter
                    luint tCount = 0;

                    // collect elements on inverse aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnInverseAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( p )( e ) )
                                ->collect_active_descendants( mActivePattern, tElements, tCount );
                    }

                    // assign matrix of indices to be sent
                    tIndexListSend( p ).set_size( tCount, 1 );

                    // loop over all elements
                    for ( luint k = 0; k < tCount; ++k )
                    {
                        // copy local index into matrix
                        tIndexListSend( p )( k ) = tElements( k )->get_hmr_index( mActivePattern );
                    }
                }
            } /* end loop over all procs */

            // communicate matrices
            communicate_mats( mMyProcNeighbors,
                    tIndexListSend,
                    tIndexListReceive );

            // loop over all received lists
            for ( uint p = 0; p < tNumberOfNeighbors; ++p )
            {
                if ( mMyProcNeighbors( p ) != gNoProcNeighbor && mMyProcNeighbors( p ) != par_rank() )
                {
                    // get number of coarsest elements on aura
                    luint tNumberOfCoarsestElementsOnAura = mCoarsestAura( p ).length();

                    // get number of active elements on aura
                    luint tNumberOfElementsOnAura = tIndexListReceive( p ).length();

                    // assign memory for cell of elements
                    Vector< Background_Element_Base* > tElements( tNumberOfElementsOnAura, nullptr );

                    // initialize counter
                    luint tCount = 0;

                    // collect active elements on aura
                    for ( luint e = 0; e < tNumberOfCoarsestElementsOnAura; ++e )
                    {
                        mCoarsestElementsIncludingAura( mCoarsestAura( p )( e ) )
                                ->collect_active_descendants( mActivePattern, tElements, tCount );
                    }

                    // make sure that number of elements is correct
                    MORIS_ERROR( tCount == tNumberOfElementsOnAura,
                            "synchronize_local_indices: Number of elements on aura does not match." );

                    // loop over all elements on aura
                    for ( luint k = 0; k < tNumberOfElementsOnAura; ++k )
                    {
                        // write index into element
                        tElements( k )->set_domain_index( mActivePattern, tIndexListReceive( p )( k ) );
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_active_elements_from_aura(
            uint                              aProcNeighbor,
            uint                              aMode,
            Vector< Background_Element_Base* >& aElementList )
    {
        // clear element list
        aElementList.clear();

        // only do something in parallel mode
        if ( par_size() > 1 )
        {
            // only do something of proc neighbor exists and is not myself
            if ( mMyProcNeighbors( aProcNeighbor ) != gNoProcNeighbor && mMyProcNeighbors( aProcNeighbor ) != par_rank() )
            {
                // initialize element counter
                luint tCount = 0;

                if ( aMode == 1 || aMode == 2 )
                {
                    // get number of coarsest elements of inverse aura
                    luint tNumberOfElements = mCoarsestInverseAura( aProcNeighbor ).length();

                    // count elements on coarsest inverse aura
                    for ( luint e = 0; e < tNumberOfElements; ++e )
                    {
                        // count children of this element
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( aProcNeighbor )( e ) )
                                ->get_number_of_active_descendants( mActivePattern, tCount );
                    }
                }
                if ( aMode == 0 || aMode == 2 )
                {
                    // get number of coarsest elements of aura
                    luint tNumberOfElements = mCoarsestAura( aProcNeighbor ).length();

                    // count elements on coarsest inverse aura
                    for ( luint e = 0; e < tNumberOfElements; ++e )
                    {
                        // count children of this element
                        mCoarsestElementsIncludingAura( mCoarsestAura( aProcNeighbor )( e ) )
                                ->get_number_of_active_descendants( mActivePattern, tCount );
                    }
                }

                // assign memory for output cell
                aElementList.resize( tCount, nullptr );

                // reset counter
                tCount = 0;

                if ( aMode == 1 || aMode == 2 )
                {
                    // get number of coarsest elements of inverse aura
                    luint tNumberOfElements = mCoarsestInverseAura( aProcNeighbor ).length();

                    // loop over elements on coarsest inverse aura
                    for ( luint e = 0; e < tNumberOfElements; ++e )
                    {
                        // count children of this element
                        mCoarsestElementsIncludingAura( mCoarsestInverseAura( aProcNeighbor )( e ) )
                                ->collect_active_descendants( mActivePattern, aElementList, tCount );
                    }
                }
                if ( aMode == 0 || aMode == 2 )
                {
                    // loop over elements on coarsest  aura
                    luint tNumberOfElements = mCoarsestAura( aProcNeighbor ).length();

                    // count elements on coarsest inverse aura
                    for ( luint e = 0; e < tNumberOfElements; ++e )
                    {
                        // count children of this element
                        mCoarsestElementsIncludingAura( mCoarsestAura( aProcNeighbor )( e ) )
                                ->collect_active_descendants( mActivePattern, aElementList, tCount );
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::create_staircase_buffer_for_element(
            Background_Element_Base* aElement,
            luint&                   aElementCounter,
            uint                     aHalfBuffer )
    {
        // cell containing neighbors of each parent
        Vector< Background_Element_Base* > tNeighbors;

        // get level of this element
        auto tLevel = aElement->get_level();

        // only do something if level > 0
        if ( tLevel > 0 )
        {
            const luint* tElIJK = aElement->get_ijk();

            Background_Element_Base* tParent = aElement->get_parent();

            // get neighbors
            tParent->get_neighbors_from_same_level( aHalfBuffer, tNeighbors );

            // calculate length of cell
            uint tNumberOfNeighbors = tNeighbors.size();

            // check neighbors of element
            for ( uint iElem = 0; iElem < tNumberOfNeighbors; ++iElem )
            {
                // get neighbor
                Background_Element_Base* tNeighbor = tNeighbors( iElem );

                // test all neighbors
                // test if neighbor is active and was not flagged
                if ( tNeighbor->is_active( mActivePattern ) && !tNeighbor->is_queued_for_refinement() )
                {
                    // test position of neighbor
                    const luint* tNeighborIJK = tNeighbor->get_ijk();

                    // max ijk distance between both elements
                    luint tMaxDistance = 0;

                    // calculate distance to element to refine
                    for ( uint i = 0; i < mNumberOfDimensions; ++i )
                    {
                        luint tNeighborI = 2 * tNeighborIJK[ i ];

                        // the distance in i-j-k direction
                        luint tDistance = 0;

                        // we must be a bit careful with luint
                        // therefore, do a size test
                        if ( tNeighborI >= tElIJK[ i ] )
                        {
                            tDistance = tNeighborI - tElIJK[ i ];
                        }
                        else
                        {
                            tDistance = tElIJK[ i ] - tNeighborI - 1;
                        }
                        if ( tDistance > tMaxDistance )
                        {
                            tMaxDistance = tDistance;
                        }
                    }

                    // flag this element if it is close
                    if ( tMaxDistance <= mBufferSize )
                    {
                        // flag this neighbor
                        tNeighbor->put_on_refinement_queue();

                        // increment element counter
                        ++aElementCounter;

                        // create staircase buffer for neighbor
                        this->create_staircase_buffer_for_element( tNeighbor,
                                aElementCounter,
                                aHalfBuffer );
                    }
                }
            } // end for: each neighboring element
        } // end if: is refined level
    } // end function: Background_Mesh_Base::create_staircase_buffer_for_element()

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::create_staircase_buffer()
    {
        // only do something if mBufferSize is set
        if ( mBufferSize > 0 )
        {
            // // update refinement queue
            // this->collect_refinement_queue();

            // report on this operation
            MORIS_LOG_INFO( "Creating staircase buffer of width %u.", (unsigned int)mBufferSize );

            // element counter
            luint tElementCounter = 0;

            // only need half buffer
            uint tHalfBuffer = ceil( 0.5 * (real)mBufferSize );

            // create staircase buffer
            for ( auto tElement : mRefinementQueue )
            {
                this->create_staircase_buffer_for_element( 
                        tElement,
                        tElementCounter,
                        tHalfBuffer );
            }
        }
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::collect_coarsest_padding_elements()
    {
        // number of all coarsest elements
        luint tNumberOfElements = mCoarsestElementsIncludingAura.size();

        // initialize counter
        luint tCount = 0;

        // count padding elements
        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // test if this is a padding element
            if ( mCoarsestElementsIncludingAura( k )->is_padding() )
            {
                // increment counter
                ++tCount;
            }
        }

        // assign memory
        mCoarsestPaddingElements.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all elements
        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // test if this is a padding element
            if ( mCoarsestElementsIncludingAura( k )->is_padding() )
            {
                // copy element
                mCoarsestPaddingElements( tCount++ ) = mCoarsestElementsIncludingAura( k );
            }
        }
    }

    //--------------------------------------------------------------------------------

    Background_Element_Base*
    Background_Mesh_Base::decode_pedigree_path(
            luint                   aAncestorID,
            const Matrix< DDUMat >& aPedigreeList,
            luint&                  aCounter )
    {
        // pick ancestor element
        auto tIndex = this->calc_subdomain_id_from_global_id( 0, aAncestorID );

        Background_Element_Base* aElement = mCoarsestElementsIncludingAura( tIndex );

        // global bitset
        Bitset< gBitsetSize > tBitset;

        // initialize bit counter
        uint tCount = 0;

        // get level of element
        uint tLevel = aPedigreeList( aCounter++ );

        // get number of dimensions from settings
        uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

        if ( tLevel > 0 )
        {
            // get number of chars to read
            luint tNumberOfChars = ceil( tLevel * ( (real)( tNumberOfDimensions ) ) / 32 );

            // create bitset from char
            for ( luint k = 0; k < tNumberOfChars; ++k )
            {
                Bitset< 32 > tChar( aPedigreeList( aCounter++ ) );

                // copy bits into path bitset
                for ( uint i = 0; i < 32; ++i )
                {
                    if ( tChar.test( i ) )
                    {
                        tBitset.set( tCount );
                    }

                    // increment bit counter
                    ++tCount;
                }
            }

            // reset bit counter
            tCount = 0;

            // loop over all levels
            for ( uint l = 0; l < tLevel; ++l )
            {
                // decode child index from bitset
                Bitset< 3 > tChildIndex;
                for ( uint i = 0; i < tNumberOfDimensions; ++i )
                {
                    if ( tBitset.test( tCount++ ) )
                    {
                        tChildIndex.set( i );
                    }
                }

                aElement = aElement->get_child( tChildIndex.to_ulong() );
            }
        }

        // return element
        return aElement;
    }

    //--------------------------------------------------------------------------------

    void
    Background_Mesh_Base::save_to_vtk( const std::string& aFilePath )
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "Background Mesh", "Save to VTK" );

        // modify filename
        std::string tFilePath = parallelize_path( aFilePath );

        // open the file
        std::ofstream tFile( tFilePath, std::ios::binary );

        // containers
        float tFChar = 0;
        int   tIChar = 0;

        tFile << "# vtk DataFile Version 3.0" << std::endl;
        tFile << "GO BUFFS!" << std::endl;
        tFile << "BINARY" << std::endl;

        // collect all elements on proc
        // total number of elements
        luint tNumberOfAllElements = this->count_all_elements_including_aura();

        // collect all elements on proc
        Vector< Background_Element_Base* > tAllElements( tNumberOfAllElements, nullptr );
        this->collect_all_elements( tAllElements );

        // initialize element counter
        luint tNumberOfElements = 0;

        // count number of elements that don't have children
        for ( auto tElement : tAllElements )
        {
            if ( !tElement->has_children() )
            {
                // increment element counter
                ++tNumberOfElements;
            }
        }

        // number of nodes per element
        uint tNumberOfNodesPerElement = std::pow( 2, mParameters->get_number_of_dimensions() );

        // count number of nodes
        uint tNumberOfNodes = tNumberOfNodesPerElement * tNumberOfElements;

        // write node data
        tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

        tFile << "POINTS " << tNumberOfNodes << " float" << std::endl;

        // temporary matrix containing corder nodes
        Matrix< DDRMat > tNodes( mParameters->get_number_of_dimensions(), tNumberOfNodesPerElement );

        // VTK cell type
        int tCellType = 0;

        if ( mParameters->get_number_of_dimensions() == 2 )
        {
            // loop over all elements
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {

                    // ask background mesh for corner nodes
                    this->calc_corner_nodes_of_element( tElement, tNodes );

                    // write node coordinates to file
                    for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
                    {
                        tFChar = swap_byte_endian( (float)tNodes( 0, k ) );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                        tFChar = swap_byte_endian( (float)tNodes( 1, k ) );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                        tFChar = swap_byte_endian( (float)0 );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                    }
                }
            }

            // VTK cell type for Quad4
            tCellType = 9;
        }
        else if ( mParameters->get_number_of_dimensions() == 3 )
        {
            // loop over all elements
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    // ask background mesh for corner nodes
                    this->calc_corner_nodes_of_element( tElement, tNodes );

                    // write node coordinates to file
                    for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
                    {
                        tFChar = swap_byte_endian( (float)tNodes( 0, k ) );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                        tFChar = swap_byte_endian( (float)tNodes( 1, k ) );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                        tFChar = swap_byte_endian( (float)tNodes( 2, k ) );
                        tFile.write( (char*)&tFChar, sizeof( float ) );
                    }
                }
            }

            // VTK cell type for Hex8
            tCellType = 12;
        }

        // create new line
        tFile << std::endl;

        // can only write element data if vtk map exists
        if ( tCellType != 0 )
        {
            // value to write in VTK file
            int tNumberOfNodesVTK = swap_byte_endian( (int)tNumberOfNodesPerElement );

            // reset node counter
            int tCount = 0;

            // write header for cells
            tFile << "CELLS " << tNumberOfElements << " "
                  << ( tNumberOfNodesPerElement + 1 ) * tNumberOfElements << std::endl;

            // loop over all elements
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    tFile.write( (char*)&tNumberOfNodesVTK, sizeof( int ) );

                    // loop over all nodes of this element
                    for ( uint k = 0; k < tNumberOfNodesPerElement; ++k )
                    {
                        // write node to mesh file
                        tIChar = swap_byte_endian( tCount++ );
                        tFile.write( (char*)&tIChar, sizeof( int ) );
                    }
                }
            }

            // write cell types
            tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
            tIChar = swap_byte_endian( tCellType );
            for ( luint k = 0; k < tNumberOfElements; ++k )
            {
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }

            // write element data
            tFile << "CELL_DATA " << tNumberOfElements << std::endl;

            // write element ID
            tFile << "SCALARS ELEMENT_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    tIChar = swap_byte_endian( (int)tElement->get_hmr_id() );
                    tFile.write( (char*)&tIChar, sizeof( int ) );
                }
            }

            // write proc owner
            tFile << "SCALARS ELEMENT_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    if ( tElement->is_padding() )
                    {
                        tIChar = swap_byte_endian( (int)-1 );
                    }
                    else
                    {
                        tIChar = swap_byte_endian( (int)tElement->get_owner() );
                    }
                    tFile.write( (char*)&tIChar, sizeof( int ) );
                }
            }
            tFile << std::endl;

            // write level
            tFile << "SCALARS ELEMENT_LEVEL int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    tIChar = swap_byte_endian( (int)tElement->get_level() );
                    tFile.write( (char*)&tIChar, sizeof( float ) );
                }
            }
            tFile << std::endl;

            // write memory index
            tFile << "SCALARS ELEMENT_MEMORY_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for ( auto tElement : tAllElements )
            {
                if ( !tElement->has_children() )
                {
                    tIChar = swap_byte_endian( (int)tElement->get_memory_index() );
                    tFile.write( (char*)&tIChar, sizeof( int ) );
                }
            }
            tFile << std::endl;

            // close the output file
            tFile.close();
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::reset_pattern( uint aPattern )
    {
        // make sure the pattern is valid
        MORIS_ERROR( aPattern < gNumberOfPatterns, "Background_Mesh_Base::reset_pattern() - Invalid Pattern." );

        // Cell containing all elements
        Vector< Background_Element_Base* > tElements;

        // collect all elements on this level including the padding/aura around it
        this->collect_elements_on_level_including_aura( 0, tElements );

        for ( auto tElement : tElements )
        {
            // do only of not a padding element
            if ( tElement->is_padding() )
            {
                // padding elements are always refined
                tElement->set_refined_flag( aPattern );
            }
            else
            {
                // activate this element
                tElement->set_active_flag( aPattern );
            }
        }

        // for all higher levels
        for ( uint l = 1; l <= mMaxLevel; ++l )
        {
            this->collect_elements_on_level_including_aura( l, tElements );
            for ( auto tElement : tElements )
            {
                if ( tElement->is_padding() )
                {
                    // padding elements are always refined
                    tElement->set_refined_flag( aPattern );
                    // tElement->deactivate( aPattern );
                }
                else
                {
                    // deactivate this element
                    tElement->deactivate( aPattern );
                }
            }
        }

        // report on this operation
        MORIS_LOG_INFO( "Reset pattern #%i on background mesh.", aPattern );

    } // end function: Background_Mesh_Base::reset_pattern()

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::clone_pattern(
            uint aSource,
            uint aTarget )
    {
        MORIS_ERROR( aSource < gNumberOfPatterns, "Source pattern invalid." );
        MORIS_ERROR( aTarget < gNumberOfPatterns, "Target pattern invalid." );

        tic tTimer;

        // reset target pattern
        this->reset_pattern( aTarget );

        // remember original pattern
        auto tOldPattern = mActivePattern;

        // select target pattern
        mActivePattern = aTarget;

        for ( uint l = 0; l < mMaxLevel; ++l )
        {
            // cell containing elements on current level
            Vector< Background_Element_Base* > tElements;

            // get elements from level that belong to myself
            this->collect_elements_on_level( l, tElements );

            // if this is the first level, activate all elements
            for ( auto tElement : tElements )
            {
                if ( tElement->is_refined( aSource ) )
                {
                    tElement->put_on_refinement_queue();
                }
            }

            // perform refinement
            this->perform_refinement( aTarget );
        }

        // get back to old pattern
        mActivePattern = tOldPattern;

        // stop timer
        real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

        MORIS_LOG_INFO( "%s Cloned refinement pattern.",
                proc_string().c_str() );
        MORIS_LOG_INFO( "Cloning took %5.3f seconds.",
                (double)tElapsedTime / 1000 );
        MORIS_LOG_INFO( " " );
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::unite_patterns(
            uint aSourceA,
            uint aSourceB,
            uint aTarget )
    {
        MORIS_ERROR( aSourceA < gNumberOfPatterns, "Source pattern A invalid." );
        MORIS_ERROR( aSourceB < gNumberOfPatterns, "Source pattern B invalid." );

        MORIS_ERROR( aTarget < gNumberOfPatterns, "Target pattern invalid." );

        // reset target pattern
        this->reset_pattern( aTarget );

        // remember original pattern
        // auto tOldPattern = mActivePattern;

        // mActivePattern = aTarget ;
        this->set_activation_pattern( aTarget );

        for ( uint l = 0; l < mMaxLevel; ++l )
        {
            // cell containing elements on current level
            Vector< Background_Element_Base* > tElements;

            // get elements from level that belong to myself
            this->collect_elements_on_level( l, tElements );

            // if this is the first level, activate all elements
            for ( auto tElement : tElements )
            {
                if ( tElement->is_refined( aSourceA ) || tElement->is_refined( aSourceB ) )
                {
                    // flag this element for refinement
                    tElement->put_on_refinement_queue();
                }
            }

            // perform refinement
            this->perform_refinement( aTarget );
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::unite_patterns(
            const moris::Vector< uint >& aSourcePattern,
            const uint                 aTarget )
    {
        for ( uint Il = 0; Il < aSourcePattern.size(); ++Il )
        {
            MORIS_ERROR( aSourcePattern( Il ) < gNumberOfPatterns, "Source pattern %-5i invalid.", Il );
        }

        MORIS_ERROR( aTarget < gNumberOfPatterns, "Target pattern invalid." );

        // reset target pattern
        this->reset_pattern( aTarget );

        this->set_activation_pattern( aTarget );

        for ( uint Ii = 0; Ii < mMaxLevel; ++Ii )
        {
            // cell containing elements on current level
            Vector< Background_Element_Base* > tElements;

            // get elements from level that belong to myself
            this->collect_elements_on_level( Ii, tElements );

            // if this is the first level, activate all elements
            for ( auto tElement : tElements )
            {
                bool tIsRefined = false;

                // check if one of the pattern is refined for this element
                for ( uint Ia = 0; Ia < aSourcePattern.size(); ++Ia )
                {
                    tIsRefined = tIsRefined || tElement->is_refined( aSourcePattern( Ia ) );
                }

                if ( tIsRefined )
                {
                    // flag this element for refinement
                    tElement->put_on_refinement_queue();
                }
            }

            // perform refinement
            this->perform_refinement( aTarget );
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::copy_pattern(
            uint aSource,
            uint aTarget )
    {
        if ( aSource != aTarget )
        {
            Vector< Background_Element_Base* > tElementList;

            // collect all elements
            this->collect_all_elements( tElementList );

            // loop over all elements
            for ( auto tElement : tElementList )
            {
                // test if element is refined on source pattern
                if ( tElement->is_active( aSource ) )
                {
                    // activate element on output pattern
                    tElement->set_active_flag( aTarget );
                }
                else if ( tElement->is_refined( aSource ) )
                {
                    // refine element on output pattern
                    tElement->set_refined_flag( aTarget );
                }
                else
                {
                    // deactivate element
                    tElement->deactivate( aTarget );
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::create_facets()
    {
        // report on this operation
        MORIS_LOG_INFO( "Creating facets on HMR background mesh" );

        // uint tPattern = mParameters->get_lagrange_output_pattern();
        // this->set_activation_pattern( tPattern );

        // loop over all levels
        for ( uint l = 0; l <= mMaxLevel; ++l )
        {
            // list of all elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level
            this->collect_elements_on_level_including_aura( l, tElementList );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElementList )
            {
                // create faces
                tElement->create_facets();
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::create_faces_and_edges()
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "Background Mesh", "Create faces and edges" );

        // // select output pattern
        // uint tPattern = mParameters->get_lagrange_output_pattern();
        // this->set_activation_pattern( tPattern );

        // loop over all levels
        for ( uint iLevel = 0; iLevel <= mMaxLevel; ++iLevel )
        {
            // list of all elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level
            this->collect_elements_on_level_including_aura( iLevel, tElementList );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElementList )
            {
                // create faces
                tElement->create_facets();

                // create edges
                tElement->create_edges();
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::delete_faces()
    {
        // loop over all levels
        for ( uint l = 0; l <= mMaxLevel; ++l )
        {
            // list of all elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level
            this->collect_elements_on_level_including_aura( l, tElementList );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElementList )
            {
                // create faces
                tElement->delete_facets();
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::delete_edges()
    {
        // loop over all levels
        for ( uint l = 0; l <= mMaxLevel; ++l )
        {
            // list of all elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level
            this->collect_elements_on_level_including_aura( l, tElementList );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElementList )
            {
                // create edges
                tElement->delete_edges();

                tElement->init_edge_container();
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Background_Mesh_Base::reset_neighbors()
    {
        // loop over all levels
        for ( uint l = 0; l <= mMaxLevel; ++l )
        {
            // list of all elements on this level
            Vector< Background_Element_Base* > tElementList;

            // collect all elements on this level
            this->collect_elements_on_level_including_aura( l, tElementList );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElementList )
            {
                // create edges
                tElement->reset_neighbors();
            }
        }
    }

    // -----------------------------------------------------------------------------

    /**
     * updates the database according to selected pattern
     */
    void
    Background_Mesh_Base::update_database()
    {
        // report this operation
        MORIS_LOG_INFO( "Updating HMR Database" );
        //Tracer tTracer( "HMR", "Database", "Update" );

        this->collect_active_elements();
        this->collect_active_elements_including_aura();
        this->update_element_indices();
        this->collect_neighbors();
    }

    // -----------------------------------------------------------------------------

    /**
     * Collect background elements on side set.
     * Side set numbers see Exodus II : A Finite Element Data Model, p. 13
     */
    void
    Background_Mesh_Base::collect_side_set_elements(
            uint                              aPattern,
            uint                              aSideOrdinal,
            Vector< Background_Element_Base* >& aElements )
    {
        luint tElementCounter = 0;

        // step 1: collect elements on coarsest level
        Vector< Background_Element_Base* > tCoarsestElements;
        collect_coarsest_elements_on_side( aSideOrdinal, tCoarsestElements );

        switch ( aSideOrdinal )
        {
            case ( 1 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_1(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_1( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            case ( 2 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_2(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_2( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            case ( 3 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_3(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_3( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            case ( 4 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_4(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_4( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            case ( 5 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_5(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_5( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            case ( 6 ):
            {
                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only count non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->get_number_of_active_descendants_on_side_6(
                                aPattern, tElementCounter );
                    }
                }

                // allocate output matrix
                aElements.resize( tElementCounter, nullptr );

                // reset counter
                tElementCounter = 0;

                // loop over all coarsest elements
                for ( Background_Element_Base* tElement : tCoarsestElements )
                {
                    // only add descendants of non padding elements
                    if ( !tElement->is_padding() )
                    {
                        tElement->collect_active_descendants_on_side_6( aPattern,
                                aElements,
                                tElementCounter );
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Invalid Side set ordinal." );
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::reset_min_refinement_levels()
    {
        // loop over all levels
        for ( uint l = 0; l <= mMaxLevel; ++l )
        {
            // collect elements from this level
            moris::Vector< Background_Element_Base* > tElements;

            this->collect_elements_on_level_including_aura( l, tElements );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElements )
            {
                // reset the value
                tElement->set_min_refinement_level( 0 );
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh_Base::clear_refinement_queue()
    {
        // loop over all levels
        for ( uint iLevel = 0; iLevel <= mMaxLevel; ++iLevel )
        {
            // collect elements from this level
            moris::Vector< Background_Element_Base* > tElements;

            this->collect_elements_on_level_including_aura( iLevel, tElements );

            // loop over all elements
            for ( Background_Element_Base* tElement : tElements )
            {
                // reset the value
                tElement->remove_from_refinement_queue();
            }
        }
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
