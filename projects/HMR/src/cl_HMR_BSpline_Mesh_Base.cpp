/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Mesh_Base.cpp
 *
 */

#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src

#include <fstream>

#include "HMR_Tools.hpp"       //HMR/src
#include "cl_Stopwatch.hpp"    //CHR/src
#include "cl_Matrix.hpp"       //LINALG/src
#include "fn_unique.hpp"       //LINALG/src
#include "cl_Map.hpp"
#include "cl_Tracer.hpp"
#include "fn_sum.hpp"

namespace moris::hmr
{
    //------------------------------------------------------------------------------
    //   public:
    //------------------------------------------------------------------------------

    BSpline_Mesh_Base::BSpline_Mesh_Base(
            const Parameters*     aParameters,
            Background_Mesh_Base* aBackgroundMesh,
            uint                  aOrder,
            uint                  aActivationPattern,
            uint                  aNumberOfBasesPerElement )
            : Mesh_Base( aParameters,
                    aBackgroundMesh,
                    aOrder,
                    aActivationPattern,
                    aNumberOfBasesPerElement )
    {
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::update_mesh()
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "B-Spline Mesh #" + std::to_string( this->get_index() ), "Update" );

        // activate pattern on background mesh
        this->select_activation_pattern();

        // tidy up memory
        this->delete_pointers();

        // create B-Spline Elements from Background Elements
        this->create_elements();

        // initialize basis and associate them to the elements
        this->create_basis();

        // this->use_only_basis_in_frame();
        this->calculate_basis_ids();

        // and determine the node ownership to be the smallest proc ID
        // of all elements connected to this node
        this->guess_basis_ownership();

        // Make sure that node ownership is correct. Correct otherwise.
        this->confirm_basis_ownership();

        // write all active basis into a container
        this->collect_active_and_refined_basis();

        // #ifdef MORIS_HAVE_DEBUG
        //             MORIS_LOG_WARNING("Sanity check for Bspline basis Ids and ownership will be performed. This might slow down the execution significantly. \n");
        //             this->sanity_check_for_ids_and_ownership();
        // #endif

        // determine indices of active and flagged basis
        // fixme: try Lagrange to B-Spline distance > 1 works if this is uncommented
        // this->calculate_basis_indices();

        // update element indices ( not needed so far )
        // this->update_element_indices();

    }    // end function: BSpline_Mesh_Base::update_mesh()

    //------------------------------------------------------------------------------

    bool
    BSpline_Mesh_Base::test_sanity()
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "B-Spline Mesh #" + std::to_string( this->get_index() ), "Perform sanity test" );

        this->calculate_basis_coordinates();

        // get parents for each BF
        this->link_bases_to_parents();

        // statement 0 : a BF can not be active and refined at the same time
        bool tTestForStateContradiction = true;

        // statement 1 : BF on top level must be active or refined
        bool tTestTopLevelState = true;

        // statement 2 : a BF that is active must have at least one refined parent
        bool tHaveRefinedParent = true;

        // FIXME: tests 2 and 3 are not sufficient in parallel
        // statement 3 : a BF must be deactivated if all parents are active
        bool tDeactivatedTest = true;

        // statement 4 : a BF that is refined must have at least one active descendant
        bool tRefinedHasActiveChild = true;

        // loop over all basis
        for ( auto tBasisFunction : mAllBasisOnProc )
        {
            // the statements
            if ( tBasisFunction->is_active() and tBasisFunction->is_refined() )
            {
                // contradiction is detected
                tTestForStateContradiction = false;
            }

            // the next steps only make sense if the basis function is actually used
            if ( tBasisFunction->is_used() )
            {
                // test level of basis
                if ( tBasisFunction->get_level() == 0 )
                {
                    // on the top level, only active or refined basis are allowed
                    tTestTopLevelState = tTestTopLevelState and ( tBasisFunction->is_active() or tBasisFunction->is_refined() );

                    /* if( par_rank() == 0 )
                    {
                        const real * tXY= tBasisFunction->get_xyz();
                        std::cout << "Basis ( " << tXY[ 0 ] << ", " << tXY[ 1 ] << " ) ["
                                << tBasisFunction->get_level() << "] "
                                << tBasisFunction->is_active()  << " " << tBasisFunction->is_refined() << " "
                               << tBasisFunction->get_memory_index() << std::endl;
                    } */
                }
                else
                {
                    // parent tests can only be done for higher level basis
                    uint tNumberOfParents = tBasisFunction->get_number_of_parents();

                    // this flag  is needed for statement 2
                    bool tRefinedParentFlag = false;

                    // this statement is needed for statement 3
                    bool tAllParentsAreActive = true;

                    // loop over all parents
                    for ( uint k = 0; k < tNumberOfParents; ++k )
                    {
                        // get pointer to parent
                        Basis* tParent = tBasisFunction->get_parent( k );

                        // only test if parent is relevant for this mesh
                        if ( tParent->is_used() )
                        {
                            // test if parent is refined
                            if ( tParent->is_refined() )
                            {
                                // set refined parent flag for statement 2
                                tRefinedParentFlag = true;
                            }

                            // set active parent flag for statement 3
                            tAllParentsAreActive = tAllParentsAreActive and tParent->is_active();
                        }
                    }

                    // test for statement 2
                    if ( tBasisFunction->is_active() )
                    {
                        tHaveRefinedParent = tHaveRefinedParent and tRefinedParentFlag;
                    }

                    // test for statement 3
                    if ( tAllParentsAreActive )
                    {
                        tDeactivatedTest = tDeactivatedTest and ( !tBasisFunction->is_active() and !tBasisFunction->is_refined() );
                        if ( !( !tBasisFunction->is_active() and !tBasisFunction->is_refined() ) )
                        {
                            const real* tXY = tBasisFunction->get_xyz();

                            std::cout << "Active Basis: " << tBasisFunction->get_level() << " "
                                      << tXY[ 0 ] << " " << tXY[ 1 ] << std::endl;
                        }
                    }
                }

                // test for statement 4
                if ( tBasisFunction->is_refined() )
                {
                    // needed for descendant counter
                    tBasisFunction->flag_descendants();

                    // test how many descendants exist
                    luint tDescendantCounter = tBasisFunction->count_descendants();

                    // initialize container of descendants
                    Cell< Basis* > tChildren( tDescendantCounter, nullptr );

                    // reset basis counter
                    tDescendantCounter = 0;

                    // collect  descendants
                    tBasisFunction->collect_descendants( tChildren, tDescendantCounter );

                    // reset found flag
                    bool tFoundActiveChild = false;

                    // loop over all children
                    for ( auto tChild : tChildren )
                    {
                        // test if child is active
                        if ( tChild->is_active() )
                        {
                            // set flag
                            tFoundActiveChild = true;

                            // break loop
                            break;
                        }
                    }

                    // set statement 4
                    tRefinedHasActiveChild = tRefinedHasActiveChild and tFoundActiveChild;
                }
            }
        }

        // tidy up flag table
        this->unflag_all_basis();

        bool aPassedTest = tTestForStateContradiction and tTestTopLevelState and tHaveRefinedParent and tDeactivatedTest and tRefinedHasActiveChild;

        if ( !aPassedTest )
        {
            std::cout << "Failed sanity test results: "
                      << tTestForStateContradiction << " "
                      << tTestTopLevelState << " "
                      << tHaveRefinedParent << " "
                      << tDeactivatedTest << " "
                      << tRefinedHasActiveChild << std::endl;
        }

        return aPassedTest;
    }

    //------------------------------------------------------------------------------
    // private:
    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::create_basis()
    {
        // report on this operation
        MORIS_LOG_INFO( "Creating basis functions" );

        // basis on first level are created separately
        this->create_basis_on_level_zero();

        // connect basis to coarsest elements
        this->link_basis_to_elements_on_level_zero();

        // reset max level that has active basis
        mMaxLevel = 0;

        // ask background mesh for number of levels
        luint tNumberOfLevels = mBackgroundMesh->get_max_level();

        for ( uint iLevel = 0; iLevel <= tNumberOfLevels; ++iLevel )
        {
            // refinement of basis on this level if the corresponding element is refined.
            this->process_level( iLevel );
        }

        this->collect_basis();

        for ( Basis* tBasis : mAllBasisOnProc )
        {
            tBasis->delete_neighbor_container();
        }

        // this->calculate_basis_coordinates();
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::collect_active_and_refined_elements_from_level(
            uint              aLevel,
            Cell< Element* >& aElements )
    {
        // cell containing background elements on this level
        Cell< Background_Element_Base* > tBackgroundElements;

        // ask background mesh about elements on this level
        mBackgroundMesh->collect_elements_on_level_including_aura( aLevel,
                tBackgroundElements );

        // count Elements
        luint tElementCount = 0;

        for ( Background_Element_Base* tBackElement : tBackgroundElements )
        {
            if ( !tBackElement->is_neither_active_nor_refined( mActivationPattern ) )
            {
                tElementCount++;
            }
        }

        // allocate cell
        aElements.resize( tElementCount, nullptr );

        // reset counter
        tElementCount = 0;
        for ( Background_Element_Base* tBackElement : tBackgroundElements )
        {
            if ( !tBackElement->is_neither_active_nor_refined( mActivationPattern ) )
            {
                aElements( tElementCount++ ) = mAllElementsOnProc( tBackElement->get_memory_index() );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::process_level( uint aLevel )
    {
        Cell< Element* > tElementsOnThisLevel;
        this->collect_active_and_refined_elements_from_level( aLevel, tElementsOnThisLevel );

        // initialize list of all basis functions on the current refinement level
        Cell< Basis* > tBFsOnThisLevel;

        // collect basis from given level
        this->preprocess_bases_from_level(
                tElementsOnThisLevel,
                tBFsOnThisLevel );

        // determine state of each basis
        this->determine_basis_state( tBFsOnThisLevel );

        // refine B-Spline mesh if this is not the last level
        if ( aLevel < mBackgroundMesh->get_max_level() )
        {
            for ( auto tElement : tElementsOnThisLevel )
            {
                // test if background element has children and is refined on pattern
                Background_Element_Base* tBgElem = tElement->get_background_element();
                if ( tBgElem->has_children() && tElement->is_refined() )
                {
                    // refine B-Spline element
                    auto tElemMemIndex = tElement->get_memory_index();
                    mAllElementsOnProc( tElemMemIndex )->refine( mAllElementsOnProc );
                }
            }
        }

        // delete list of basis functions // TODO: why?, there shouldn't be a new called underneath
        tBFsOnThisLevel.clear();
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::collect_basis()
    {
        // loop over all coarsest basis
        for ( auto tBasis : mAllCoarsestBasisOnProc )
        {
            // collect descendants
            tBasis->flag_descendants();
        }

        // loop over all coarsest basis
        mNumberOfAllBasis = 0;
        for ( auto tBasis : mAllCoarsestBasisOnProc )
        {
            // collect descendants
            mNumberOfAllBasis += tBasis->count_descendants();
        }

        // assign memory for cell
        mAllBasisOnProc.resize( mNumberOfAllBasis, nullptr );

        // reset counter
        luint tDescendantCount = 0;

        // loop over all coarsest basis
        for ( auto tBasis : mAllCoarsestBasisOnProc )
        {
            // collect descendants
            tBasis->collect_descendants( mAllBasisOnProc, tDescendantCount );
        }

        // reset counter
        luint tBasisIndex = 0;

        // unflag all basis and set index
        for ( auto tBasis : mAllBasisOnProc )
        {
            // unflag basis
            tBasis->unflag();

            // set memory index
            tBasis->set_memory_index( tBasisIndex++ );
        }

        // This function seems not to be needed, however for debugging purposes the basis can be stored in a vtk file.
        // In this case the information is needed
        // calculate basis coordinates
        //            this->calculate_basis_coordinates();
    }

    //------------------------------------------------------------------------------

    /* void
    BSpline_Mesh_Base::use_only_basis_in_frame()
    {
        // unuse all basis
        for( auto tBasis : mAllBasisOnProc )
        {
            tBasis->unuse();
        }

        uint tMyRank = par_rank();

        for( auto tElement : mAllElementsOnProc )
        {
            if ( tElement->get_owner() == tMyRank )
            {
                for( uint k=0; k<mNumberOfBasesPerElement; ++k )
                {
                    tElement->get_basis( k )->use();
                }
            }
        }
    } */

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::calculate_basis_ids()
    {
        // loop over all basis
        for ( auto tBasis : mAllBasisOnProc )
        {
            // get position of basis
            const luint* tIJK = tBasis->get_ijk();

            // calc id and write into basis
            uint  tLevel   = tBasis->get_level();
            luint tBasisId = this->calculate_basis_id( tLevel, tIJK );
            tBasis->set_domain_id( tBasisId );
        }
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::calculate_basis_indices( const Matrix< IdMat >& aCommTable )
    {
        // report on this operation
        MORIS_LOG_INFO( "B-Spline Mesh #%i: Computing basis function indices", this->get_index() );

        // get my rank
        moris_id tMyRank = par_rank();

        // get number of ranks
        uint tNumberOfProcs = par_size();

        // clean container
        mIndexedBasis.clear();

        // special function for multigrid
        this->flag_refined_basis_of_owned_elements();

        // reset all indices
        for ( Basis* tBasis : mAllBasisOnProc )
        {
            tBasis->set_local_index( gNoEntityID );
            tBasis->set_domain_index( gNoEntityID );
        }

        // synchronize flags if we are in parallel mode
        if ( tNumberOfProcs > 1 )
        {
            // Step 1: make sure that flags are synchronized
            this->synchronize_flags( aCommTable );
        }

        // reset all indices on the proc
        for ( Basis* tBasis : mAllBasisOnProc )
        {
            tBasis->set_local_index( gNoIndex );
            tBasis->set_domain_index( gNoID );
        }

        // counter for basis
        luint tBasisIndex = 0;

        // set local index of basis
        for ( Basis* tBasis : mActiveBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                // set index of basis
                tBasis->set_local_index( tBasisIndex++ );
            }
        }

        if ( mParameters->use_multigrid() )
        {
            for ( Basis* tBasis : mRefinedBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    // set index of basis
                    tBasis->set_local_index( tBasisIndex++ );
                }
            }
        }

        // allocate container
        mIndexedBasis.resize( tBasisIndex, nullptr );

        // reset counter
        tBasisIndex = 0;

        // copy indexed basis into container
        for ( Basis* tBasis : mActiveBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                mIndexedBasis( tBasisIndex++ ) = tBasis;
            }
        }

        if ( mParameters->use_multigrid() )
        {
            for ( Basis* tBasis : mRefinedBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    mIndexedBasis( tBasisIndex++ ) = tBasis;
                }
            }
        }

        if ( tNumberOfProcs == 1 )
        {
            // reset counter
            tBasisIndex = 0;

            for ( Basis* tBasis : mActiveBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    // set index of basis
                    tBasis->set_domain_index( tBasisIndex++ );
                }
            }

            if ( mParameters->use_multigrid() )
            {
                for ( Basis* tBasis : mRefinedBasisOnProc )
                {
                    if ( tBasis->is_flagged() )
                    {
                        // set index of basis
                        tBasis->set_domain_index( tBasisIndex++ );
                    }
                }
            }
        }    // end serial
        else
        {
            // Step 3: count flagged basis that are owned

            // reset counters
            luint tActiveCount  = 0;
            luint tRefinedCount = 0;

            // domain indices (= MTK IDs) loop over all basis
            for ( Basis* tBasis : mActiveBasisOnProc )
            {
                // test if basis is active, flagged and owned
                if ( tBasis->get_owner() == tMyRank and tBasis->is_flagged() )
                {
                    tBasis->set_domain_index( tActiveCount++ );
                }
            }

            if ( mParameters->use_multigrid() )
            {
                for ( Basis* tBasis : mRefinedBasisOnProc )
                {
                    // test if basis is active, flagged and owned
                    if ( tBasis->get_owner() == tMyRank and tBasis->is_flagged() )
                    {
                        tBasis->set_domain_index( tRefinedCount++ );    // FIXME should this be active count too?
                    }
                }
            }
            // - - - - - - - - - - - - - - - -

            // Step 4: communicate offset and add to domain index

            // communicate number of owned and active basis with other procs
            Matrix< DDLUMat > tActiveBasisCount;

            comm_gather_and_broadcast( tActiveCount, tActiveBasisCount );

            // get my offset
            moris_id tMyActiveOffset = 0;

            for ( moris_id p = 1; p <= tMyRank; ++p )
            {
                tMyActiveOffset += tActiveBasisCount( p - 1 );
            }

            Matrix< DDLUMat > tRefinedBasisCount;
            moris_id          tMyRefinedOffset = sum( tActiveBasisCount );

            if ( mParameters->use_multigrid() )
            {
                comm_gather_and_broadcast( tRefinedCount, tRefinedBasisCount );
                for ( moris_id p = 1; p <= tMyRank; ++p )
                {
                    tMyRefinedOffset += tRefinedBasisCount( p - 1 );
                }
            }

            // reset basis counter
            tActiveBasisCount.fill( 0 );

            // loop over all basis
            for ( Basis* tBasis : mActiveBasisOnProc )
            {
                // test if basis is flagged
                if ( tBasis->is_flagged() )
                {
                    // get owner of basis
                    auto tOwner = tBasis->get_owner();

                    // test if basis is mine
                    if ( tOwner == tMyRank )
                    {
                        tBasis->set_domain_index( tBasis->get_hmr_index() + tMyActiveOffset );
                    }
                    else
                    {
                        // increment basis counter per proc
                        ++tActiveBasisCount( tOwner );
                    }
                }
            }

            if ( mParameters->use_multigrid() )
            {
                tRefinedBasisCount.fill( 0 );
                for ( Basis* tBasis : mRefinedBasisOnProc )
                {
                    // test if basis is flagged
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is mine
                        if ( tOwner == tMyRank )
                        {
                            tBasis->set_domain_index( tBasis->get_hmr_index() + tMyRefinedOffset );
                        }
                        else
                        {
                            // increment basis counter per proc
                            ++tRefinedBasisCount( tOwner );
                        }
                    }
                }
            }

            // - - - - - - - - - - - - - - - -

            // Step 5: create map for communication
            Matrix< DDUMat > tProcIndices( tNumberOfProcs, 1, tNumberOfProcs );

            uint tCommLength = aCommTable.length();

            for ( uint k = 0; k < tCommLength; ++k )
            {
                tProcIndices( aCommTable( k ) ) = k;
            }

            // - - - - - - - - - - - - - - - -

            // Step 6: allocate memory for communication lists

            // dummy matrices for cells to send
            Matrix< DDLUMat > tEmptyLuint;
            Matrix< DDUMat >  tEmptyUint;

            // create cells for basis and element indices to send
            Cell< Matrix< DDLUMat > > tSendIndex( tCommLength, tEmptyLuint );
            Cell< Matrix< DDUMat > >  tSendBasis( tCommLength, tEmptyUint );
            Cell< Matrix< DDUMat > >  tSendPedigree( tCommLength, tEmptyUint );

            // assign memory for Index and Basis
            for ( uint p = 0; p < tCommLength; ++p )
            {
                luint tNumberOfBasis;

                if ( mParameters->use_multigrid() )
                {
                    tNumberOfBasis = tActiveBasisCount( aCommTable( p ) ) + tRefinedBasisCount( aCommTable( p ) );
                }
                else
                {
                    // get number of basis
                    tNumberOfBasis = tActiveBasisCount( aCommTable( p ) );
                }
                if ( tNumberOfBasis > 0 )
                {
                    ++tBasisIndex;
                    tSendIndex( p ).set_size( tNumberOfBasis, 1 );
                    tSendBasis( p ).set_size( tNumberOfBasis, 1 );
                }
            }

            // - - - - - - - - - - - - - - - -

            // Step 7: create lists with basis of which index is requested

            // reset counter
            tBasisIndex = 0;

            // reset counter
            Matrix< DDLUMat > tProcCount( tCommLength, 1, 0 );

            // loop over all basis
            for ( Basis* tBasis : mActiveBasisOnProc )
            {
                // test if basis is active
                if ( tBasis->is_flagged() )
                {
                    // get owner of basis
                    auto tOwner = tBasis->get_owner();

                    // test if basis is not mine
                    if ( tOwner != tMyRank )
                    {
                        // get index of owner
                        uint tProcIndex = tProcIndices( tOwner );

                        // pointer to element
                        this->get_reference_element_of_basis( tBasis,
                                tSendIndex( tProcIndex )( tProcCount( tProcIndex ) ),
                                tSendBasis( tProcIndex )( tProcCount( tProcIndex ) ) );

                        // increment counter
                        ++tProcCount( tProcIndex );
                    }
                }
            }

            if ( mParameters->use_multigrid() )
            {
                // loop over all basis
                for ( Basis* tBasis : mRefinedBasisOnProc )
                {
                    // test if basis is active
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is not mine
                        if ( tOwner != tMyRank )
                        {
                            // get index of owner
                            uint tIndex = tProcIndices( tOwner );

                            // pointer to element
                            this->get_reference_element_of_basis( tBasis,
                                    tSendIndex( tIndex )( tProcCount( tIndex ) ),
                                    tSendBasis( tIndex )( tProcCount( tIndex ) ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }
            }
            // local basis IDs received by other procs
            Cell< Matrix< DDUMat > > tReceiveBasis( tCommLength, tEmptyUint );

            // communicate local basis indices to request
            communicate_mats( aCommTable,
                    tSendBasis,
                    tReceiveBasis );

            // free memory
            tSendBasis.clear();

            // now we need to determine the memory needed for the
            // element pedigree paths

            // reset counter
            tProcCount.fill( 0 );

            // determine memory for pedigree path
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tSendIndex( p ).length();

                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    tProcCount( p ) += mAllElementsOnProc( tSendIndex( p )( k ) )
                                               ->get_background_element()
                                               ->get_length_of_pedigree_path();
                }
            }

            // encode pedigree paths
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tSendIndex( p ).length();

                // assign memory for path to send
                tSendPedigree( p ).set_size( tProcCount( p ), 1 );

                // reset counter
                tBasisIndex = 0;

                // loop over all elements
                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    // get pointer to element
                    Background_Element_Base* tElement = mAllElementsOnProc( tSendIndex( p )( k ) )
                                                                ->get_background_element();

                    // encode path and overwrite tSendElement with Ancestor Index
                    tElement->encode_pedigree_path( tSendIndex( p )( k ),
                            tSendPedigree( p ),
                            tBasisIndex );
                }
            }

            Cell< Matrix< DDLUMat > > tReceiveIndex( tCommLength, tEmptyLuint );
            Cell< Matrix< DDUMat > >  tReceivePedigree( tCommLength, tEmptyUint );

            // communicate ancestor IDs
            communicate_mats( aCommTable,
                    tSendIndex,
                    tReceiveIndex );

            // communicate pedigree paths
            communicate_mats( aCommTable,
                    tSendPedigree,
                    tReceivePedigree );

            // clear memory
            tSendPedigree.clear();

            // now we loop over all elements and determine the index of the requested basis
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tReceiveIndex( p ).length();

                // resize send index
                tSendIndex( p ).set_size( tNumberOfElements, 1 );

                // reset counter
                luint tPedigreeCount = 0;

                // loop over all elements
                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    // decode path and get pointer to element
                    Element* tElement = mAllElementsOnProc( mBackgroundMesh->decode_pedigree_path(
                                                                                   tReceiveIndex( p )( k ),
                                                                                   tReceivePedigree( p ),
                                                                                   tPedigreeCount )
                                                                    ->get_memory_index() );

                    // write index of requested basis into matrix
                    tSendIndex( p )( k ) = tElement->get_basis( tReceiveBasis( p )( k ) )
                                                   ->get_hmr_index();
                }
            }

            // clear memory
            tReceivePedigree.clear();
            tReceiveBasis.clear();
            tReceiveIndex.clear();

            // communicate requested indices back to original proc
            communicate_mats(
                    aCommTable,
                    tSendIndex,
                    tReceiveIndex );

            // clear memory
            tSendIndex.clear();

            // finally, we can set the indices of the unknown basis

            // reset counter
            tProcCount.fill( 0 );

            // loop over all basis
            for ( auto tBasis : mActiveBasisOnProc )
            {
                // test if basis is flagged
                if ( tBasis->is_flagged() )
                {
                    // get owner of basis
                    auto tOwner = tBasis->get_owner();

                    // test if basis is mine
                    if ( tOwner != tMyRank )
                    {
                        // get index of owner
                        uint tIndex = tProcIndices( tOwner );

                        // get counter
                        tBasisIndex = tProcCount( tIndex );

                        // write index into Communication List as is
                        tBasis->set_domain_index( tReceiveIndex( tIndex )( tBasisIndex ) );

                        // increment counter
                        ++tProcCount( tIndex );
                    }
                }
            }

            if ( mParameters->use_multigrid() )
            {
                for ( auto tBasis : mRefinedBasisOnProc )
                {
                    // test if basis is flagged
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is mine
                        if ( tOwner != tMyRank )
                        {
                            // get index of owner
                            uint tIndex = tProcIndices( tOwner );

                            // get counter
                            tBasisIndex = tProcCount( tIndex );

                            // write index into ba Communication List as is
                            tBasis->set_domain_index( tReceiveIndex( tIndex )( tBasisIndex ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }
            }
            // perform a small sanity test :
            tBasisIndex = 0;

            // loop over all basis
            for ( auto tBasis : mAllBasisOnProc )
            {
                // test if basis is used, active and has no id
                if ( tBasis->is_flagged()
                        and tBasis->is_active()
                        and tBasis->get_hmr_index() == gNoEntityID )
                {
                    std::cout << par_rank() << " bad basis " << tBasis->get_hmr_id() << " " << tBasis->get_owner() << std::endl;

                    // increment counter
                    ++tBasisIndex;
                }
            }

            MORIS_ERROR( tBasisIndex == 0, "%s ERROR.\n               Could not identify indices of %lu basis.\n               This might happen if a proc uses an active basis that does not belong to\n               itself or any direct neighbor. Suggestion: use denser mesh on top level.\n\n", proc_string().c_str(), (long unsigned int)tBasisIndex );
        }    // end if parallel

             // insert parents if we are in multigrid
        if ( mParameters->use_multigrid() )
        {
            // get parents for each basis
            this->link_bases_to_parents();
        }

        /*
        #if defined(DEBUG)
                // Test sanity #CHRISTIAN
                luint tNumberOfBSplines = this->get_number_of_active_basis_on_proc();
                Matrix< DDLUMat > tIDs( tNumberOfBSplines, 1 );
                for( uint k=0; k<tNumberOfBSplines; ++k )
                {
                    tIDs( k ) = this->get_active_basis( k )->get_hmr_id();
                }

                Matrix< DDLUMat > tIDsUnique;
                unique( tIDs, tIDsUnique );
                std::cout << "B-Spline mesh " <<
                        this->get_order() << " " << this->get_activation_pattern()
                        << " " << tNumberOfBSplines << std::endl;

                MORIS_ERROR( tIDsUnique.length() == tNumberOfBSplines, "Duplicate Basis created" );
                MORIS_ERROR( tNumberOfActiveBasis == tNumberOfBSplines, "Number of Basis does not match" );
        #endif
        */
    } // end function::BSpline_Mesh_Base::calculate_basis_indices()

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::synchronize_flags( const Matrix< IdMat >& aCommTable )
    {
        // get number of ranks
        uint tNumberOfProcs = par_size();

        if ( tNumberOfProcs > 1 )
        {
            // get my rank
            moris_id tMyRank = par_rank();

            // length of communication table
            uint tCommLength = aCommTable.length();

            // count number of basis per proc
            Matrix< DDUMat > tBasisCount( tNumberOfProcs, 1, 0 );

            // loop over all active basis
            for ( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    // increment basis counter
                    ++tBasisCount( tBasis->get_owner() );
                }
            }

            // matrix to check if communication table makes sense
            Matrix< DDUMat > tBasisCommCheck( tBasisCount );

            // make sure that communication table is sane
            for ( uint Ik = 0; Ik < tCommLength; ++Ik )
            {
                tBasisCommCheck( aCommTable( Ik ) ) = 0;
            }

            // reset my own counter
            tBasisCommCheck( par_rank() ) = 0;

            if ( tBasisCommCheck.max() != 0 )
            {
                std::cout << "Processor " << par_rank() << std::endl;
                print( aCommTable, "CommTable" );
                print( tBasisCommCheck, "CommCheck" );
            }

            MORIS_ERROR( tBasisCommCheck.max() == 0, "synchronize_flags: error in communication table" );

            // dummy matrices for cells to send
            Matrix< DDLUMat > tEmptyLuint;
            Matrix< DDUMat >  tEmptyUint;

            // create cells for basis and element indices to send
            Cell< Matrix< DDLUMat > > tSendIndex( tCommLength, tEmptyLuint );
            Cell< Matrix< DDUMat > >  tSendBasis( tCommLength, tEmptyUint );
            Cell< Matrix< DDUMat > >  tSendPedigree( tCommLength, tEmptyUint );

            // assign memory for Index and Basis
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of basis
                luint tNumberOfBasis = tBasisCount( aCommTable( p ) );

                if ( tNumberOfBasis > 0 )
                {
                    tSendIndex( p ).set_size( tNumberOfBasis, 1 );
                    tSendBasis( p ).set_size( tNumberOfBasis, 1 );
                }
            }

            // reset counter
            Matrix< DDLUMat > tProcCount( tCommLength, 1, 0 );

            // this table converts the proc id to an index
            Matrix< DDUMat > tProcIndices( tNumberOfProcs, 1, tNumberOfProcs );
            for ( uint k = 0; k < tCommLength; ++k )
            {
                tProcIndices( aCommTable( k ) ) = k;
            }

            // loop over all basis
            for ( auto tBasis : mAllBasisOnProc )
            {
                // test if basis is active
                if ( tBasis->is_flagged() )
                {
                    // get owner of basis
                    auto tOwner = tBasis->get_owner();

                    // test if basis is not mine
                    if ( tOwner != tMyRank )
                    {
                        // get index of owner
                        uint tProcIndex = tProcIndices( tOwner );

                        // pointer to element
                        this->get_reference_element_of_basis( tBasis,
                                tSendIndex( tProcIndex )( tProcCount( tProcIndex ) ),
                                tSendBasis( tProcIndex )( tProcCount( tProcIndex ) ) );

                        // increment counter
                        ++tProcCount( tProcIndex );
                    }
                }
            }

            // local basis IDs received by other procs
            Cell< Matrix< DDUMat > > tReceiveBasis( tCommLength, tEmptyUint );

            // communicate local basis IDs to request
            communicate_mats( aCommTable,
                    tSendBasis,
                    tReceiveBasis );

            // free memory
            tSendBasis.clear();

            // now we need to determine the memory needed for the
            // element pedigree paths

            // reset counter
            tProcCount.fill( 0 );

            // determine memory for pedigree path
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tSendIndex( p ).length();

                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    tProcCount( p ) += mAllElementsOnProc( tSendIndex( p )( k ) )
                                               ->get_background_element()
                                               ->get_length_of_pedigree_path();
                }
            }

            // encode pedigree paths
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tSendIndex( p ).length();

                // assign memory for path to send
                tSendPedigree( p ).set_size( tProcCount( p ), 1 );

                // reset counter
                luint tPedigreeCount = 0;

                // loop over all elements
                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    // get pointer to element
                    Background_Element_Base* tElement = mAllElementsOnProc( tSendIndex( p )( k ) )
                                                                ->get_background_element();

                    // encode path and overwrite tSendElement with Ancestor Index
                    tElement->encode_pedigree_path( tSendIndex( p )( k ),
                            tSendPedigree( p ),
                            tPedigreeCount );
                }
            }

            Cell< Matrix< DDLUMat > > tReceiveIndex( tCommLength, tEmptyLuint );
            Cell< Matrix< DDUMat > >  tReceivePedigree( tCommLength, tEmptyUint );

            // communicate ancestor IDs
            communicate_mats( aCommTable,
                    tSendIndex,
                    tReceiveIndex );

            // communicate pedigree paths
            communicate_mats( aCommTable,
                    tSendPedigree,
                    tReceivePedigree );

            // clear memory
            tSendPedigree.clear();

            // now we loop over all elements and determine the index of the requested basis
            for ( uint p = 0; p < tCommLength; ++p )
            {
                // get number of elements
                luint tNumberOfElements = tReceiveIndex( p ).length();

                // resize send index
                tSendIndex( p ).set_size( tNumberOfElements, 1 );

                // reset counter
                luint tPedigreeCount = 0;

                // loop over all elements
                for ( luint k = 0; k < tNumberOfElements; ++k )
                {
                    // decode path and get pointer to element
                    Element* tElement = mAllElementsOnProc( mBackgroundMesh->decode_pedigree_path(
                                                                                   tReceiveIndex( p )( k ),
                                                                                   tReceivePedigree( p ),
                                                                                   tPedigreeCount )
                                                                    ->get_memory_index() );

                    // now we flag this basis
                    tElement->get_basis( tReceiveBasis( p )( k ) )->flag();
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::collect_active_and_refined_basis()
    {
        // report on this operation
        MORIS_LOG_INFO( "Collect active and refined basis functions" );

        // reset counter
        mNumberOfActiveBasisOnProc  = 0;
        mNumberOfRefinedBasisOnProc = 0;
        mMaxLevel                   = 0;

        // collect additional basis functions if multi-grid is being used, if not see below
        if ( mParameters->use_multigrid() )
        {
            // count active basis on proc
            for ( auto tBasis : mAllBasisOnProc )
            {
                // reset index
                tBasis->set_active_index( gNoEntityID );

                // count basis
                if ( tBasis->is_used() )
                {
                    if ( tBasis->is_active() )
                    {
                        ++mNumberOfActiveBasisOnProc;
                        mMaxLevel = std::max( tBasis->get_level(), mMaxLevel );
                    }
                    else if ( tBasis->is_refined() )
                    {
                        ++mNumberOfRefinedBasisOnProc;
                    }
                }
            }

            // reserve memory
            mActiveBasisOnProc.resize( mNumberOfActiveBasisOnProc, nullptr );
            mRefinedBasisOnProc.resize( mNumberOfRefinedBasisOnProc, nullptr );

            // reset counters
            mNumberOfActiveBasisOnProc  = 0;
            mNumberOfRefinedBasisOnProc = 0;

            // count active basis on proc
            for ( auto tBasis : mAllBasisOnProc )
            {
                // count basis
                if ( tBasis->is_used() )
                {
                    if ( tBasis->is_active() )
                    {
                        tBasis->set_active_index( mNumberOfActiveBasisOnProc );    // FIXME

                        mActiveBasisOnProc( mNumberOfActiveBasisOnProc++ ) = tBasis;
                    }
                    else if ( tBasis->is_refined() )
                    {
                        mRefinedBasisOnProc( mNumberOfRefinedBasisOnProc++ ) = tBasis;
                    }
                }
            }
        }    // end if: geometric multi-grid is used

        // if: multi-grid is NOT used
        else
        {
            // count active basis on proc
            for ( auto tBasis : mAllBasisOnProc )
            {
                // reset index
                tBasis->set_active_index( gNoEntityID );

                // count basis
                if ( tBasis->is_used() )
                {
                    if ( tBasis->is_active() )
                    {
                        ++mNumberOfActiveBasisOnProc;
                        mMaxLevel = std::max( tBasis->get_level(), mMaxLevel );
                    }
                }
            }

            // reserve memory
            mActiveBasisOnProc.resize( mNumberOfActiveBasisOnProc, nullptr );

            // initialize counter
            mNumberOfActiveBasisOnProc = 0;

            // populate container
            for ( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_active() and tBasis->is_used() )
                {
                    tBasis->set_active_index( mNumberOfActiveBasisOnProc );

                    mActiveBasisOnProc( mNumberOfActiveBasisOnProc++ ) = tBasis;
                }
            }

        }    // end if: multi-grid is NOT used

    }        // end function: hmr::BSpline_Mesh_Base::collect_active_and_refined_basis()

    //------------------------------------------------------------------------------
    bool
    BSpline_Mesh_Base::test_for_double_basis()
    {
        luint tNumberOfBasis = mAllBasisOnProc.size();

        Matrix< DDLUMat > tBasisIDs( tNumberOfBasis, 1 );

        // initialize counter
        luint tBasisCount = 0;

        // populate container
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_used() and ( tBasis->is_active() or tBasis->is_refined() ) )
            // still have doubled basis if they are not used
            // by proc however, that does not matter since
            // they are no DOFs
            {
                tBasisIDs( tBasisCount++ ) = tBasis->get_hmr_id();
            }
        }

        tBasisIDs.resize( tBasisCount, 1 );

        // make basis unique
        Matrix< DDLUMat > tBasisUniqueIDs;
        unique( tBasisIDs, tBasisUniqueIDs );

        return tBasisUniqueIDs.length() == tBasisCount;
    }

    //------------------------------------------------------------------------------

    void
    BSpline_Mesh_Base::save_to_vtk( const std::string& aFilePath )
    {
        // log & trace this operation
        Tracer tTracer( "HMR", "B-Spline Mesh #" + std::to_string( this->get_index() ), "Save to VTK" );

        this->calculate_basis_coordinates();

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

        // unflag all bases
        this->unflag_all_basis();

        // get my rank
        moris_id tMyRank = par_rank();

        // Flag all bases of elements on this proc
        for ( auto tElement : mAllElementsOnProc )
        {
            if ( tElement->get_owner() == tMyRank )
            {
                tElement->flag_all_bases();
            }
        }

        // count flagged basis
        uint tNumberOfBasis = 0;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                // increment counter
                ++tNumberOfBasis;
            }
        }

        // write node data
        tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        tFile << "POINTS " << tNumberOfBasis << " float" << std::endl;

        // ask settings for number of dimensions
        auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

        if ( tNumberOfDimensions == 2 )
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    // get coordinate from basis
                    const real* tXY = tBasis->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float)tXY[ 0 ] );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                    tFChar = swap_byte_endian( (float)tXY[ 1 ] );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                    tFChar = swap_byte_endian( (float)0 );
                    // tFChar = swap_byte_endian( (float) tBasis->get_level() );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                }
            }
        }
        else if ( tNumberOfDimensions == 3 )
        {
            for ( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_flagged() )
                {
                    // get coordinate from node
                    const real* tXYZ = tBasis->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float)tXYZ[ 0 ] );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                    tFChar = swap_byte_endian( (float)tXYZ[ 1 ] );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                    tFChar = swap_byte_endian( (float)tXYZ[ 2 ] );
                    tFile.write( (char*)&tFChar, sizeof( float ) );
                }
            }
        }

        tFile << std::endl;

        // write each basis as its own element
        tFile << "CELLS " << tNumberOfBasis << " "
              << 2 * tNumberOfBasis << std::endl;

        int tOne = swap_byte_endian( (int)1 );

        // reset counter
        int tBasisCount = 0;

        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( tBasisCount );
                tFile.write( (char*)&tOne, sizeof( int ) );
                tFile.write( (char*)&tIChar, sizeof( int ) );

                ++tBasisCount;
            }
        }

        // write cell types
        tFile << "CELL_TYPES " << tNumberOfBasis << std::endl;
        tIChar = swap_byte_endian( (int)2 );
        for ( luint k = 0; k < tNumberOfBasis; ++k )
        {
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        // write node data
        tFile << "POINT_DATA " << tNumberOfBasis << std::endl;

        // write state
        tFile << "SCALARS STATE int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                // state flag
                int tState = 0;

                // get state from basis
                if ( tBasis->is_active() )
                {
                    tState = 2;
                }
                else if ( tBasis->is_refined() )
                {
                    tState = 1;
                }
                tIChar = swap_byte_endian( tState );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write basis ID
        tFile << "SCALARS ID int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_id() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write internal basis ID
        tFile << "SCALARS HMR_ID int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_hmr_id() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write basis level
        tFile << "SCALARS LEVEL int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_level() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write basis owner
        tFile << "SCALARS OWNER int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_owner() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write memory index
        tFile << "SCALARS MEMORY_INDEX int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_memory_index() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write mtk index
        tFile << "SCALARS MTK_INDEX int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                tIChar = swap_byte_endian( (int)tBasis->get_index() );

                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << std::endl;

        // write active index
        tFile << "SCALARS ACTIVE_INDEX int" << std::endl;
        tFile << "LOOKUP_TABLE default" << std::endl;
        for ( auto tBasis : mAllBasisOnProc )
        {
            if ( tBasis->is_flagged() )
            {
                if ( tBasis->is_active() )
                {
                    tIChar = swap_byte_endian( (int)tBasis->get_active_index() );
                }
                else
                {
                    tIChar = swap_byte_endian( (int)-1 );
                }
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }

        tFile << std::endl;

        // close the output file
        tFile.close();

        // unflag all bases
        this->unflag_all_basis();
    }

    //------------------------------------------------------------------------------

    Matrix< DDSMat >
    BSpline_Mesh_Base::get_children_ind_for_basis( const moris::sint aParentBasisIndex )
    {
        // get basis pointer
        Basis* tBasis = this->get_basis_by_index( aParentBasisIndex );

        // get child indices
        Matrix< DDSMat > tBasisLocalChildInds;

        tBasis->get_basis_local_child_inds( tBasisLocalChildInds );

        uint tNumberOfChildren = tBasisLocalChildInds.length();

        Matrix< DDSMat > tChildInds( tNumberOfChildren, 1 );

        for ( uint k = 0; k < tNumberOfChildren; ++k )
        {
            tChildInds( k ) = tBasis->get_child( tBasisLocalChildInds( k ) )->get_index();
        }

        return tChildInds;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    BSpline_Mesh_Base::get_children_weights_for_parent( const moris::sint aParentBasisIndex )
    {
        // get basis pointer
        Basis* tBasis = this->get_basis_by_index( aParentBasisIndex );

        // get child indices
        Matrix< DDSMat > tBasisLocalChildInds;

        tBasis->get_basis_local_child_inds( tBasisLocalChildInds );

        uint tNumberOfChildren = tBasisLocalChildInds.length();

        // create weights
        Matrix< DDRMat > tWeights( tNumberOfChildren, 1 );

        for ( uint k = 0; k < tNumberOfChildren; ++k )
        {
            tWeights( k ) = mChildStencil( tBasisLocalChildInds( k ) );
        }
        return tWeights;
    }

    //------------------------------------------------------------------------------

    uint
    BSpline_Mesh_Base::get_number_of_basis_connected_to_basis( const moris_index aIndex )
    {
        // get basis pointer
        Basis* tBasis = mActiveBasisOnProc( aIndex );

        // step 1: unflag all connected basis from connected elements

        // get number of connected elements
        uint tNumberOfElements = tBasis->get_element_counter();

        // loop over all elements of this basis
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            // get Cell of connected vertices
            moris::Cell< mtk::Vertex* > tVertices = tBasis->get_element( e )->get_vertex_pointers();

            // unflag these vertices
            for ( mtk::Vertex* tVertex : tVertices )
            {
                tVertex->unflag();
            }
        }

        // step 2: count connected vertices

        // reset counter
        uint tBasisCount = 0;

        // loop over all elements of this basis
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            // get Cell of connected vertices
            moris::Cell< mtk::Vertex* > tVertices = tBasis->get_element( e )->get_vertex_pointers();

            // unflag these vertices
            for ( mtk::Vertex* tVertex : tVertices )
            {
                if ( !tVertex->is_flagged() )
                {
                    tVertex->unflag();

                    ++tBasisCount;
                }
            }
        }

        // return counter
        return tBasisCount;
    }

}    // namespace moris::hmr
