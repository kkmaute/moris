/*
 * cl_HMR_BSpline_Mesh_Base.cpp
 *
 *  Created on: Okt 12, 2018
 *      Author: Schmidt
 */
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src

#include <fstream>

#include "HMR_Tools.hpp" //HMR/src
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "fn_unique.hpp" //LINALG/src
#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------
//   public:
//------------------------------------------------------------------------------

        BSpline_Mesh_Base::BSpline_Mesh_Base (
                const Parameters      * aParameters,
                Background_Mesh_Base  * aBackgroundMesh,
                const uint            & aOrder,
                const uint            & aActivationPattern ) : Mesh_Base( aParameters,
                                                                          aBackgroundMesh,
                                                                          aOrder,
                                                                          aActivationPattern ),
                                                               mNumberOfChildrenPerBasis( std::pow( aOrder + 2,
                                                                                          aParameters->get_number_of_dimensions() ) ),
                                                               mNumberOfElementsPerBasis( std::pow( aOrder+1,
                                                                                          aParameters->get_number_of_dimensions() ) )
        {
            this->calculate_child_stencil();

            //FIXME this function is equal to calculate_lookup_tables(). This function is called in the child class
            //this->calculate_basis_level_offset();
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::update_mesh()
        {
            // start timer
            tic tTimer;

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

//#ifdef DEBUG
//            MORIS_LOG_WARNING("Sanity check for Bspline basis Ids and ownership will be performed. This might slow down the execution significantly. \n");
//            this->sanity_check_for_ids_and_ownership();
//#endif

            // determine indices of active and flagged basis
            // fixme: try Lagrange to B-Spline distance > 1 works if this is uncommented
            //this->calculate_basis_indices();

            // update element indices ( not needed so far )
            // this->update_element_indices();

            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            MORIS_LOG_INFO( "%s Created B-Spline mesh of order %u on pattern %u.",
                    proc_string().c_str(),
                    ( unsigned int )      mOrder,
                    ( unsigned int )      mActivationPattern);

            MORIS_LOG_INFO( "Mesh has %lu  Elements and %lu basis in total.",
                    ( long unsigned int ) mNumberOfAllElementsOnProc,
                    ( long unsigned int ) mNumberOfAllBasis );

            MORIS_LOG_INFO( "Mesh uses %lu basis on proc.",
                    ( long unsigned int ) mNumberOfAllBasis );

            MORIS_LOG_INFO( "Mesh has %lu active basis on proc.",
                    ( long unsigned int ) mNumberOfActiveBasisOnProc );

            MORIS_LOG_INFO( "Creation took %5.3f seconds.",
                    ( double ) tElapsedTime / 1000 );
            MORIS_LOG_INFO( " " );


//            MORIS_LOG_INFO( "%s Created B-Spline mesh of order %u on pattern %u."
//                            "Mesh has %lu  Elements and %lu basis in total."
//                            "Mesh uses %lu basis on proc."
//                            "Mesh has %lu active basis on proc."
//                            "Creation took %5.3f seconds.",
//                    proc_string().c_str(),
//                    ( unsigned int )      mOrder,
//                    ( unsigned int )      mActivationPattern,
//                    ( long unsigned int ) mNumberOfAllElementsOnProc,
//                    ( long unsigned int ) mNumberOfAllBasis,
//                    ( long unsigned int ) mNumberOfBasis,
//                    ( long unsigned int ) mNumberOfActiveBasisOnProc,
//                    ( double ) tElapsedTime / 1000 );
        }

//------------------------------------------------------------------------------

        Basis * BSpline_Mesh_Base::get_coarsest_basis_by_ij( const luint & aI,
                                                             const luint & aJ )
        {
            MORIS_ASSERT( aI < mNumberOfCoarsestBasisOnProc[ 0 ] && aJ < mNumberOfCoarsestBasisOnProc[ 1 ],
                    "get_coarsest_basis_by_ij(), requested basis outside of the domain limits");

//            if ( aI < mNumberOfCoarsestBasisOnProc[ 0 ] && aJ < mNumberOfCoarsestBasisOnProc[ 1 ] )
//            {
                return mAllCoarsestBasisOnProc( aI + aJ*mNumberOfCoarsestBasisOnProc[ 0 ] );
//            }
//            else
//            {
//                return nullptr;
//            }
        }

//------------------------------------------------------------------------------

        Basis * BSpline_Mesh_Base::get_coarsest_basis_by_ijk( const luint & aI,
                                                              const luint & aJ,
                                                              const luint & aK )
        {
            MORIS_ASSERT( aI < mNumberOfCoarsestBasisOnProc[ 0 ] && aJ < mNumberOfCoarsestBasisOnProc[ 1 ] && aK < mNumberOfCoarsestBasisOnProc[ 2 ],
                    "get_coarsest_basis_by_ij(), requested basis outside of the domain limits");

//            if (    aI < mNumberOfCoarsestBasisOnProc[ 0 ] &&
//                    aJ < mNumberOfCoarsestBasisOnProc[ 1 ] &&
//                    aK < mNumberOfCoarsestBasisOnProc[ 2 ] )
//            {
                return mAllCoarsestBasisOnProc( aI + mNumberOfCoarsestBasisOnProc[ 0 ] *
                                              ( aJ + aK * mNumberOfCoarsestBasisOnProc[ 1 ] ) );
//            }
//            else
//            {
//                return nullptr;
//            }
        }

//------------------------------------------------------------------------------

        /**
         *
         * This test returns true if the activation pattern of the basis
         * seems correct. The rule is as follows
         *
         *  A basis is:
         *
         *      - active, if all connected elements are either active,
         *                refined or padding, and at least one element is active
         *
         *      - refined, if all connected elements are refined or padding
         *
         *      - deactive, if at least one element is deactive ( or does not exist)
         *
         *
         *  Form these rules, some statements ( see comments in source ) have been developed.
         *  Note that fulfilling the checked statements are NECESSARY, BUT NOT SUFFICIENT conditions.
         *
         *  This test is not meant to be run during runtime.
         */
        bool BSpline_Mesh_Base::test_sanity()
        {
            this->calculate_basis_coordinates();

            // start clock
            tic tTimer;

            // get parents for each basis
            this->link_basis_to_parents();

            // statement 0 : a basis can not be active and refined at the same time
            bool  tTestForStateContratiction = true;

            // statement 1 : basis on top level must be active or refined
            bool  tTestTopLevelState = true;

            // statement 2 : a basis that is active must have at least one refined parent
            bool tHaveRefinedParent = true;

            // FIXME: tests 2 and 3 are not sufficient in parallel
            // statement 3 : a basis must be deactive if all parents are active
            bool tDeactiveTest = true;

            // statement 4 : a basis that is refined must have at least one active descendant
            bool tRefinedHasActiveChild = true;

            // loop over all basis
            for( auto tBasis : mAllBasisOnProc )
            {
                // the statements
                if ( tBasis->is_active() && tBasis->is_refined() )
                {
                    // contradiciton is detected
                    tTestForStateContratiction = false;
                }

                // the next steps only make sense if the basis is actually used
                if ( tBasis->is_used() )
                {
                    // test level of basis
                    if( tBasis->get_level() == 0 )
                    {
                        // on the top level, only active or refined basis are allowed
                        tTestTopLevelState = tTestTopLevelState && ( tBasis->is_active() || tBasis->is_refined() );

                        /* if( par_rank() == 0 )
                        {
                            const real * tXY= tBasis->get_xyz();
                            std::cout << "Basis ( " << tXY[ 0 ] << ", " << tXY[ 1 ] << " ) ["
                                    << tBasis->get_level() << "] "
                                    << tBasis->is_active()  << " " << tBasis->is_refined() << " "
                                   << tBasis->get_memory_index() << std::endl;
                        } */
                    }
                    else
                    {
                        // parent tests can only be done for higher level basis
                        uint tNumberOfParents = tBasis->get_number_of_parents();

                        // this flag  is needed for statement 2
                        bool tRefinedParetFlag = false;

                        // this statement is needed for statement 3
                        bool tAllParentsAreActive = true;

                        // loop over all parents
                        for( uint k=0; k<tNumberOfParents; ++k )
                        {
                            // get pointer to parent
                            Basis* tParent = tBasis->get_parent( k );

                            // only test if parent is relevant for this mesh
                            if ( tParent->is_used() )
                            {
                                // test if parent is refined
                                if( tParent->is_refined() )
                                {
                                    // set refined parent flag for statement 2
                                    tRefinedParetFlag = true;
                                }

                                // set active parent flag for statement 3
                                tAllParentsAreActive = tAllParentsAreActive && tParent->is_active();
                            }
                        }

                        // test for statement 2
                        if ( tBasis->is_active() )
                        {
                            tHaveRefinedParent =  tHaveRefinedParent && tRefinedParetFlag;
                        }

                        // test for statement 3
                        if ( tAllParentsAreActive )
                        {
                            tDeactiveTest  = tDeactiveTest && ( ! tBasis->is_active() && ! tBasis->is_refined() );
                            if ( ! ( ! tBasis->is_active() && ! tBasis->is_refined() ) )
                            {
                                const real* tXY = tBasis->get_xyz();

                                std::cout << "Active Basis: " << tBasis->get_level() << " "
                                        << tXY[ 0 ] << " " << tXY[ 1 ] << std::endl;
                            }
                        }
                    }

                    // test for statement 4
                    if ( tBasis->is_refined() )
                    {

                        // reset basis counter
                        luint tCounter = 0;

                        // needed for descendant counter
                        tBasis->flag_descendants();

                        // test how many descendants exist
                        tBasis->count_descendants( tCounter );

                        // initialize container of descendants
                        Cell< Basis* > tChildren( tCounter, nullptr );

                        // reset basis counter
                        tCounter = 0;

                        // collect  descendants
                        tBasis->collect_descendants( tChildren, tCounter );

                        // reset foun flag
                        bool tFoundActiveChild = false;

                        // loop over all children
                        for( auto tChild: tChildren )
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
                        tRefinedHasActiveChild = tRefinedHasActiveChild && tFoundActiveChild;
                    }
                }
            }

            // tidy up flag table
            this->unflag_all_basis();

            bool aPassedTest = tTestForStateContratiction &&
                               tTestTopLevelState &&
                               tHaveRefinedParent &&
                               tDeactiveTest &&
                               tRefinedHasActiveChild ;

           // stop timer
           real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

           if ( aPassedTest )
           {
               MORIS_LOG_INFO( "%s Tested basis activation sanity.",
                       proc_string().c_str() );
               MORIS_LOG_INFO( "Test took %5.3f seconds.",
                       ( double ) tElapsedTime / 1000 );
               MORIS_LOG_INFO( "All tests passed.");
               MORIS_LOG_INFO( " " );

           }
           else
           {
               MORIS_LOG_INFO("%s Tested basis activation sanity.",
                       proc_string().c_str() );
               MORIS_LOG_INFO("Test took %5.3f seconds.",
                        ( double ) tElapsedTime / 1000 );
               MORIS_LOG_INFO( "AT LEAST ONE TEST FAILED.");
               MORIS_LOG_INFO( " " );

                   std::cout << "Test result: "
                             <<  tTestForStateContratiction << " "
                             << tTestTopLevelState << " "
                             << tHaveRefinedParent << " "
                             << tDeactiveTest << " "
                             << tRefinedHasActiveChild  << std::endl;

           }

            return aPassedTest;
        }

//------------------------------------------------------------------------------
// private:
//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::calculate_basis_level_offset()
        {
            // first level is zero
            mBasisLevelOffset[ 0 ] = 0;

            // get number of elements from Background Mesh
            Matrix< DDLUMat > tNumberOfElements = mBackgroundMesh->get_number_of_elements_per_direction();

            // get basis per direction
            luint tBasisPerDirection[ 3 ];

            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                tBasisPerDirection[ k ] = tNumberOfElements( k, 0 ) + mOrder;
            }

            // loop over all higher levels
            for( uint l=0; l<gMaxNumberOfLevels-1; ++l )
            {
                // number of basis on this level
                luint tBasisOnLastLevel = 1;

                // loop over all dimensions
                for (  uint k=0; k<mNumberOfDimensions; ++k )
                {
                    tBasisOnLastLevel       *= tBasisPerDirection[ k ];

                    tBasisPerDirection[ k ] *= 2;
                    tBasisPerDirection[ k ] += mOrder;
                }

                // increment delta

                // add offset to table
                mBasisLevelOffset[ l+1 ] = mBasisLevelOffset[ l ] + tBasisOnLastLevel;
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::create_basis()
        {
            // basis on first level are created separately
            this->create_basis_on_level_zero();

            // connect basis to coarsest elements
            this->link_basis_to_elements_on_level_zero();

            // reset max level that has active basis
            mMaxLevel = 0;

            // ask background mesh for number of levels
            luint tNumberOfLevels = mBackgroundMesh->get_max_level();

            for( uint l = 0; l <= tNumberOfLevels; ++l )
            {
                // refinement of basis on this level if the corresponding element is refined.
                this->process_level( l );
            }

            this->collect_basis();

            for( Basis * tBasis : mAllBasisOnProc )
            {
                tBasis->delete_neighbor_container();
            }

            //this->calculate_basis_coordinates();
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::create_basis_on_level_zero()
        {
            // ask mesh for relevant ijk positions
            Matrix< DDLUMat > tIJK = mBackgroundMesh ->get_number_of_elements_per_direction_on_proc();

            // initialize basis counter
            luint tCount = 0;

            if( mNumberOfDimensions == 2)
            {
                // unroll min and max i and j
                mNumberOfCoarsestBasisOnProc[ 0 ] = tIJK( 0, 0 ) + mOrder;

                mNumberOfCoarsestBasisOnProc[ 1 ] = tIJK( 1, 0 ) + mOrder ;

                // initialize array
                mAllCoarsestBasisOnProc.resize( mNumberOfCoarsestBasisOnProc[ 0 ] * mNumberOfCoarsestBasisOnProc[ 1 ],
                                                nullptr );

                // container with position to be passed to new basis
                luint tIJ[ 2 ];

                // loop over all j
                for( luint j = 0; j < mNumberOfCoarsestBasisOnProc[ 1 ]; ++j )
                {
                    // save j-position
                    tIJ[ 1 ] = j;

                    // loop over all i
                    for( luint i = 0; i < mNumberOfCoarsestBasisOnProc[ 0 ]; ++i )
                    {
                        // save i-position
                        tIJ[ 0 ] = i;

                        // create new basis
                        mAllCoarsestBasisOnProc( tCount++ ) = this->create_basis( tIJ,
                                                                                  0,
                                                                                  gNoProcOwner );
                    }
                }

            }
            else if( mNumberOfDimensions == 3)
            {
                // unroll min and max i and j
                mNumberOfCoarsestBasisOnProc[ 0 ] = tIJK( 0, 0 ) + mOrder;

                mNumberOfCoarsestBasisOnProc[ 1 ] = tIJK( 1, 0 ) + mOrder ;

                mNumberOfCoarsestBasisOnProc[ 2 ] = tIJK( 2, 0 ) + mOrder ;

                // initialize array
                mAllCoarsestBasisOnProc.resize( mNumberOfCoarsestBasisOnProc[ 0 ]
                                                *mNumberOfCoarsestBasisOnProc[ 1 ]
                                                *mNumberOfCoarsestBasisOnProc[ 2 ],
                                                nullptr );

                // container with position to be passed to new basis
                luint tIJK[ 3 ];

                // loop over all k
                for( luint k = 0; k < mNumberOfCoarsestBasisOnProc[ 2 ]; ++k )
                {
                    tIJK[ 2 ] = k;
                    // loop over all j
                    for( luint j = 0; j < mNumberOfCoarsestBasisOnProc[ 1 ]; ++j )
                    {
                        // save j-position
                        tIJK[ 1 ] = j;

                        // loop over all i
                        for( luint i = 0; i < mNumberOfCoarsestBasisOnProc[ 0 ]; ++i )
                        {
                            // save i-position
                            tIJK[ 0 ] = i;

                            // create new basis
                            mAllCoarsestBasisOnProc( tCount++ ) = this->create_basis( tIJK,
                                                                                      0,
                                                                                      gNoProcOwner );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::collect_active_and_refined_elements_from_level( const uint              & aLevel,
                                                                                      Cell< Element * > & aElements )
        {
            // cell containing background elements on this level
            Cell< Background_Element_Base* > tBackgroundElements;

            // ask background mesh about elements on this level
            mBackgroundMesh->collect_elements_on_level_including_aura( aLevel,
                                                                       tBackgroundElements );

            // count Elements
            luint tCount = 0;

            for( Background_Element_Base* tBackElement : tBackgroundElements )
            {
                if( ! tBackElement->is_deactive( mActivationPattern ) )
                {
                    ++tCount;
                }
            }

            // allocate cell
            aElements.resize( tCount, nullptr );

            // reset counter
            tCount = 0;
            for( Background_Element_Base* tBackElement : tBackgroundElements )
            {
                if( ! tBackElement->is_deactive( mActivationPattern ) )
                {
                    aElements( tCount++ ) = mAllElementsOnProc( tBackElement->get_memory_index() );
                }
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::process_level( const uint & aLevel )
        {
            Cell< Element* > tElementsOnThisLevel;
            this->collect_active_and_refined_elements_from_level( aLevel, tElementsOnThisLevel );

            Cell< Basis* > tBasisOnThisLevel;

            // collect basis from given level
            this->preprocess_basis_from_level( tElementsOnThisLevel,
                                               tBasisOnThisLevel );

            // determine state of each basis
            this->determine_basis_state( tBasisOnThisLevel );

            // refine B-Spline mesh if this is not the last level
            if ( aLevel < mBackgroundMesh->get_max_level() )
            {
                // initialize basis counter
                luint tCount = 0;

                for( auto tElement : tElementsOnThisLevel )
                {
                    // test if background element has children and is refined on pattern
                    if ( tElement->get_background_element()->has_children() && tElement->is_refined() )
                    {
                        // refine B-Spline element
                        mAllElementsOnProc( tElement->get_memory_index() )->refine( mAllElementsOnProc, tCount );
                    }
                }
            }

            tBasisOnThisLevel.clear();
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::preprocess_basis_from_level( Cell< Element * > & aElements,
                                                             Cell< Basis * >   & aBasis )
        {
            // reset flags for basis
            for( Element * tElement : aElements )
            {
                // loop over all basis from this element
                for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis * tBasis = tElement->get_basis( k );

                    if( tBasis != NULL )
                    {
                        // unflag this basis
                        tBasis->unflag();

                        // unuse this basis
                        tBasis->unuse();

                        // counts this element
                        tBasis->increment_element_counter();
                    }
                }
            }

           // initialize basis counter
           luint tCount = 0;

           // get my rank
           moris_id tMyRank = par_rank();

           for( Element * tElement : aElements )
           {
               // loop over all basis from this element
               for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
               {
                   // get pointer to basis
                   Basis * tBasis = tElement->get_basis( k );

                   if( tBasis != NULL )
                   {
                       // test if basis has been counted
                       if ( ! tBasis->is_flagged() )
                       {
                           // count this basis
                           ++tCount;

                           // flag this basis
                           tBasis->flag();
                       }
                   }
               }

               // use basis if element is owned
               if( tElement->get_owner() == tMyRank )
               {
                   // loop over all basis from this element
                   for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                   {
                       // get pointer to basis
                       Basis * tBasis = tElement->get_basis( k );
                       if( tBasis != NULL )
                       {
                           tBasis->use();
                       }
                   }
               }
           }

           // assign memory for basis container
           aBasis.resize( tCount, nullptr );

           // reset counter
           tCount = 0;

           // loop over all elements
           for( Element * tElement : aElements )
           {
               // loop over all basis from this element
               for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
               {
                   // get pointer to basis
                   Basis * tBasis = tElement->get_basis( k );

                   if( tBasis != NULL )
                   {
                       // test if this basis has been processed already
                       if ( tBasis->is_flagged() )
                       {
                           // copy pointer to basis
                           aBasis( tCount++ ) = tBasis;

                           // unflag this basis
                           tBasis->unflag();
                       }
                   }
               }

               // init neighbor container
               if( tElement->is_refined() )
               {
                   // loop over all basis
                   for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                   {
                       // get pointer to basis
                       Basis * tBasis = tElement->get_basis( k );

                       if( tBasis != NULL )
                       {
                           // assign memory of neighbor container for this basis
                           tBasis->init_neighbor_container();
                       }
                   }
               }
           } // end loop over elements

           // link elements to basis
           for( auto tBasis : aBasis )
           {
               // initialize element container
               tBasis->init_element_container();
           }

           // loop over all elements
           for( Element * tElement : aElements )
           {
               // loop over all basis from this element
               for( uint k = 0; k < mNumberOfBasisPerElement; ++k )
               {
                   Basis * tBasis = tElement->get_basis( k );

                   if( tBasis != NULL )
                   {
                       // insert this element into basis
                       tBasis->insert_element( tElement );
                   }
               }
           }

           // delete_unused_basis (nice feature, not sure if worth the effort)
           //this->delete_unused_basis( aLevel, aBackgroundElements, aBasis );              // FIXME Saves Memory

           // link basis with neighbors
           for( Element * tElement : aElements )
           {
               // calculate basis neighborship
               if ( tElement->is_refined() )
               {
                   // determine basis neighborship
                   tElement->link_basis_with_neighbors( mAllElementsOnProc );
               }
           }

           // reset flag of all basis
           for( auto tBasis : aBasis )
           {
               // initialize element container
               tBasis->unflag();
           }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::delete_unused_basis( const uint                             & aLevel,
                                                           Cell< Background_Element_Base* > & aBackgroundElements,
                                                           Cell< Basis* >                   & aBasis )
        {
            // start timer
            tic tTimer;

            // reset counter for debugging output
            luint tDeleteCount = 0;

            if ( aLevel > 0 )
            {
                // step 1: remove basis from parents

                // collect basis from upper level
                Cell< Basis* > tParents;

                this->collect_basis_from_level( aLevel -1, tParents );

                // loop over all parents
                for( auto tParent : tParents )
                {
                    // loop over all children of this parent
                    for( uint k=0; k<mNumberOfChildrenPerBasis; ++k )
                    {
                        // get pointer to child
                        Basis* tChild = tParent->get_child( k );

                        // test if child exists
                        if( tChild != NULL )
                        {
                            // test if basis is not used
                            if ( ! tChild->is_used() )
                            {
                                // write null into child
                                tParent->insert_child( k, nullptr );
                            }
                        }
                    }
                }

                // tidy up memory
                tParents.clear();

                // step 2: remove basis from basis neighbors
                for( auto tBasis : aBasis )
                {
                    // loop over all neighbors
                    // ( number of neighbors per basis = mNumberOfNeighborsPerElement)
                    for( uint k=0; k<mNumberOfNeighborsPerElement; ++k )
                    {
                       // get pointer to neighbor
                       Basis* tNeighbor = tBasis->get_neighbor( k );

                       // test if neighbor exists
                       if ( tNeighbor != NULL )
                       {
                           // test if neighbor is not used
                           if( ! tNeighbor->is_used() )
                           {
                               // write nullptr into neighbor
                               tBasis->insert_neighbor( k, nullptr );
                           }
                       }
                    }
                }

                // step 3: remove basis from elements
                for( auto tBackElement : aBackgroundElements )
                {
                    // get pointer to element
                    Element* tElement = mAllElementsOnProc( tBackElement->get_memory_index() );

                    // loop over all basis
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

                        // test if basis exists
                        if( tBasis != NULL )
                        {
                            // test if basis is not used
                            if ( ! tBasis->is_used() )
                            {
                                // write null into element
                                tElement->insert_basis( k, nullptr );
                            }
                        }
                    }
                }

                // step 4: count active basis

                // init counter
                luint tCount = 0;
                for( auto tBasis: aBasis )
                {
                    // test if basis is used
                    if ( tBasis->is_used() )
                    {
                        // increment counter
                        ++tCount;
                    }
                }

                // initialize output array
                Cell< Basis* > tBasisOut( tCount, nullptr );

                // reset counter
                tCount = 0;

                for( auto tBasis: aBasis )
                {
                    // test if basis is used
                    if ( tBasis->is_used() )
                    {
                        // increment counter
                        tBasisOut( tCount++ ) = tBasis;
                    }
                    else
                    {
                        delete tBasis;
                        ++tDeleteCount;
                    }
                }

                // remember number of basis for verbose output
                luint tNumberOfAllBasis = aBasis.size();

                // move output
                aBasis.clear();
                aBasis = std::move( tBasisOut );

                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "%s Deleted %lu unused basis of %lu total on level %u.",
                        proc_string().c_str(),
                        ( long unsigned int ) tDeleteCount,
                        ( long unsigned int ) tNumberOfAllBasis,
                        ( unsigned int )      aLevel);

                MORIS_LOG_INFO( "Deleting took %5.3f seconds.",
                        ( double ) tElapsedTime / 1000 );
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::collect_basis_from_level( const uint           & aLevel,
                                                                Cell< Basis* > & aBasis )
        {
            Cell< Element * > tElements;

            this->collect_active_and_refined_elements_from_level( aLevel, tElements );

            // reset flags for basis
            for( Element * tElement : tElements )
            {
                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );

                    // test if basis exists
                    if ( tBasis != NULL )
                    {
                        // unflag this basis
                        tBasis->unflag();
                    }
                }
            }

            // initialize basis counter
            luint tCount = 0;

            // count basis on this level
            for( Element * tElement : tElements )
            {
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );

                    // test if basis exists
                    if ( tBasis != NULL )
                    {
                        // test if basis has been counted
                        if ( ! tBasis->is_flagged() )
                        {
                            // count this basis
                            ++tCount;

                            // flag this basis
                            tBasis->flag();
                        }
                    }
                }
            }

            // assign memory for basis container
            aBasis.resize( tCount, nullptr );

            // reset basis counter
            tCount = 0;

            // add basis to container
            for( Element * tElement : tElements )
            {
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );

                    // test if basis exists
                    if ( tBasis != NULL )
                    {
                        // test if basis has been added
                        if ( tBasis->is_flagged() )
                        {
                            // count this basis
                            aBasis( tCount++ ) = tBasis;

                            // flag this basis
                            tBasis->unflag();
                        }
                    }
                }
            }
        }
//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::determine_basis_state( Cell< Basis* > & aBasis )
        {
            // loop over all basis
            for( Basis * tBasis : aBasis )
            {
                // only process basis that are used by this proc
                if ( tBasis->is_used() )
                {
                    // test number of elements per basis
                    uint tNumberOfElements = tBasis->get_element_counter();

                    // apply deactive lemma
                    if ( tNumberOfElements < mNumberOfElementsPerBasis )
                    {
                        // mark this basis as deactive
                        tBasis->set_deactive_flag();
                    }
                    else
                    {
                        // check if any connected element is deactive
                        bool tHasDeactiveElement = false;

                        for( uint k=0; k<mNumberOfElementsPerBasis; ++k )
                        {
                            if ( tBasis->get_element( k )->is_deactive() )
                            {
                                tHasDeactiveElement = true;
                                break;
                            }
                        }

                        if ( tHasDeactiveElement )
                        {
                            tBasis->set_deactive_flag();
                        }
                        else
                        {
                            bool tIsActive = false;

                            // loop over all basis and check if one is active
                            for( uint k=0; k<mNumberOfElementsPerBasis; ++k )
                            {
                                if ( tBasis->get_element( k )->is_active() )
                                {
                                    tIsActive = true;

                                    // break loop
                                    break;
                                }
                            }

                            if ( tIsActive )
                            {
                                // flag this basis as active
                                tBasis->set_active_flag();
                            }
                            else
                            {
                                // flag this basis as refined
                                tBasis->set_refined_flag();
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::link_basis_to_elements_on_level_zero()
        {
            if ( mNumberOfDimensions == 2 )
            {
                // loop over all elements
                for ( auto tElement : mAllCoarsestElementsOnProc )
                {
                    // init basis container
                    tElement->init_basis_container();

                    // loop over all basis of this element
                    for ( uint k = 0; k < mNumberOfBasisPerElement; ++k )
                    {
                        // get IJ position of this basis
                        luint tIJ[ 2 ];

                        tElement->get_ijk_of_basis( k, tIJ );

                        // insert pointer to basis into element
                        tElement->insert_basis( k, this->get_coarsest_basis_by_ij( tIJ[ 0 ],
                                                                                   tIJ[ 1 ] ) );
                    }
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                {
                    // loop over all elements
                    for ( auto tElement : mAllCoarsestElementsOnProc )
                    {
                        // init basis container
                        tElement->init_basis_container();

                        // loop over all basis of this element
                        for ( uint k = 0; k<mNumberOfBasisPerElement; ++k )
                        {
                            // get IJK position of this basis
                            luint tIJK[ 3 ];
                            tElement->get_ijk_of_basis( k, tIJK );

                            // insert pointer to basis into element
                            tElement->insert_basis( k, this->get_coarsest_basis_by_ijk( tIJK[ 0 ],
                                                                                        tIJK[ 1 ],
                                                                                        tIJK[ 2 ] ) );
                        }
                    }
                }
            }
            else
            {
                MORIS_ERROR( false, "BSpline Mesh: unknown number of dimensions");
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::collect_basis()
        {
            // loop over all coarsest basis
            for( auto tBasis : mAllCoarsestBasisOnProc )
            {
                // collect descendants
                tBasis->flag_descendants();
            }

            // initialize counter
            mNumberOfAllBasis = 0;

            // loop over all coarsest basis
            for( auto tBasis : mAllCoarsestBasisOnProc )
            {
                // collect descendants
                tBasis->count_descendants( mNumberOfAllBasis );
            }

            // assign memory for cell
            mAllBasisOnProc.resize( mNumberOfAllBasis, nullptr );

            // reset counter
            luint tCount = 0;

            // loop over all coarsest basis
            for( auto tBasis : mAllCoarsestBasisOnProc )
            {
                // collect descendants
                tBasis->collect_descendants( mAllBasisOnProc, tCount );
            }

            // reset counter
            tCount = 0;

            // unflag all basis and set index
            for( auto tBasis : mAllBasisOnProc )
            {
                // unflag basis
               tBasis->unflag();

                // set memory index
                tBasis->set_memory_index( tCount++ );
            }


            // This function seems not to be needed, however for debugging purposes the basis can be stored in a vtk file.
            // In this case the information is needed
            // calculate basis coordinates
//            this->calculate_basis_coordinates();

        }
//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::link_basis_to_parents()
        {
            // ask background mesh for max number of levels
            uint tMaxLevel = mBackgroundMesh->get_max_level();

            // Cell containing children
            Cell< Basis* > tChildren;

            // start with level 0

            // loop over all levels but the last
            for( uint l=0; l<tMaxLevel; ++l )
            {
                // make children from last step to parents
                Cell< Basis* > tBasis;
                this->collect_basis_from_level( l, tBasis );

                // loop over all basis on this level
                for( auto tParent: tBasis )
                {
                    // test if parent has children
                    if ( tParent->has_children() )
                    {
                        // loop over all children
                        for( uint k=0; k<mNumberOfChildrenPerBasis; ++k )
                        {
                            Basis * tChild = tParent->get_child( k );

                            // pointer may be null because we deleted unused basis
                            if( tChild != NULL )
                            {
                                // increment parent counter for child
                                tChild->increment_parent_counter();
                            }
                        }
                    }
                }

                // loop over all basis on this level
                for( auto tParent: tBasis )
                {
                    // test if parent has children
                    if ( tParent->has_children() )
                    {
                        // loop over all children
                        for( uint k=0; k<mNumberOfChildrenPerBasis; ++k )
                        {
                            Basis * tChild = tParent->get_child( k );

                            // pointer may be null because we deleted unused basis
                            if( tChild != NULL )
                            {
                                // copy pointer of parent to child
                                tParent->get_child( k )->insert_parent( tParent );
                            }
                        }
                    }
                }
            }
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
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        tElement->get_basis( k )->use();
                    }
                }
            }
        } */

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::calculate_basis_ids()
        {
            if ( mNumberOfDimensions == 2 )
            {
                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // get position of basis
                    const luint * tIJ = tBasis->get_ijk();

                    // calc id and write into basis
                    tBasis->set_domain_id( this->calculate_basis_id( tBasis->get_level(),
                                                                     tIJ[ 0 ],
                                                                     tIJ[ 1 ] ) );
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // get position of basis
                    const luint * tIJK = tBasis->get_ijk();

                    // calc id and write into basis
                    tBasis->set_domain_id( this->calculate_basis_id( tBasis->get_level(),
                                                                     tIJK[ 0 ],
                                                                     tIJK[ 1 ],
                                                                     tIJK[ 2 ] ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::calculate_basis_indices( const Matrix< IdMat > & aCommTable )
        {
            tic tTimer;
            // get my rank
            moris_id tMyRank = par_rank();

            // get number of ranks
            uint tNumberOfProcs = par_size();

            // clean container
            mIndexedBasis.clear();

            // special function for multigrid
            this->flag_refined_basis_of_owned_elements();

            // reset all indices
            for( Basis * tBasis: mAllBasisOnProc )
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
            for ( Basis * tBasis : mAllBasisOnProc )
            {
                tBasis->set_local_index( gNoIndex );
                tBasis->set_domain_index( gNoID );
            }

            // counter for basis
            luint tCount = 0;

            // set local index of basis
            for ( Basis * tBasis : mActiveBasisOnProc )
            {
                if( tBasis->is_flagged() )
                {
                    // set index of basis
                    tBasis->set_local_index( tCount++ );
                }
            }

            if( mParameters->use_multigrid() )
            {
                for ( Basis * tBasis : mRefinedBasisOnProc )
                {
                    if( tBasis->is_flagged() )
                    {
                        // set index of basis
                        tBasis->set_local_index( tCount++ );
                    }
                }
            }

            // allocate container
            mIndexedBasis.resize( tCount, nullptr );

            // reset counter
            tCount = 0;

            // copy indexed basis into container
            for ( Basis * tBasis : mActiveBasisOnProc )
            {
                if( tBasis->is_flagged() )
                {
                    mIndexedBasis( tCount++ ) = tBasis;
                }
            }

            if( mParameters->use_multigrid() )
            {
                for ( Basis * tBasis : mRefinedBasisOnProc )
                {
                    if( tBasis->is_flagged() )
                    {
                        mIndexedBasis( tCount++ ) = tBasis;
                    }
                }
            }

            if ( tNumberOfProcs == 1 )
            {
                // reset counter
                tCount = 0;

                for ( Basis * tBasis : mActiveBasisOnProc )
                {
                    if( tBasis->is_flagged() )
                    {
                        // set index of basis
                        tBasis->set_domain_index( tCount++ );
                    }
                }
               // tNumberOfActiveBasis = tCount;

                if( mParameters->use_multigrid() )
                {
                    for ( Basis * tBasis : mRefinedBasisOnProc )
                    {
                        if( tBasis->is_flagged() )
                        {
                            // set index of basis
                            tBasis->set_domain_index( tCount++ );
                        }
                    }
                }
            } // end serial
            else
            {
                // Step 3: count flagged basis that are owned

                // reset counters
                luint tActiveCount = 0;
                luint tRefinedCount = 0;

                // domain indices (= MTK IDs) loop over all basis
                for ( Basis * tBasis : mActiveBasisOnProc )
                {
                    // test if basis is active, flagged and owned
                    if ( tBasis->get_owner() == tMyRank && tBasis->is_flagged() )
                    {
                        tBasis->set_domain_index( tActiveCount++ );
                    }
                }

                // tNumberOfActiveBasis = tCount;

                if( mParameters->use_multigrid() )
                {
                    for ( Basis * tBasis : mRefinedBasisOnProc )
                    {
                        // test if basis is active, flagged and owned
                        if ( tBasis->get_owner() == tMyRank && tBasis->is_flagged() )
                        {
                            tBasis->set_domain_index( tRefinedCount++ );                          //FIXME should this be active count too?
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

                for( moris_id p=1; p<=tMyRank; ++p )
                {
                    tMyActiveOffset += tActiveBasisCount( p-1 );
                }

                Matrix< DDLUMat > tRefinedBasisCount;
                moris_id tMyRefinedOffset = sum( tActiveBasisCount );

                if( mParameters->use_multigrid() )
                {
                    comm_gather_and_broadcast( tRefinedCount, tRefinedBasisCount );
                    for( moris_id p=1; p<=tMyRank; ++p )
                    {
                        tMyRefinedOffset += tRefinedBasisCount( p-1 );
                    }
                }

                // reset basis counter
                tActiveBasisCount.fill( 0 );

                // loop over all basis
                for( Basis * tBasis : mActiveBasisOnProc )
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

                if( mParameters->use_multigrid() )
                {
                    tRefinedBasisCount.fill( 0 );
                    for( Basis * tBasis : mRefinedBasisOnProc )
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
                Matrix< DDUMat > tProcIndex( tNumberOfProcs, 1, tNumberOfProcs );

                uint tCommLength = aCommTable.length();

                for( uint k=0; k<tCommLength; ++k )
                {
                    tProcIndex( aCommTable( k ) ) = k;
                }

                // - - - - - - - - - - - - - - - -

                // Step 6: allocate memory for communication lists

                // dummy matrces for cells to send
                Matrix< DDLUMat > tEmptyLuint;
                Matrix< DDUMat > tEmptyUint;

                // create cells for basis and element indices to send
                Cell< Matrix< DDLUMat > > tSendIndex   ( tCommLength, tEmptyLuint );
                Cell< Matrix< DDUMat > >  tSendBasis   ( tCommLength, tEmptyUint  );
                Cell< Matrix< DDUMat > >  tSendPedigree( tCommLength, tEmptyUint  );

                // assign memory for Index and Basis
                for( uint p=0; p<tCommLength; ++p )
                {
                    luint tNumberOfBasis;

                    if( mParameters->use_multigrid() )
                    {
                        tNumberOfBasis = tActiveBasisCount( aCommTable( p ) ) + tRefinedBasisCount( aCommTable( p )  );
                    }
                    else
                    {
                        // get number of basis
                        tNumberOfBasis = tActiveBasisCount( aCommTable( p ) );
                    }
                    if ( tNumberOfBasis  > 0 )
                    {
                        ++tCount;
                        tSendIndex( p ).set_size( tNumberOfBasis, 1 );
                        tSendBasis( p ).set_size( tNumberOfBasis, 1 );
                    }
                }

                // - - - - - - - - - - - - - - - -

                // Step 7: create lists with basis of which index is requested

                // reset counter
                tCount = 0;

                // reset counter
                Matrix< DDLUMat > tProcCount( tCommLength, 1, 0 );

                // loop over all basis
                for( Basis* tBasis : mActiveBasisOnProc )
                {
                    // test if basis is active
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is not mine
                        if( tOwner != tMyRank )
                        {
                            // get index of owner
                            uint tIndex = tProcIndex( tOwner );

                            // get counter
                            uint tCount = tProcCount( tIndex );

                            // pointer to element
                            this->get_reference_element_of_basis( tBasis,
                                                                  tSendIndex( tIndex )( tCount ),
                                                                  tSendBasis( tIndex )( tCount ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }

                if( mParameters->use_multigrid() )
                {
                    // loop over all basis
                    for( Basis* tBasis : mRefinedBasisOnProc )
                    {
                        // test if basis is active
                        if ( tBasis->is_flagged() )
                        {
                            // get owner of basis
                            auto tOwner = tBasis->get_owner();

                            // test if basis is not mine
                            if( tOwner != tMyRank )
                            {
                                // get index of owner
                                uint tIndex = tProcIndex( tOwner );

                                // get counter
                                uint tCount = tProcCount( tIndex );

                                // pointer to element
                                this->get_reference_element_of_basis( tBasis,
                                                                      tSendIndex( tIndex )( tCount ),
                                                                      tSendBasis( tIndex )( tCount ) );

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
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tSendIndex( p ).length();

                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        tProcCount( p ) += mAllElementsOnProc( tSendIndex( p )( k ) )
                                                               ->get_background_element()
                                                               ->get_length_of_pedigree_path();
                    }
                }

                // encode pedigree paths
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tSendIndex( p ).length();

                    // assign memory for path to send
                    tSendPedigree( p ).set_size( tProcCount( p ), 1 );

                    // reset counter
                    tCount = 0;

                    // loop over all elements
                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        // get pointer to element
                        Background_Element_Base* tElement = mAllElementsOnProc( tSendIndex( p )( k ) )
                                                                                ->get_background_element();

                        // encode path and overwrite tSendElement with Ancestor Index
                        tElement->endcode_pedigree_path( tSendIndex( p )( k ),
                                                         tSendPedigree( p ),
                                                         tCount );
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
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tReceiveIndex( p ).length();

                    // resize send index
                    tSendIndex( p ).set_size( tNumberOfElements, 1 );

                    // reset counter
                    tCount = 0;

                    // loop over all elements
                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        // decode path and get pointer to element
                        Element* tElement = mAllElementsOnProc( mBackgroundMesh->decode_pedigree_path(
                                                       tReceiveIndex( p )( k ),
                                                       tReceivePedigree( p ),
                                                       tCount )->get_memory_index() );

                        // write index of requested basis into matrix
                        tSendIndex( p )( k )= tElement->get_basis( tReceiveBasis( p )( k ) )
                                                                   ->get_hmr_index();
                    }
                }

                // clear memory
                tReceivePedigree.clear();
                tReceiveBasis   .clear();
                tReceiveIndex   .clear();

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
                for( auto tBasis : mActiveBasisOnProc )
                {
                    // test if basis is flagged
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is mine
                        if( tOwner != tMyRank )
                        {
                            // get index of owner
                            uint tIndex = tProcIndex( tOwner );

                            // get counter
                            tCount = tProcCount( tIndex );

                            // write index into baCommunicationListasis
                            tBasis->set_domain_index( tReceiveIndex( tIndex )( tCount ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }

                if( mParameters->use_multigrid() )
                {
                    for( auto tBasis : mRefinedBasisOnProc )
                    {
                        // test if basis is flagged
                        if ( tBasis->is_flagged() )
                        {
                            // get owner of basis
                            auto tOwner = tBasis->get_owner();

                            // test if basis is mine
                            if( tOwner != tMyRank )
                            {
                                // get index of owner
                                uint tIndex = tProcIndex( tOwner );

                                // get counter
                                tCount = tProcCount( tIndex );

                                // write index into baCommunicationListasis
                                tBasis->set_domain_index( tReceiveIndex( tIndex )( tCount ) );

                                // increment counter
                                ++tProcCount( tIndex );
                            }
                        }
                    }
                }
                // perform a small sanity test :
                tCount = 0;

                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // test if basis is used, active and has no id
                    if (       tBasis->is_flagged()
                            && tBasis->is_active()
                            && tBasis->get_hmr_index() == gNoEntityID )
                    {
                        std::cout << par_rank() << " bad basis " << tBasis->get_hmr_id() << " " << tBasis->get_owner() << std::endl;

                        // increment counter
                        ++tCount;
                    }
                }

                MORIS_ERROR( tCount == 0, "%s ERROR.\n               Could not identify indices of %lu basis.\n               This might happen if a proc uses an active basis that does not belong to\n               itself or any direct neighbor. Suggestion: use denser mesh on top level.\n\n",
                                          proc_string().c_str(), ( long unsigned int ) tCount );
            } // end if parallel

           // insert parents if we are in multigrid
            if( mParameters->use_multigrid() )
            {
                // get parents for each basis
                this->link_basis_to_parents();
            }
            
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( " Calculate basis indices for B-Spline mesh of order %u on pattern %u.",
                    ( unsigned int ) mOrder,
                    ( unsigned int ) mActivationPattern);

            MORIS_LOG_INFO( "Calculation took %5.3f seconds.",
                    ( double ) tElapsedTime / 1000 );
            MORIS_LOG_INFO( " " );
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
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::flag_refined_basis_of_owned_elements()
        {
            if( mParameters->use_multigrid() )
            {
                auto tMyRank = par_rank();

                // loop over all elements
                for( Element * tElement : mAllElementsOnProc )
                {
                    // element must be neither padding or deactive
                    if( tElement->get_owner() == tMyRank )
                    {
                        if( tElement->is_refined() && ! tElement->is_padding() )
                        {
                            // loop over all basis of this element
                            for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                            {
                                // get pointer to basis
                                Basis * tBasis = tElement->get_basis( k );

                                // flag basis if it is refined
                                if( tBasis->is_refined() )
                                {
                                    tBasis->flag();
                                }
                            }
                        }
                    }
                }
            }
        }
//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::synchronize_flags( const Matrix< IdMat > & aCommTable )
        {
            // get number of ranks
            uint tNumberOfProcs = par_size();

            if( tNumberOfProcs > 1 )
            {
                // get my rank
                moris_id tMyRank = par_rank();

                // length of communication table
                uint tCommLength = aCommTable.length();

                // count number of basis per proc
                Matrix< DDUMat > tBasisCount( tNumberOfProcs, 1, 0 );


                // loop over all active basis
                for( auto tBasis : mAllBasisOnProc )
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
                for( uint Ik = 0; Ik < tCommLength; ++Ik )
                {
                    tBasisCommCheck( aCommTable( Ik ) ) = 0;
                }

                // reset my own counter
                tBasisCommCheck( par_rank() ) = 0;

                if( tBasisCommCheck.max() != 0 )
                {
                    std::cout<< "Processor "<< par_rank()<<std::endl;
                    print( aCommTable, "CommTable" );
                    print( tBasisCommCheck, "CommCheck" );
                }

                MORIS_ERROR( tBasisCommCheck.max() == 0, "synchronize_flags: error in communication table" );


                // dummy matrices for cells to send
                Matrix< DDLUMat > tEmptyLuint;
                Matrix< DDUMat >  tEmptyUint;

                // create cells for basis and element indices to send
                Cell< Matrix< DDLUMat > > tSendIndex   ( tCommLength, tEmptyLuint );
                Cell< Matrix< DDUMat > >  tSendBasis   ( tCommLength, tEmptyUint  );
                Cell< Matrix< DDUMat > >  tSendPedigree( tCommLength, tEmptyUint  );

                // reset counter
                luint tCount = 0;

                // assign memory for Index and Basis
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of basis
                    luint tNumberOfBasis = tBasisCount( aCommTable( p ) );

                    if (  tNumberOfBasis  > 0 )
                    {
                        ++tCount;
                        tSendIndex( p ).set_size( tNumberOfBasis, 1 );
                        tSendBasis( p ).set_size( tNumberOfBasis, 1 );
                    }
                }

                // reset counter
                Matrix< DDLUMat > tProcCount( tCommLength, 1, 0 );

                // this table converts the proc id to an index
                Matrix< DDUMat > tProcIndex( tNumberOfProcs, 1, tNumberOfProcs );
                for( uint k=0; k<tCommLength; ++k )
                {
                    tProcIndex( aCommTable( k ) ) = k;
                }

                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // test if basis is active
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is not mine
                        if( tOwner != tMyRank )
                        {
                            // get index of owner
                            uint tIndex = tProcIndex( tOwner );

                            // get counter
                            uint tCount = tProcCount( tIndex );

                            // pointer to element
                            this->get_reference_element_of_basis( tBasis,
                                                                  tSendIndex( tIndex )( tCount ),
                                                                  tSendBasis( tIndex )( tCount ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }

                // local basis IDs received by other procs
                Cell< Matrix< DDUMat > >  tReceiveBasis( tCommLength, tEmptyUint );

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
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tSendIndex( p ).length();

                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        tProcCount( p ) += mAllElementsOnProc(  tSendIndex( p )( k ) )
                                          ->get_background_element()
                                          ->get_length_of_pedigree_path();
                    }
                }

                // encode pedigree paths
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tSendIndex( p ).length();

                    // assign memory for path to send
                    tSendPedigree( p ).set_size( tProcCount( p ), 1 );

                    // reset counter
                    tCount = 0;

                    // loop over all elements
                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        // get pointer to element
                        Background_Element_Base* tElement =  mAllElementsOnProc(  tSendIndex( p )( k ) )
                                                                     ->get_background_element();

                        // encode path and overwrite tSendElement with Ancestor Index
                        tElement->endcode_pedigree_path( tSendIndex( p )( k ),
                                                         tSendPedigree( p ),
                                                         tCount );
                    }
                }

                Cell< Matrix< DDLUMat > > tReceiveIndex   ( tCommLength, tEmptyLuint );
                Cell< Matrix< DDUMat > >  tReceivePedigree( tCommLength, tEmptyUint  );

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
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of elements
                    luint tNumberOfElements = tReceiveIndex( p ).length();

                    // resize send index
                    tSendIndex( p ).set_size( tNumberOfElements, 1 );

                    // reset counter
                    tCount = 0;

                    // loop over all elements
                    for( luint k=0; k<tNumberOfElements; ++k )
                    {
                        // decode path and get pointer to element
                        Element * tElement = mAllElementsOnProc( mBackgroundMesh->decode_pedigree_path(
                                                                     tReceiveIndex( p )( k ),
                                                                     tReceivePedigree( p ),
                                                                     tCount )->get_memory_index() );

                        // now we flag this basis
                        tElement->get_basis( tReceiveBasis( p )( k ) )->flag();
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::collect_active_and_refined_basis()
        {
            // reset counter
            mNumberOfActiveBasisOnProc  = 0;
            mNumberOfRefinedBasisOnProc = 0;
            mNumberOfBasis              = 0;
            mMaxLevel                   = 0;

            if( mParameters->use_multigrid() )
            {
                // count active basis on proc
                for( auto tBasis : mAllBasisOnProc )
                {
                    // reset index
                    tBasis->set_active_index( gNoEntityID );

                    // count basis
                    if( tBasis->is_used() )
                    {
                        ++mNumberOfBasis;
                        if ( tBasis->is_active() )
                        {
                            ++mNumberOfActiveBasisOnProc;
                            mMaxLevel = std::max( tBasis->get_level(), mMaxLevel );
                        }
                        else if( tBasis->is_refined() )
                        {
                            ++mNumberOfRefinedBasisOnProc;
                        }
                    }
                }

                // reserve memory
                mActiveBasisOnProc .resize( mNumberOfActiveBasisOnProc , nullptr );
                mRefinedBasisOnProc.resize( mNumberOfRefinedBasisOnProc, nullptr );

                // reset counters
                mNumberOfActiveBasisOnProc  = 0;
                mNumberOfRefinedBasisOnProc = 0;

                // count active basis on proc
                for( auto tBasis : mAllBasisOnProc )
                {
                    // count basis
                    if( tBasis->is_used() )
                    {
                        ++mNumberOfBasis;
                        if ( tBasis->is_active() )
                        {
                            tBasis->set_active_index( mNumberOfActiveBasisOnProc );                       //FIXME

                            mActiveBasisOnProc( mNumberOfActiveBasisOnProc++ ) = tBasis;
                        }
                        else if( tBasis->is_refined() )
                        {
                            mRefinedBasisOnProc( mNumberOfRefinedBasisOnProc++ ) = tBasis;
                        }
                    }
                }
            }
            else
            {
                // count active basis on proc
                for( auto tBasis : mAllBasisOnProc )
                {
                    // reset index
                    tBasis->set_active_index( gNoEntityID );

                    // count basis
                    if( tBasis->is_used() )
                    {
                        ++mNumberOfBasis;
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
                for( auto tBasis : mAllBasisOnProc )
                {
                    if ( tBasis->is_active() && tBasis->is_used() )
                    {
                        tBasis->set_active_index( mNumberOfActiveBasisOnProc );

                        mActiveBasisOnProc( mNumberOfActiveBasisOnProc++ ) = tBasis;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        bool
        BSpline_Mesh_Base::test_for_double_basis()
        {
            luint tNumberOfBasis = mAllBasisOnProc.size();

            Matrix< DDLUMat > tBasisIDs( tNumberOfBasis, 1 );

            // initialize counter
            luint tCount = 0;

            // populate container
            for( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_used() && ( tBasis->is_active() || tBasis->is_refined() ) )
                    // still have doubled basis if they are not used
                    // by proc however, that does not matter since
                    // they are no DOFs
                {
                    tBasisIDs( tCount++ ) = tBasis->get_hmr_id();
                }
            }

            tBasisIDs.resize( tCount, 1 );

            // make basis unique
            Matrix< DDLUMat > tBasisUniqueIDs;
            unique( tBasisIDs, tBasisUniqueIDs );

            return tBasisUniqueIDs.length() == tCount;
        }

//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::save_to_vtk( const std::string & aFilePath )
        {
            this->calculate_basis_coordinates();

            // start timer
            tic tTimer;

            // modify filename
            std::string tFilePath =  parallelize_path( aFilePath );

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

           for( auto tElement : mAllElementsOnProc )
           {
               if( tElement->get_owner() == tMyRank )
               {
                   // flag all basis of this element
                   for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                   {
                       tElement->get_basis( k )->flag();
                   }
               }
           }

           // count flagged basis
           uint tNumberOfBasis = 0;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   // increment counter
                   ++tNumberOfBasis;
               }
           }

           // write node data
           tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
           tFile << "POINTS " << tNumberOfBasis << " float"  << std::endl;

           // ask settings for numner of dimensions
           auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

           if ( tNumberOfDimensions == 2 )
           {
               for( auto tBasis : mAllBasisOnProc )
               {
                   if( tBasis->is_flagged() )
                   {
                       // get coordinate from basis
                       const real* tXY = tBasis->get_xyz();

                       // write coordinates to mesh
                       tFChar = swap_byte_endian( (float) tXY[ 0 ] );
                       tFile.write( (char*) &tFChar, sizeof(float));
                       tFChar = swap_byte_endian( (float) tXY[ 1 ] );
                       tFile.write( (char*) &tFChar, sizeof(float));
                       tFChar = swap_byte_endian( (float) 0 );
                       //tFChar = swap_byte_endian( (float) tBasis->get_level() );
                       tFile.write( (char*) &tFChar, sizeof(float));
                   }
               }
           }
           else if ( tNumberOfDimensions == 3 )
           {
               for( auto tBasis : mAllBasisOnProc )
               {
                   if( tBasis->is_flagged() )
                   {
                       // get coordinate from node
                       const real* tXYZ = tBasis->get_xyz();

                       // write coordinates to mesh
                       tFChar = swap_byte_endian( (float) tXYZ[ 0 ] );
                       tFile.write( (char*) &tFChar, sizeof(float));
                       tFChar = swap_byte_endian( (float) tXYZ[ 1 ] );
                       tFile.write( (char*) &tFChar, sizeof(float));
                       tFChar = swap_byte_endian( (float) tXYZ[ 2 ] );
                       tFile.write( (char*) &tFChar, sizeof(float));
                   }
               }
           }

           tFile << std::endl;

           // write each basis as its own element
           tFile << "CELLS " << tNumberOfBasis << " "
                   << 2*tNumberOfBasis << std::endl;

           int tOne = swap_byte_endian( (int) 1 );

           // reset counter
           int tCount = 0;

           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( tCount );
                   tFile.write( ( char* ) &tOne, sizeof(int));
                   tFile.write( ( char *) &tIChar, sizeof(int));

                   ++tCount;
               }
           }

           // write cell types
           tFile << "CELL_TYPES " << tNumberOfBasis << std::endl;
           tIChar = swap_byte_endian( (int) 2 );
           for ( luint k = 0; k < tNumberOfBasis; ++k)
           {
               tFile.write( (char*) &tIChar, sizeof( int ) );
           }

           // write node data
           tFile << "POINT_DATA " << tNumberOfBasis << std::endl;

           // write state
           tFile << "SCALARS STATE int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
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

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write basis ID
           tFile << "SCALARS ID int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( (int) tBasis->get_id() );

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write internal basis ID
           tFile << "SCALARS HMR_ID int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( (int) tBasis->get_hmr_id() );

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write basis level
           tFile << "SCALARS LEVEL int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( (int) tBasis->get_level() );

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write basis owner
           tFile << "SCALARS OWNER int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( (int) tBasis->get_owner() );

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write memory index
           tFile << "SCALARS MEMORY_INDEX int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   tIChar = swap_byte_endian( (int) tBasis->get_memory_index() );

                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }
           tFile << std::endl;

           // write active index
           tFile << "SCALARS ACTIVE_INDEX int" << std::endl;
           tFile << "LOOKUP_TABLE default" << std::endl;
           for( auto tBasis : mAllBasisOnProc )
           {
               if( tBasis->is_flagged() )
               {
                   if( tBasis->is_active() )
                   {
                       tIChar = swap_byte_endian( (int) tBasis->get_active_index() );
                   }
                   else
                   {
                       tIChar = swap_byte_endian( (int) -1 );
                   }
                   tFile.write( ( char *) &tIChar, sizeof(int));
               }
           }

           tFile << std::endl;

           // close the output file
           tFile.close();

           // unflag all bases
           this->unflag_all_basis();

           // stop timer
           real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

           // print output
           MORIS_LOG_INFO( "%s Created VTK debug file.",
                   proc_string().c_str() );

           MORIS_LOG_INFO( "Mesh has %lu basis.",
                   ( long unsigned int ) tNumberOfBasis );

           MORIS_LOG_INFO( "Creation took %5.3f seconds.",
                   ( double ) tElapsedTime / 1000 );
           MORIS_LOG_INFO( " " );
        }
//------------------------------------------------------------------------------

        void BSpline_Mesh_Base::calculate_child_stencil()
        {
            // number of children in nd
            uint tNumberOfChildren = this->get_number_of_children_per_basis();

            // get order
            uint tOrder = this->get_order();

            uint tNumberOfChildrenPerDirection = tOrder + 2;

            // allocate matrix
            mChildStencil.set_size( tNumberOfChildren, 1 );

            uint tChild = 0;

            switch ( mParameters->get_number_of_dimensions() )
            {
                case( 2 ) :
                {
                    for ( uint j=0; j< tNumberOfChildrenPerDirection ; ++j )
                    {
                        for( uint i=0; i<tNumberOfChildrenPerDirection; ++i )
                        {
                            mChildStencil( tChild++ ) =  nchoosek( tOrder+1, i ) * nchoosek( tOrder+1, j ) / std::pow( 2, tOrder ) / std::pow( 2, tOrder );
                        }
                    }
                    break;
                }
                case( 3 ) :
                {
                    for( uint k=0; k< tNumberOfChildrenPerDirection; ++k )
                    {
                        for ( uint j=0; j< tNumberOfChildrenPerDirection ; ++j )
                        {
                            for( uint i=0; i<tNumberOfChildrenPerDirection; ++i )
                            {
                                mChildStencil( tChild++ ) = nchoosek( tOrder+1, i ) * nchoosek( tOrder+1, j ) * nchoosek( tOrder+1, k ) / std::pow( 2, tOrder )
                                                                                                                                        / std::pow( 2, tOrder )
                                                                                                                                        / std::pow( 2, tOrder );
                            }
                        }
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Ivalid dimension.");
                    break;
                }
            }

            //print(mChildStencil,"mChildStencil");
            //mChildStencil = mChildStencil / std::pow( 2, tOrder );
        }

//------------------------------------------------------------------------------

        Matrix< DDSMat > BSpline_Mesh_Base::get_children_ind_for_basis( const moris::sint aParentBasind )
        {
            // get basis pointer
            Basis * tBasis = this->get_basis_by_index( aParentBasind );

            // get child indices
            Matrix< DDSMat > tBasisLocalChildInds;

            tBasis->get_basis_local_child_inds( tBasisLocalChildInds  );

            uint tNumberOfChildren = tBasisLocalChildInds.length();

            Matrix< DDSMat > tChildInds( tNumberOfChildren, 1 );

            for( uint k=0; k<tNumberOfChildren; ++k )
            {
                tChildInds( k ) = tBasis->get_child( tBasisLocalChildInds( k ) )->get_index();
            }

            return tChildInds;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > BSpline_Mesh_Base::get_children_weights_for_parent( const moris::sint aParentBasind )
        {
            // get basis pointer
            Basis * tBasis = this->get_basis_by_index( aParentBasind );

            // get child indices
            Matrix< DDSMat > tBasisLocalChildInds;

            tBasis->get_basis_local_child_inds( tBasisLocalChildInds  );

            uint tNumberOfChildren = tBasisLocalChildInds.length();

            // create weights
            Matrix< DDRMat > tWeights( tNumberOfChildren, 1 );

            for( uint k=0; k<tNumberOfChildren; ++k )
            {
                tWeights( k ) = mChildStencil ( tBasisLocalChildInds( k ) );
            }
            return tWeights;
        }

//------------------------------------------------------------------------------

        uint BSpline_Mesh_Base::get_number_of_basis_connected_to_basis( const moris_index aIndex )
        {
            // get basis pointer
            Basis * tBasis = mActiveBasisOnProc( aIndex );

            // step 1: unflag all connected basis from connected elements

            // get number of connected elements
            uint tNumberOfElements = tBasis->get_element_counter();

            // loop over all elements of this basis
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get Cell of connected vertices
                moris::Cell< mtk::Vertex* > tVertices = tBasis->get_element( e )->get_vertex_pointers();

                // unflag these vertices
                for( mtk::Vertex* tVertex : tVertices )
                {
                    tVertex->unflag();
                }
            }

            // step 2: count connected vertices

            // reset counter
            uint tCount = 0;

            // loop over all elements of this basis
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get Cell of connected vertices
                moris::Cell< mtk::Vertex* > tVertices = tBasis->get_element( e )->get_vertex_pointers();

                // unflag these vertices
                for( mtk::Vertex* tVertex : tVertices )
                {
                    if ( ! tVertex->is_flagged() )
                    {
                        tVertex->unflag();

                        ++tCount;
                    }
                }
            }

            // return counter
            return tCount;
        }

    } /* namespace hmr */
} /* namespace moris */
