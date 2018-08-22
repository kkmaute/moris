/*
 * cl_HMR_BSpline_Mesh_Base.cpp
 *
 *  Created on: Jun 12, 2018
 *      Author: messe
 */
#include <fstream>
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_Mat.hpp" //LNA/src
#include "fn_unique.hpp" //LNA/src
#include "HMR_Tools.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------
//   public:
//------------------------------------------------------------------------------

        BSpline_Mesh_Base::BSpline_Mesh_Base (
                const Parameters       * aParameters,
                Background_Mesh_Base * aBackgroundMesh,
                const uint           & aOrder ) :
                Mesh_Base(
                    aParameters,
                    aBackgroundMesh,
                    aOrder ),
                    mNumberOfChildrenPerBasis(
                            std::pow( aOrder + 2,
                            aParameters->get_number_of_dimensions() ) ),
                    mNumberOfElementsPerBasis( std::pow( aOrder+1 ,
                            aParameters->get_number_of_dimensions() ) )

        {
            this->calculate_basis_level_offset();
        }

//------------------------------------------------------------------------------


        void
        BSpline_Mesh_Base::update_mesh()
        {
            // start timer
            tic tTimer;

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
            this->collect_active_basis();

            // determine indices of active and flagged basis
            //this->calculate_basis_indices();

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created B-Spline mesh.\n               Mesh has %lu  Elements and %lu basis in total.\n               Mesh uses %lu basis on proc.\n               Mesh has %lu active basis on proc.\n               Creation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) mNumberOfAllElementsOnProc,
                        ( long unsigned int ) mNumberOfAllBasis,
                        ( long unsigned int ) mNumberOfBasis,
                        ( long unsigned int ) mNumberOfActiveBasisOnProc,
                        ( double ) tElapsedTime / 1000 );
            }

        }

//------------------------------------------------------------------------------

        Basis*
        BSpline_Mesh_Base::get_coarsest_basis_by_ij(
                            const luint & aI,
                            const luint & aJ )
        {
            if ( aI < mNumberOfCoarsestBasisOnProc[ 0 ] &&
                 aJ < mNumberOfCoarsestBasisOnProc[ 1 ] )
            {
                return mAllCoarsestBasisOnProc(
                        aI + aJ*mNumberOfCoarsestBasisOnProc[ 0 ] );
            }
            else
            {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------

        Basis*
        BSpline_Mesh_Base::get_coarsest_basis_by_ijk(
                const luint & aI,
                const luint & aJ,
                const luint & aK )
        {
            if (    aI < mNumberOfCoarsestBasisOnProc[ 0 ] &&
                    aJ < mNumberOfCoarsestBasisOnProc[ 1 ] &&
                    aK < mNumberOfCoarsestBasisOnProc[ 2 ] )
            {
                return mAllCoarsestBasisOnProc(
                        aI + mNumberOfCoarsestBasisOnProc[ 0 ] *
                        ( aJ + aK*mNumberOfCoarsestBasisOnProc[ 1 ] ) );
            }
            else
            {
                return nullptr;
            }

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
        bool
        BSpline_Mesh_Base::test_sanity()
        {

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
                        tTestTopLevelState = tTestTopLevelState
                                && ( tBasis->is_active() || tBasis->is_refined() );

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

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                if ( aPassedTest )
                {
                    // print output
                    std::fprintf( stdout,"%s Tested basis activation sanity.\n               Test took %5.3f seconds.\n               All tests passed.\n\n",
                        proc_string().c_str(),
                        ( double ) tElapsedTime / 1000 );
                }
                else
                {
                        // print output
                        std::fprintf( stdout,"%s Tested basis activation sanity.\n               Test took %5.3f seconds.\n               AT LEAST ONE TEST FAILED.\n\n",
                        proc_string().c_str(),
                                ( double ) tElapsedTime / 1000 );

                        std::cout << "Test result: "
                                  <<  tTestForStateContratiction << " "
                                  << tTestTopLevelState << " "
                                  << tHaveRefinedParent << " "
                                  << tDeactiveTest << " "
                                  << tRefinedHasActiveChild  << std::endl;
                }

            }

            return aPassedTest;
        }

//------------------------------------------------------------------------------
// private:
//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::calculate_basis_level_offset()
        {
            // first level is zero
            mBasisLevelOffset[ 0 ] = 0;

            // get number of elements from Background Mesh
            Mat< luint > tNumberOfElements
                = mBackgroundMesh->get_number_of_elements_per_direction();


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
                    tBasisOnLastLevel *= tBasisPerDirection[ k ];
                    tBasisPerDirection[ k ] *= 2;
                    tBasisPerDirection[ k ] += mOrder;
                }

                // increment delta

                // add offset to table
                mBasisLevelOffset[ l+1 ] =  mBasisLevelOffset[ l ]
                                            + tBasisOnLastLevel;
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::create_basis()
        {
            // basis on first level are created separately
            this->create_basis_on_level_zero();

            // connect basis to coarsest elements
            this->link_basis_to_elements_on_level_zero();

            // ask background mesh for number of levels
            luint tNumberOfLevels = mBackgroundMesh->get_max_level();
            for( uint l=0; l<=tNumberOfLevels; ++l )
            {
                this->process_level( l );
            }

            this->collect_basis();

            //this->calculate_basis_coordinates();
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::create_basis_on_level_zero()
        {
            // ask mesh for relevant ijk positions
            Mat< luint > tIJK =
                    mBackgroundMesh
                    ->get_number_of_elements_per_direction_on_proc();

            // initialize basis counter
            luint tCount = 0;

            if( mNumberOfDimensions == 2)
            {
                // unroll min and max i and j
                mNumberOfCoarsestBasisOnProc[ 0 ]
                    = tIJK( 0, 0 ) + mOrder;

                mNumberOfCoarsestBasisOnProc[ 1 ]
                    = tIJK( 1, 0 ) + mOrder ;

                // initialize array
                mAllCoarsestBasisOnProc.resize(
                        mNumberOfCoarsestBasisOnProc[ 0 ]
                       *mNumberOfCoarsestBasisOnProc[ 1 ],
                       nullptr );



                // container with position to be passed to new basis
                luint tIJ[ 2 ];

                // loop over all j
                for( luint j=0; j<mNumberOfCoarsestBasisOnProc[ 1 ]; ++j )
                {
                    // save j-position
                    tIJ[ 1 ] = j;

                    // loop over all i
                    for( luint i=0; i<mNumberOfCoarsestBasisOnProc[ 0 ]; ++i )
                    {
                        // save i-position
                        tIJ[ 0 ] = i;

                        // create new basis
                        mAllCoarsestBasisOnProc( tCount++ )
                            = this->create_basis( tIJ, 0, gNoProcOwner );
                    }
                }

            }
            else if( mNumberOfDimensions == 3)
            {
                // unroll min and max i and j
                mNumberOfCoarsestBasisOnProc[ 0 ]
                    = tIJK( 0, 0 ) + mOrder;

                mNumberOfCoarsestBasisOnProc[ 1 ]
                    = tIJK( 1, 0 ) + mOrder ;

                mNumberOfCoarsestBasisOnProc[ 2 ]
                    = tIJK( 2, 0 ) + mOrder ;

                // initialize array
                mAllCoarsestBasisOnProc.resize(
                        mNumberOfCoarsestBasisOnProc[ 0 ]
                       *mNumberOfCoarsestBasisOnProc[ 1 ]
                       *mNumberOfCoarsestBasisOnProc[ 2 ],
                       nullptr );

                // container with position to be passed to new basis
                luint tIJK[ 3 ];

                // loop over all k
                for( luint k=0; k<mNumberOfCoarsestBasisOnProc[ 2 ]; ++k )
                {
                    tIJK[ 2 ] = k;
                    // loop over all j
                    for( luint j=0; j<mNumberOfCoarsestBasisOnProc[ 1 ]; ++j )
                    {
                        // save j-position
                        tIJK[ 1 ] = j;

                        // loop over all i
                        for( luint i=0; i<mNumberOfCoarsestBasisOnProc[ 0 ]; ++i )
                        {
                            // save i-position
                            tIJK[ 0 ] = i;

                            // create new basis
                            mAllCoarsestBasisOnProc( tCount++ )
                                = this->create_basis( tIJK, 0,  gNoProcOwner );
                        }
                    }
                }

            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::process_level( const uint& aLevel )
        {
            // cell containing background elements on this level
            Cell< Background_Element_Base* > tBackgroundElements;

            // ask background mesh about elements on this level
            mBackgroundMesh->collect_elements_on_level_including_aura(
                    aLevel,
                    tBackgroundElements );

            // cell containing basis on level
            Cell< Basis* > tBasisOnThisLevel;

            // collect basis from given level
            this->preprocess_basis_from_level(
                    aLevel,
                    tBackgroundElements,
                    tBasisOnThisLevel );


            // determine state of each basis
            this->determine_basis_state( tBasisOnThisLevel );

            // refine B-Spline mesh if this is not the last level
            if ( aLevel < mBackgroundMesh->get_max_level() )
            {
                // initialize basis counter
                luint tCount = 0;

                for( auto tElement : tBackgroundElements )
                {


                    // test if backgroud element has children
                    if ( tElement->has_children() )
                    {
                        // refine B-Spline element
                        mAllElementsOnProc( tElement->get_memory_index() )
                                        ->refine( mAllElementsOnProc, tCount );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::preprocess_basis_from_level(
                const uint                       & aLevel,
                Cell< Background_Element_Base* > & aBackgroundElements,
                Cell< Basis* >                   & aBasis)
        {

            // reset flags for basis
            for ( auto tBackElement : aBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );

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
            uint tMyRank = par_rank();

            // count basis on this level
            for ( auto tBackElement : aBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );



                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );
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

                // get owner of element
                auto tOwner = tElement->get_owner();

                // use basis if element is owned
                if( tOwner == tMyRank )
                {
                    // loop over all basis from this element
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

                        tBasis->use();
                    }
                }
            }

            // assign memory for basis container
            aBasis.resize( tCount, nullptr );

            // reset basis counter
            tCount = 0;

            // loop over all elements
            for ( auto tBackElement : aBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );
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

                // init neighbor container if element is refined
                if( tBackElement->is_refined() )
                {
                    // loop over all basis
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to basis
                        Basis* tBasis = tElement->get_basis( k );

                        if( tBasis != NULL )
                        {
                            // assign memory of neighbor container for this basis
                            tBasis->init_neighbor_container();
                        }
                    }
                }
            }

            // link elements to basis
            for( auto tBasis : aBasis )
            {
                // initialize element container
                tBasis->init_element_container();
            }

            // loop over all elements
            for ( auto tBackElement : aBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    Basis* tBasis = tElement->get_basis( k );

                    if( tBasis != NULL )
                    {
                        // insert this element into basis
                        tBasis->insert_element( tElement );
                    }

                }

            }

            // delete_unused_basis (nice feature, not sure if worth the effort)
            this->delete_unused_basis( aLevel, aBackgroundElements, aBasis );

            // link basis with neighbors
            for ( auto tBackElement : aBackgroundElements )
            {


                // calculate basis neighborship
                if ( tBackElement->is_refined() )
                {
                    // get pointer to element
                    Element* tElement = mAllElementsOnProc(
                            tBackElement->get_memory_index() );

                    // determine basis neighborship
                    tElement->link_basis_with_neighbors( mAllElementsOnProc );
                }

            }

            // rest flag of all basis
            for( auto tBasis : aBasis )
            {
                // initialize element container
                tBasis->unflag();
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::delete_unused_basis(
                const uint                       & aLevel,
                Cell< Background_Element_Base* > & aBackgroundElements,
                Cell< Basis* >                   & aBasis )
        {
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
                    Element* tElement = mAllElementsOnProc(
                            tBackElement->get_memory_index() );

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
                    }
                }

                // move output
                aBasis.clear();
                aBasis = std::move( tBasisOut );
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::collect_basis_from_level(
                const uint     & aLevel,
                Cell< Basis* > & aBasis )
        {

            Cell< Background_Element_Base* > tBackgroundElements;

            // get elements from background mesh
            mBackgroundMesh
                ->collect_elements_on_level_including_aura(
                        aLevel,
                        tBackgroundElements );

            // reset flags for basis
            for ( auto tBackElement : tBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

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
            for ( auto tBackElement : tBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

                // loop over all basis from this element
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

            // loop over all elements
            for ( auto tBackElement : tBackgroundElements )
            {
                // get pointer to element
                Element* tElement = mAllElementsOnProc(
                        tBackElement->get_memory_index() );

                // loop over all basis from this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    // get pointer to basis
                    Basis* tBasis = tElement->get_basis( k );
                    // test if basis exists
                    if ( tBasis != NULL )
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
            }
        }
//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::determine_basis_state( Cell< Basis* > & aBasis )
        {
            // loop over all basis
            for( auto tBasis : aBasis )
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
                        // flag this basis as refined
                        tBasis->set_refined_flag();

                        // loop over all basis and check if one is active
                        for( uint k=0; k<mNumberOfElementsPerBasis; ++k )
                        {
                            if ( tBasis->get_element( k )->is_active() )
                            {
                                // flag this basis as active
                                tBasis->set_active_flag();

                                // break loop
                                break;
                            }
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::link_basis_to_elements_on_level_zero()
        {
            if ( mNumberOfDimensions == 2 )
            {
                // loop over all elements
                for ( auto tElement : mAllCoarsestElementsOnProc )
                {
                    // loop over all basis of this element
                    for ( uint k = 0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get IJ position of this basis
                        luint tIJ[ 2 ];
                        tElement->get_ijk_of_basis( k, tIJ );

                        // insert pointer to basis into element
                        tElement->insert_basis( k,
                                this->get_coarsest_basis_by_ij(
                                        tIJ[ 0 ],
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
                        // loop over all basis of this element
                        for ( uint k = 0; k<mNumberOfBasisPerElement; ++k )
                        {
                            // get IJK position of this basis
                            luint tIJK[ 3 ];
                            tElement->get_ijk_of_basis( k, tIJK );

                            // insert pointer to basis into element
                            tElement->insert_basis( k,
                                    this->get_coarsest_basis_by_ijk(
                                            tIJK[ 0 ],
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

        void
        BSpline_Mesh_Base::collect_basis()
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

            // calculate basis coordinates
            this->calculate_basis_coordinates();


        }
//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::link_basis_to_parents()
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
                            // increment parent counter for child
                            tParent->get_child( k )->increment_parent_counter();
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
                            // copy pointer of parent to child
                            tParent->get_child( k )->insert_parent( tParent );
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

        void
        BSpline_Mesh_Base::calculate_basis_ids()
        {
            if ( mNumberOfDimensions == 2 )
            {
                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // get position of basis
                    const luint * tIJ = tBasis->get_ijk();

                    // calc id and write into basis
                    tBasis->set_domain_id(
                            this->calculate_basis_id(
                                    tBasis->get_level(),
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
                    tBasis->set_domain_id(
                            this->calculate_basis_id(
                                    tBasis->get_level(),
                                    tIJK[ 0 ],
                                    tIJK[ 1 ],
                                    tIJK[ 2 ] ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::calculate_basis_indices()
        {
            // get my rank
            uint tMyRank = par_rank();

            // get number of ranks
            uint tNumberOfProcs = par_size();

            // counter for basis
            luint tCount = 0;

            if ( tNumberOfProcs == 1 )
            {
                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // test of basis is active and flagged
                    if ( tBasis->is_active() && tBasis->is_flagged() )
                    {
                        // set index of basis
                        tBasis->set_domain_index( tCount++ );
                    }
                    else
                    {
                        tBasis->set_domain_index( gNoEntityID );
                    }
                }
            }
            else
            {
                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // test if basis is active, flagged and owned
                    if ( tBasis->is_active() && tBasis->get_owner() == tMyRank && tBasis->is_flagged() )
                    {
                        tBasis->set_domain_index( tCount++ );
                    }
                    else
                    {
                        tBasis->set_domain_index( gNoEntityID );
                    }
                }

                // communicate number of owned and activa basis with other procs
                Mat< luint > tBasisCount = comm_gather_and_broadcast( tCount );

                // get my offset
                luint tMyOffset = 0;

                for( uint p=1; p<=tMyRank; ++p )
                {
                    tMyOffset += tBasisCount( p-1 );
                }

                // reset basis counter
                tBasisCount.fill( 0 );

                // loop over all basis
                for( auto tBasis : mActiveBasisOnProc )
                {
                    // test if basis is flagged
                    if ( tBasis->is_flagged() )
                    {
                        // get owner of basis
                        auto tOwner = tBasis->get_owner();

                        // test if basis is mine
                        if ( tOwner == tMyRank )
                        {
                            tBasis->set_domain_index(
                                    tBasis->get_domain_index()
                                    + tMyOffset );
                        }
                        else
                        {
                            // increment basis counter per proc
                            ++tBasisCount( tOwner );
                        }

                    }
                }

                // tell each proc how many basis are requested :
                // communication list for all procs
                Mat< uint > tAllProcs( tNumberOfProcs, 1 );
                for( uint k = 0; k<tNumberOfProcs; ++k )
                {
                    tAllProcs( k ) = k;
                }

                // list of basis requested by other procs
                Mat< luint > tBasisRequest;

                communicate_scalars(
                        tAllProcs,
                        tBasisCount,
                        tBasisRequest );

                // create communication list
                Mat< uint > tCommunicationList( tNumberOfProcs, 1, 0 );

                // index for each proc
                Mat< uint > tProcIndex( tNumberOfProcs, 1, tNumberOfProcs );

                // reset counter
                tCount = 0;

                // loop over all procs and check if a communication is necessary
                for( uint p=0; p<tNumberOfProcs; ++p )
                {
                    if ( tBasisCount( p ) != 0 || tBasisRequest( p ) != 0 )
                    {
                        tProcIndex( p ) = tCount;
                        tCommunicationList( tCount++ ) = p;
                    }
                }

                // remember length of communication list
                uint tCommLength = tCount;

                // rescale communication list
                tCommunicationList.resize( tCommLength, 1 );

                // dummy matrces for cells to send
                Mat< luint > tEmptyLuint;
                Mat< uint > tEmptyUint;

                // create cells for basis and element indices to send
                Cell< Mat< luint > > tSendIndex( tCommLength, tEmptyLuint );
                Cell< Mat<  uint > > tSendBasis( tCommLength, tEmptyUint );
                Cell< Mat<  uint > > tSendPedigree( tCommLength, tEmptyUint );

                // reset counter
                tCount = 0;

                // assign memory for Index and Basis
                for( uint p=0; p<tCommLength; ++p )
                {
                    // get number of basis
                    luint tNumberOfBasis = tBasisCount( tCommunicationList( p ) );

                    if (  tNumberOfBasis  > 0 )
                    {
                        ++tCount;
                        tSendIndex( p ).set_size( tNumberOfBasis, 1 );
                        tSendBasis( p ).set_size( tNumberOfBasis, 1 );
                    }
                }

                // reset counter
                Mat< luint > tProcCount( tCommLength, 1, 0 );

                // loop over all basis
                for( auto tBasis : mActiveBasisOnProc )
                {
                    // test if basis is active
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
                            uint tCount = tProcCount( tIndex );

                            // pointer to element

                            this->get_reference_element_of_basis(
                                    tBasis,
                                    tSendIndex( tIndex )( tCount ),
                                    tSendBasis( tIndex )( tCount ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }

                // local basis IDs received by other procs
                Cell< Mat< uint > >  tReceiveBasis( tCommLength, tEmptyUint );

                // communicate local basis IDs to request
                communicate_mats(
                        tCommunicationList,
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
                        tProcCount( p ) +=
                                mAllElementsOnProc(  tSendIndex( p )( k ) )
                                ->get_background_element()->
                                get_length_of_pedigree_path();

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
                        Background_Element_Base* tElement
                        =  mAllElementsOnProc(  tSendIndex( p )( k ) )
                        ->get_background_element();

                        // encode path and overwrite tSendElement with Ancestor Index
                        tElement->endcode_pedigree_path(
                                tSendIndex( p )( k ),
                                tSendPedigree( p ),
                                tCount );
                    }
                }

                Cell< Mat< luint > > tReceiveIndex( tCommLength, tEmptyLuint );
                Cell< Mat< uint > >  tReceivePedigree( tCommLength, tEmptyUint );

                // communicate ancestor IDs
                communicate_mats(
                        tCommunicationList,
                        tSendIndex,
                        tReceiveIndex );

                // communicate pedigree paths
                communicate_mats(
                        tCommunicationList,
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
                        Element*
                        tElement =
                                mAllElementsOnProc(
                                        mBackgroundMesh->decode_pedigree_path(
                                                tReceiveIndex( p )( k ),
                                                tReceivePedigree( p ),
                                                tCount )->get_memory_index() );

                        // write index of requested basis into matrix
                        tSendIndex( p )( k )
                                               = tElement->get_basis( tReceiveBasis( p )( k ) )
                                               ->get_domain_index();
                    }
                }

                // clear memory
                tReceivePedigree.clear();
                tReceiveBasis.clear();
                tReceiveIndex.clear();

                // communicate requested indices back to original proc
                communicate_mats(
                        tCommunicationList,
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

                            // write index into basis
                            tBasis->set_domain_index(
                                    tReceiveIndex( tIndex )( tCount ) );

                            // increment counter
                            ++tProcCount( tIndex );
                        }
                    }
                }

                // perform a small sanity test :
                tCount = 0;

                // loop over all basis
                for( auto tBasis : mAllBasisOnProc )
                {
                    // test if basis is used, active and has no id
                    if (       tBasis->is_used()
                            && tBasis->is_active()
                            && tBasis->get_domain_index() == gNoEntityID )
                    {
                        // increment counter
                        ++tCount;
                    }
                }

                if ( tCount != 0 )
                {
                    /**
                     * This error might occur if a proc uses an active basis
                     * that belongs to a proc that is not a direct neighbor.
                     */
                    std::fprintf( stdout,"%s ERROR.\n               Could not identify indices of %lu basis.\n               This might happen if a proc uses an active basis that does not belong to\n               itself or any direct neighbor. Suggestion: use denser mesh on top level.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tCount );
                    exit( -1 );
                }
            }

        }


//------------------------------------------------------------------------------

        void
        BSpline_Mesh_Base::collect_active_basis()
        {
            // reset counter
            mNumberOfActiveBasisOnProc = 0;

            mNumberOfBasis = 0;

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
                    }
                }
            }

            // reserve memory
            mActiveBasisOnProc.resize( mNumberOfActiveBasisOnProc, nullptr );

            // initialize counter
            luint tCount = 0;

            // populate container
            for( auto tBasis : mAllBasisOnProc )
            {
                if ( tBasis->is_active() && tBasis->is_used() )
                {
                    tBasis->set_active_index( tCount );

                    mActiveBasisOnProc( tCount++ ) = tBasis;
                }
            }
        }

//------------------------------------------------------------------------------
        bool
        BSpline_Mesh_Base::test_for_double_basis()
        {
            luint tNumberOfBasis = mAllBasisOnProc.size();

            Mat< luint > tBasisIDs( tNumberOfBasis, 1 );

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

                    tBasisIDs( tCount++ ) = tBasis->get_domain_id();

                }
            }

            tBasisIDs.resize( tCount, 1 );

            // make basis unique
            tBasisIDs = unique( tBasisIDs );

            return tBasisIDs.length() == tCount;
        }

//------------------------------------------------------------------------------
        void
        BSpline_Mesh_Base::save_to_vtk( const std::string & aFilePath )
        {
            // start timer
            tic tTimer;

            // modify filename
            std::string tFilePath;
            if ( moris::par_size() > 1 )
            {
                auto tFileExt = aFilePath.substr(aFilePath.find_last_of("."),aFilePath.length());
                auto tBasePath = aFilePath.substr(0,aFilePath.find_last_of("."));
                tFilePath = tBasePath + "_" +  std::to_string(moris::par_rank()) + tFileExt;

            }
            else
            {
                tFilePath = aFilePath;
            }

            // open the file
            std::ofstream tFile( tFilePath, std::ios::binary );

            // containers
            float tFChar = 0;
            int   tIChar = 0;

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;

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
                    // tFChar = swap_byte_endian( (float) 0 );
                    tFChar = swap_byte_endian( (float) mAllBasisOnProc( k )->get_level() );
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



            // count number of non padding elements
            luint tNumberOfElements = 0;

            // can only write element data if vtk map exists

            int tNumberOfNodesPerElement = swap_byte_endian( (int) mNumberOfBasisPerElement );


            luint tNumberOfAllElementsOnProc = mAllElementsOnProc.size();

            uint tMyRank = par_rank();

            for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
            {
                if ( mAllElementsOnProc( k )->get_owner() == tMyRank )
                {
                    // increment element counter
                    ++tNumberOfElements;
                }
            }

            tFile << "CELLS " << tNumberOfElements << " "
                    << ( mNumberOfBasisPerElement + 1 )*tNumberOfElements  << std::endl;


            // matrix containing node indices
            Mat< luint > tNodes( mNumberOfBasisPerElement, 1 );

            for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
            {
                if ( mAllElementsOnProc( k )->get_owner() == tMyRank )
                {
                    tFile.write( (char*) &tNumberOfNodesPerElement, sizeof(int));

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
            tIChar = swap_byte_endian( (int) 2 );
            for ( luint k = 0; k < tNumberOfElements; ++k)
            {
                tFile.write( (char*) &tIChar, sizeof(int));
            }

            // write node data
            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS BASIS_STATE int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {

                // state flag
                int tState = 0;

                // get state from basis
                if ( mAllBasisOnProc( k )->is_active() )
                {
                    tState = 2;
                }
                else if (mAllBasisOnProc( k )->is_refined() )
                {
                    tState = 1;

                }
                tIChar = swap_byte_endian( tState );


                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            tFile << "SCALARS BASIS_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                tIChar = swap_byte_endian( ( int )  mAllBasisOnProc( k )->get_domain_id() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }


            tFile << "SCALARS BASIS_LEVEL int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                tIChar = swap_byte_endian( ( int )  mAllBasisOnProc( k )->get_level() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }

            tFile << std::endl;

            tFile << "SCALARS BASIS_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                tIChar = swap_byte_endian( ( int )  mAllBasisOnProc( k )->get_owner() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }


            tFile << std::endl;

            /*tFile << "SCALARS BASIS_DOF_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                if ( mAllBasisOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( ( int )   mAllBasisOnProc( k )->get_domain_index() );
                }
                else
                {
                    tIChar = swap_byte_endian( ( int ) -1 );
                }
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl; */

            tFile << "SCALARS BASIS_MEMORY_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                tIChar = swap_byte_endian( ( int )  mAllBasisOnProc( k )->get_memory_index() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            tFile << "SCALARS BASIS_ACTIVE_INDEX int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {
                if ( mAllBasisOnProc( k )->is_active() )
                {
                    tIChar = swap_byte_endian( ( int )   mAllBasisOnProc( k )->get_active_index() );
                }
                else
                {
                    tIChar = swap_byte_endian( ( int ) -1 );
                }
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            // close the output file
            tFile.close();

            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created VTK debug file.\n               Mesh has %lu basis.\n               Creation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) tNumberOfNodes,
                        ( double ) tElapsedTime / 1000 );
            }

        }

    } /* namespace hmr */
} /* namespace moris */
