/*
 * cl_HMR_Refinement_Manager.hpp
 *
 *  Created on: May 30, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_
#define SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_

#include "typedefs.hpp" //COR/src
#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Knot.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src


namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------
    class Refinement_Manager
    {

        const Parameters             *   mParameters;
        uint                           mMaxVolumeLevel;
        Background_Mesh_Base       *   mBackgroundMesh;

        //Lagrange_Mesh_Base         *   mLagrangeMesh;
//-------------------------------------------------------------------------------
    public:
//-------------------------------------------------------------------------------

        /**
         * \brief object that manages refinement
         *
         *  @param[in]      aParameters        user settings
         *  @param[in]      aBackgroudMesh   background mesh to operate on
         *  @param[in]      aLagrangeMesh    Lagrange mesh where Data are
         */
        Refinement_Manager (
                const Parameters         * aParameters,
                Background_Mesh_Base   * aBackgroudMesh ) :
                    mParameters( aParameters ),
                    mMaxVolumeLevel( aParameters->get_max_volume_level() ),
                    mBackgroundMesh( aBackgroudMesh )
        {

        }

//-------------------------------------------------------------------------------

        /* destructor for refinement manager */
        ~Refinement_Manager(){}

//-------------------------------------------------------------------------------
        void
        refine_against_function( real (*aFunction)( const real & aX, const real & aY, const real & aZ ) )
        {
            // start timer
            tic tTimer;

            // ask background mesh for number of elements
            luint tNumberOfElements
                = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // get number of dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // get number of corner nodes
            uint tNumberOfNodesPerElement = std::pow( 2, tNumberOfDimensions );

            // reset element counter
            luint tElementCounter = 0;

            Mat< real > tNodes( tNumberOfDimensions, tNumberOfNodesPerElement );

            // loop over all active elements owned by proc
            for( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to Lagrange element
                Background_Element_Base* tElement
                    = mBackgroundMesh->get_element( e );

                // get nodes of element
                mBackgroundMesh->calc_corner_nodes_of_element( tElement, tNodes );

                // a node
                Mat< real > tNode( tNumberOfDimensions, 1 );

                // copy first node
                tNode.cols( 0, 0 ) = tNodes.cols( 0, 0 );

                // calculate value
                real tValue = aFunction(
                                                tNodes( 0, 0 ),
                                                tNodes( 1, 0 ),
                                                tNodes( 2, 0 ) );

                // flag checking if value of first node is < 0
                bool tSign = tValue < 0;

                // loop over all remaining elements
                for ( uint k = 1; k<tNumberOfNodesPerElement; ++k )
                {
                    // calculate value
                    tValue = aFunction(
                            tNodes( 0, k ),
                            tNodes( 1, k ),
                            tNodes( 2, k ) );

                    bool tThisSign  = tValue < 0;

                    // check if node is intersected
                    if ( tThisSign != tSign )
                    {
                        // flag this element
                        tElement->put_on_queue();

                        // increment element counter
                        ++tElementCounter;

                        // exit this loop
                        break;
                    }
                }
            }

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                if ( tElementCounter == 1 )
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( double ) tElapsedTime / 1000 );
                }
                else
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( long unsigned int ) tElementCounter,
                            ( double ) tElapsedTime / 1000 );
                }
            }

            // refine mesh
            this->refine_mesh();
        }

//-------------------------------------------------------------------------------

        void
        refine_against_function( real (*aFunction)( const real & aX, const real & aY ) )
        {
            // start timer
            tic tTimer;

            // ask background mesh for number of elements
            luint tNumberOfElements
                = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // get number of dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // get number of corner nodes
            uint tNumberOfNodesPerElement = std::pow( 2, tNumberOfDimensions );

            // reset element counter
            luint tElementCounter = 0;

            Mat< real > tNodes( tNumberOfDimensions, tNumberOfNodesPerElement );

            // loop over all active elements owned by proc
            for( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to Lagrange element
                Background_Element_Base* tElement
                = mBackgroundMesh->get_element( e );

                // get nodes of element
                mBackgroundMesh->calc_corner_nodes_of_element( tElement, tNodes );

                // a node
                Mat< real > tNode( tNumberOfDimensions, 1 );

                // copy first node
                tNode.cols( 0, 0 ) = tNodes.cols( 0, 0 );

                // calculate value
                real tValue = aFunction(
                        tNodes( 0, 0 ),
                        tNodes( 1, 0 ) );

                // flag checking if value of first node is < 0
                bool tSign = tValue < 0;

                // loop over all remaining elements
                for ( uint k = 1; k<tNumberOfNodesPerElement; ++k )
                {
                    // calculate value
                    tValue = aFunction(
                            tNodes( 0, k ),
                            tNodes( 1, k ) );

                    bool tThisSign  = tValue < 0;

                    // check if node is intersected
                    if ( tThisSign != tSign )
                    {
                        // flag this element
                        tElement->put_on_queue();

                        // increment element counter
                        ++tElementCounter;

                        // exit this loop
                        break;
                    }
                }
            }

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                if ( tElementCounter == 1 )
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( double ) tElapsedTime / 1000 );
                }
                else
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( long unsigned int ) tElementCounter,
                            ( double ) tElapsedTime / 1000 );
                }
            }

            // refine mesh
            this->refine_mesh();
        }


//-------------------------------------------------------------------------------

        void
        refine_against_functions(
                real (*aVolumeFunction) ( const real & aX, const real & aY, const real & aZ ),
                real (*aSurfaceFunction)( const real & aX, const real & aY, const real & aZ  )
                  )
        {
            // start timer
            tic tTimer;

            // ask background mesh for number of elements
            luint tNumberOfElements
            = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // get number of dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // get number of corner nodes
            uint tNumberOfNodesPerElement = std::pow( 2, tNumberOfDimensions );

            // matrix containing node coordinates
            Mat< real > tNodes( tNumberOfDimensions, tNumberOfNodesPerElement );

            // ask settings for maximum levels
            auto tMaxVolumeLevel = mParameters->get_max_volume_level();

            //auto tMaxSurfaceLevel = mParameters->get_max_surface_level();

            // reset element counter
            luint tElementCounter = 0;

            // loop over all active elements owned by proc
            for( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to Lagrange element
                Background_Element_Base* tElement
                    = mBackgroundMesh->get_element( e );

                // get nodes of element
                mBackgroundMesh->calc_corner_nodes_of_element( tElement, tNodes );

                // get level of element
                uint tLevel = tElement->get_level();

                Mat< real > tValues( 8, 1 );

                if ( tLevel < tMaxVolumeLevel )
                    for( uint k=0; k<8; ++k )
                    {
                        tValues( k ) = aVolumeFunction(
                                tNodes( 0, k ),
                                tNodes( 1, k ),
                                tNodes( 2, k ) );

                    }
                else
                {
                    for( uint k=0; k<8; ++k )
                    {
                        tValues( k ) = aSurfaceFunction(
                                tNodes( 0, k ),
                                tNodes( 1, k ),
                                tNodes( 2, k ) );

                    }
                }

                real tMinVal = tValues.min();
                real tMaxVal = tValues.max();
                // test if element is intersected

                // test if this element is inside
                if (  tMinVal <= 0 && tMaxVal >= 0  && tLevel < tMaxVolumeLevel )
                {
                    // flag this element
                    tElement->put_on_queue();

                    // increment element counter
                    ++tElementCounter;
                }
                else if ( tMinVal <= 0 && tLevel >= tMaxVolumeLevel )
                    {
                        // flag this element
                        tElement->put_on_queue();

                        // increment element counter
                        ++tElementCounter;
                    }
                //}
            }

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                if ( tElementCounter == 1 )
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( double ) tElapsedTime / 1000 );
                }
                else
                {
                    // print output
                    std::fprintf( stdout,"%s Checked refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( long unsigned int ) tElementCounter,
                            ( double ) tElapsedTime / 1000 );
                }
            }

            // refine mesh
            this->refine_mesh();
        }

//-------------------------------------------------------------------------------

        /**
         * a demo function
         */
        void
        refine_against_knot( const Knot & aKnot )
        {
            // start timer
            tic tTimer;

            // ask background mesh for number of elements
            luint tNumberOfElements
            = mBackgroundMesh->get_number_of_active_elements_on_proc();

            // get number of dimensions
            uint tNumberOfDimensions = mParameters->get_number_of_dimensions();

            // get number of corner nodes
            uint tNumberOfNodesPerElement = std::pow( 2, tNumberOfDimensions );

            // reset element counter
            luint tElementCounter = 0;

            Mat< real > tNodes( tNumberOfDimensions, tNumberOfNodesPerElement );

            // ask settings for maximum levels
            auto tMaxVolumeLevel = mParameters->get_max_volume_level();

            auto tMaxSurfaceLevel = mParameters->get_max_surface_level();

            // loop over all active elements owned by proc
            for( luint e=0; e<tNumberOfElements; ++e )
            {

                // get pointer to Lagrange element
                Background_Element_Base* tElement
                    = mBackgroundMesh->get_element( e );

                // get nodes of element
                mBackgroundMesh->calc_corner_nodes_of_element( tElement, tNodes );

                // get level of element
                uint tLevel = tElement->get_level();

                Mat< real > tValues( 8, 1 );

                if ( tLevel < tMaxVolumeLevel )
                    for( uint k=0; k<8; ++k )
                    {
                        tValues( k ) =
                                aKnot.test_for_intersection(
                                tNodes( 0, k ),
                                tNodes( 1, k ),
                                tNodes( 2, k ) );

                    }
                else
                {
                    for( uint k=0; k<8; ++k )
                    {
                        tValues( k ) = aKnot.calculate_level_set(
                                tNodes( 0, k ),
                                tNodes( 1, k ),
                                tNodes( 2, k ) );

                    }
                }
                real tMinVal = tValues.min();
                real tMaxVal = tValues.max();
                // test if this element is inside
                if ( tMaxVal <= 0  && tLevel < tMaxVolumeLevel )
                {
                    // flag this element
                    tElement->put_on_queue();

                    // increment element counter
                    ++tElementCounter;
                }

                // only do this if element is not flagged already
                if ( ! tElement->is_queued() )
                {
                    // test if this element is intersected
                    if ( tMinVal <= 0 && tMaxVal >= 0  && tLevel < tMaxSurfaceLevel )
                    {
                        // flag this element
                        tElement->put_on_queue();

                        // increment element counter
                        ++tElementCounter;
                    }
                }
            }

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                if ( tElementCounter == 1 )
                {
                    // print output
                    std::fprintf( stdout,"%s Checked KNOT refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( double ) tElapsedTime / 1000 );
                }
                else
                {
                    // print output
                    std::fprintf( stdout,"%s Checked KNOT refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                            proc_string().c_str(),
                            ( long unsigned int ) tNumberOfElements,
                            ( long unsigned int ) tElementCounter,
                            ( double ) tElapsedTime / 1000 );
                }
            }

            // refine mesh
            this->refine_mesh();
        }

//-------------------------------------------------------------------------------
    private:
//-------------------------------------------------------------------------------
        /**
         * updates background and Lagrange mesh
         */
        void
        refine_mesh()
        {
            mBackgroundMesh->perform_refinement();
        }

        }; /* Refinement_Manager */
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_REFINEMENT_MANAGER_HPP_ */
