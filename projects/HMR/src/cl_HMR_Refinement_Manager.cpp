/*
 * cl_HMR_Refinement_Manager.cpp
 *
 *  Created on: May 30, 2018
 *      Author: messe
 */

#include "assert.hpp"
#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_HMR_Refinement_Manager.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------

        Refinement_Manager::Refinement_Manager ( Lagrange_Mesh_Base * aMesh ) :
                mMesh ( aMesh ),
                mMaxSurfaceLevel( aMesh->get_parameters()->get_max_surface_level() ),
                mMaxVolumeLevel( aMesh->get_parameters()->get_max_volume_level() )
        {

        }

//-------------------------------------------------------------------------------

    void
    Refinement_Manager::flag_against_nodal_field(
            const Mat< real > & aNodeValues,
            const real          aLowerBound,
            const real          aUpperBound )
    {
        // start timer
        tic tTimer;

        MORIS_ERROR( aLowerBound <= aUpperBound,
                     "aLowerBound must be less or equal to aUpperBound " );

        // make sure that field has correct length
        MORIS_ERROR( aNodeValues.length() == mMesh->get_number_of_nodes_on_proc(),
                "Field length does not match number of active nodes on proc" );

        // activate pattern on background mesh
        mMesh->select_activation_pattern();

        // get number of elements
        luint tNumberOfElements = mMesh->get_number_of_elements();

        // number of nodes per element
        uint tNumberOfNodes = mMesh->get_number_of_basis_per_element();

        MORIS_ERROR( aNodeValues.length() == mMesh->get_number_of_nodes_on_proc(),
                "Number of nodes does not match" );

        // matrix with nodal values
        Mat< real > tField( tNumberOfNodes, 1 );

        // get max level
        uint tMaxLevel = std::max( mMaxSurfaceLevel, mMaxVolumeLevel );

        // element counter
        luint tElementCounter = 0;

        // loop over all elements for coarsening
        for( luint e=0; e<tNumberOfElements; ++e )
        {
            // get pointer to Lagrange Element
            auto tElement = mMesh->get_element( e );

            // test if element is of interest
            if ( tElement->get_level() < tMaxLevel )
            {
                // copy nodal values into tField
                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    tField( k ) = aNodeValues( tElement->get_basis( k )->get_index() );
                }

                // flag telling if element has been flagged
                bool tElementIsProcessed = false;

                // test volume criterion
                if (  tElement->get_level() < mMaxVolumeLevel )
                {
                    if ( tField.min() <= aLowerBound )
                    {
                        tElement->get_background_element()->put_on_refinement_queue();
                        tElementIsProcessed = true;
                        ++tElementCounter;
                    }
                }

                // test surface criterion
                if (  tElement->get_level() < mMaxSurfaceLevel && ! tElementIsProcessed )
                {
                    if ( tField.min() <= aLowerBound && tField.max() >= aUpperBound )
                    {
                        tElement->get_background_element()->put_on_refinement_queue();
                        ++tElementCounter;
                    }
                }
            }

        }

        // fixme: make this also work for parallel
        /*if( par_size() == 1 )
        {
        // determine number of children
        uint tNumberOfChildrenPerElement = std::pow( 2, mMesh->get_parameters()->get_number_of_dimensions() );

        // initialize parent field
        Mat< real > tParentField( tNumberOfChildrenPerElement*tNumberOfNodes, 1 );

        // loop over all elements for refining
        for( luint e=0; e<tNumberOfElements; ++e )
        {
            // get pointer to Lagrange Element
            auto tElement = mMesh->get_element( e );

            // test if element is of interest
            if ( tElement->get_level() > 0 && ! tElement->get_background_element()->is_queued_for_coarsening() )
            {


                // get parent of element
                auto tParent = tElement->get_background_element()->get_parent();

                // check if any child is marked for refinement
                bool tChildRefine = false;
                for( uint c=0; c<tNumberOfChildrenPerElement; ++c )
                {
                    if ( tParent->get_child( c )->is_queued_for_refinement() )
                    {
                        tChildRefine = true;
                        break;
                    }
                }

                // only do this if no child is flagged for refinement
                if ( ! tChildRefine )
                {
                    // initialize counter
                    uint tCount = 0;

                    for( uint c=0; c<tNumberOfChildrenPerElement; ++c )
                    {
                        // get child
                        auto tChild = mMesh->get_element_by_memory_index( tParent->get_child( c )->get_memory_index() );

                        for( uint k=0; k<tNumberOfNodes; ++k )
                        {
                            tParentField( tCount++ ) = aField( tChild->get_basis( k )->get_index() );
                        }
                    }

                    // check if all sibling are not intersected
                    if ( tParentField.min() > aUpperBound || tParentField.max() < aLowerBound )
                    {
                        // flag all children for coarsening
                        for( uint c=0; c<tNumberOfChildrenPerElement; ++c )
                        {
                            tParent->get_child( c )->put_on_coarsening_queue();
                        }
                    }
                }
            }
        }
        } // end par_size == 1 */

        // print a debug statement if verbosity is set
        if (  mMesh->get_parameters()->is_verbose() )
        {
            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            if ( tElementCounter == 1 )
            {
                // print output
                std::fprintf( stdout,"%s Checked nodal refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) tNumberOfElements,
                        ( double ) tElapsedTime / 1000 );
            }
            else
            {
                // print output
                std::fprintf( stdout,"%s Checked nodal refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) tNumberOfElements,
                        ( long unsigned int ) tElementCounter,
                        ( double ) tElapsedTime / 1000 );
            }
        }
    }

//-------------------------------------------------------------------------------

    void
       Refinement_Manager::flag_against_elemental_field(
               const Mat< real > & aElementValues,
               const real          aLowerBound)
       {
           // start timer
           tic tTimer;

           // make sure that field has correct length
           MORIS_ERROR( aElementValues.length() == mMesh->get_number_of_elements(),
                   "Field length does not match number of active elements on proc" );

           // activate pattern on background mesh
           mMesh->select_activation_pattern();

           // get number of elements
           luint tNumberOfElements = mMesh->get_number_of_elements();

           // element counter
           luint tElementCounter = 0;

           // loop over all elements
           for( luint e=0; e<tNumberOfElements; ++e )
           {
               // get pointer to Lagrange Element
               auto tElement = mMesh->get_element( e );

               // test volume criterion
               if (  tElement->get_level() < mMaxVolumeLevel )
               {
                   if ( aElementValues( e ) <= aLowerBound )
                   {
                       tElement->get_background_element()->put_on_refinement_queue();
                       ++tElementCounter;
                   }
               }
           }

           // print a debug statement if verbosity is set
           if (  mMesh->get_parameters()->is_verbose() )
           {
               // stop timer
               real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

               if ( tElementCounter == 1 )
               {
                   // print output
                   std::fprintf( stdout,"%s Checked nodal refinement criterion for %lu elements.\n               Flagged 1 element.\n               Took %5.3f seconds.\n\n",
                           proc_string().c_str(),
                           ( long unsigned int ) tNumberOfElements,
                           ( double ) tElapsedTime / 1000 );
               }
               else
               {
                   // print output
                   std::fprintf( stdout,"%s Checked nodal refinement criterion for %lu elements.\n               Flagged %lu elements.\n               Took %5.3f seconds.\n\n",
                           proc_string().c_str(),
                           ( long unsigned int ) tNumberOfElements,
                           ( long unsigned int ) tElementCounter,
                           ( double ) tElapsedTime / 1000 );
               }
           }
       }

//-------------------------------------------------------------------------------

    }
} /* namespace moris */
