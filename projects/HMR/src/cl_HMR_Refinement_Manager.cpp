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
            const Mat< real > & aField,
            const real          aLowerBound,
            const real          aUpperBound )
    {
        // start timer
        tic tTimer;

        MORIS_ERROR( aLowerBound <= aUpperBound,
                     "aLowerBound must be less or equal to aUpperBound " );

        // make sure that field has correct length
        MORIS_ERROR( aField.length() == mMesh->get_number_of_nodes_on_proc(),
                "Field length does not match number of active nodes on proc" );

        // activate pattern on background mesh
        mMesh->select_activation_pattern();

        // get number of elements
        luint tNumberOfElements = mMesh->get_number_of_elements();

        // number of nodes per element
        uint tNumberOfNodes = mMesh->get_number_of_basis_per_element();

        // matrix with nodal values
        Mat< real > tField( tNumberOfNodes, 1 );

        // get max level
        uint tMaxLevel = std::max( mMaxSurfaceLevel, mMaxVolumeLevel );

        // element counter
        luint tElementCounter = 0;

        // loop over all elements
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
                    tField( k ) = aField( tElement->get_basis( k )->get_index() );
                }

                // flag telling if element has been flagged
                bool tElementIsProcessed = false;

                // test volume criterion
                if (  tElement->get_level() < mMaxVolumeLevel )
                {
                    if ( tField.min() <= aLowerBound )
                    {
                        tElement->get_background_element()->put_on_queue();
                        tElementIsProcessed = true;
                        ++tElementCounter;
                    }
                }

                // test surface criterion
                if (  tElement->get_level() < mMaxSurfaceLevel && ! tElementIsProcessed )
                {
                    if ( tField.min() <= aLowerBound && tField.max() >= aUpperBound )
                    {
                        tElement->get_background_element()->put_on_queue();
                        ++tElementCounter;
                    }
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

    void
       Refinement_Manager::flag_against_elemental_field(
               const Mat< real > & aField,
               const real          aLowerBound)
       {
           // start timer
           tic tTimer;

           // make sure that field has correct length
           MORIS_ERROR( aField.length() == mMesh->get_number_of_elements(),
                   "Field length does not match number of active elements on proc" );

           // activate pattern on background mesh
           mMesh->select_activation_pattern();

           // get number of elements
           luint tNumberOfElements = mMesh->get_number_of_elements();

           // number of nodes per element
           uint tNumberOfNodes = mMesh->get_number_of_basis_per_element();

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
                   if ( aField( e ) <= aLowerBound )
                   {
                       tElement->get_background_element()->put_on_queue();
                       tElementIsProcessed = true;
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
