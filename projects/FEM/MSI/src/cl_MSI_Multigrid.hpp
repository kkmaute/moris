/*
 * cl_MSI_Multigrid.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_MULTIGRID_HPP_
#define SRC_FEM_CL_MSI_MULTIGRID_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;
    }
    namespace MSI
    {
        class Model_Solver_Interface;
        class Multigrid
        {
        private:
            moris::sint mMultigridLevels = -1;

            mtk::Mesh * mMesh;

            moris::Cell< Matrix< DDUMat > > mListAdofExtIndMap;   // List of fine of coarse external index
            moris::Cell< Matrix< DDSMat > > mListAdofTypeTimeIdentifier;   // List of type time identifiers for carse and fine mesh


        public:
           // Multigrid(){};

            Multigrid( moris::MSI::Model_Solver_Interface * aModelSolverInterface,
                       moris::mtk::Mesh                   * aMesh );
//            {
//                mMultigridLevel = 3;          //FIXME
//
//                moris::Cell < Adof * > tOwnedAdofList= aModelSolverInterface->get_dof_manager()->get_owned_adofs();
//
//                moris::uint tNumOwnedAdofs = tOwnedAdofList.size();
//
//                mListAdofExtIndMap.resize( mMultigridLevels + 1 );
//
//                mListAdofExtIndMap( 0 ).set_size( tNumOwnedAdofs, 1 );
//
//                moris::uint Ik = 0;
//                for ( Adof* tAdof : tOwnedAdofList )
//                {
//                    mListAdofExtIndMap( 0 )( Ik, 0) = tAdof->get_adof_external_ind();
//
//                    mListAdofTypeTimeIdentifier( 0 )( Ik++, 0) = tAdof->get_adof_type_time_identifier();
//                }
//
//                MORIS_ASSERT( mListAdofTypeTimeIdentifier( 0 ).min() != -1, "moris::MSI::Multigrid: Type/time identifier not specified");
//
//            };

            ~Multigrid(){};

            void create_multigrid_level_dof_ordering();
//            {
//                // Gets the maximal mesh level
//                //moris::uint tMaxMeshLevel = ;
//
//                // Loop over all multigrid levels
//                for ( moris::sint Ik = 0; Ik < mMultigridLevels; Ik++ )
//                {
//                    moris::uint tNumDofsOnLevel = mListAdofExtIndMap( Ik ).size();
//                    moris::uint tCounter = 0;
//                    moris::uint tCounterTooFine = 0;
//
//                    mListAdofExtIndMap( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );
//                    mListAdofTypeTimeIdentifier( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );
//
//                    Matrix< DDSMat > tEntryOfTooFineDofs( tNumDofsOnLevel, 1, -1);
//
//                    // Loop over all dofs on this level
//                    for ( moris::uint Ii = 0; Ii < tNumDofsOnLevel; Ii++ )
//                    {
//                        // Ask mesh for the level of this mesh index
//                        //moris::uint tDofLevel =                       ( mListAdofExtIndMap( Ik )( Ii, 0 ) );
//
//                        // If Index is inside of the set of dofs on this multigrid level, than add it to list.
//                        if( tDofLevel <= tMaxMeshLevel - Ik )
//                        {
//                            mListAdofExtIndMap( Ik + 1 )( tCounter ) = mListAdofExtIndMap( Ik )( Ii, 0 );
//
//                            mListAdofTypeTimeIdentifier( Ik + 1 )( tCounter++ ) = mListAdofTypeTimeIdentifier( Ik )( Ii, 0 );
//                        }
//                        else
//                        {
//                            tEntryOfTooFineDofs( tCounterTooFine++, 0 ) = Ii;
//                        }
//                    }
//
//                    // Is this nessecary
//                    //tEntryOfTooFineDofs.resize( tCounterTooFine, 1 );
//
//
//                    // Loop over all refined dofs on this level
//                    for ( moris::uint Ii = 0; Ii < tCounterTooFine; Ii++ )
//                    {
//                        // Ask for refined dofs on this level
//                    }
//
//                    mListAdofExtIndMap( Ik + 1 ).resize( tCounter, 1 );
//                    mListAdofTypeTimeIdentifier( Ik + 1 ).resize( tCounter, 1 );
//                }
//            };

        };
    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_MULTIGRID_HPP_ */
