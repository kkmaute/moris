/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Multigrid.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

#include "fn_print.hpp"

namespace moris
{
    namespace MSI
    {
    Multigrid::Multigrid( moris::MSI::Model_Solver_Interface * aModelSolverInterface,
                          moris::mtk::Mesh                   * aMesh ) : mMesh( aMesh )
    {
        mMultigridLevels = 3;                                                                   //FIXME Input

        moris::Cell < Adof * > tOwnedAdofList = aModelSolverInterface->get_dof_manager()->get_owned_adofs();

        moris::uint tNumOwnedAdofs = tOwnedAdofList.size();

        std::cout<< tNumOwnedAdofs << " tNumOwnedAdofs"<<std::endl;

        mListAdofExtIndMap.resize( mMultigridLevels + 1 );

        mListAdofExtIndMap( 0 ).set_size( tNumOwnedAdofs, 1 );

        moris::uint Ik = 0;
        for ( Adof* tAdof : tOwnedAdofList )
        {
            mListAdofExtIndMap( 0 )( Ik, 0) = tAdof->get_adof_external_ind();

            mListAdofTypeTimeIdentifier( 0 )( Ik++, 0) = tAdof->get_adof_type_time_identifier();
        }

        MORIS_ASSERT( mListAdofTypeTimeIdentifier( 0 ).min() != -1, "moris::MSI::Multigrid: Type/time identifier not specified");
    }

    void Multigrid::create_multigrid_level_dof_ordering()
    {
        // Gets the maximal mesh level
        moris::uint tMaxMeshLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )->get_max_level();

        // Loop over all multigrid levels
        for ( moris::sint Ik = 0; Ik < mMultigridLevels; Ik++ )
        {
            moris::uint tNumDofsOnLevel = mListAdofExtIndMap( Ik ).length();
            moris::uint tCounter = 0;
            moris::uint tCounterTooFine = 0;

            mListAdofExtIndMap( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );
            mListAdofTypeTimeIdentifier( Ik + 1 ).set_size( tNumDofsOnLevel, 1 );

            Matrix< DDSMat > tEntryOfTooFineDofs( tNumDofsOnLevel, 1, -1);

            // Loop over all dofs on this level
            for ( moris::uint Ii = 0; Ii < tNumDofsOnLevel; Ii++ )
            {
                // Ask mesh for the level of this mesh index
                moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )
                                                                 ->get_basis_by_memory_index( mListAdofExtIndMap( Ik )( Ii, 0 ) )
                                                                 ->get_level();

                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
                if( tDofLevel <= tMaxMeshLevel - Ik )
                {
                    mListAdofExtIndMap( Ik + 1 )( tCounter ) = mListAdofExtIndMap( Ik )( Ii, 0 );

                    mListAdofTypeTimeIdentifier( Ik + 1 )( tCounter++ ) = mListAdofTypeTimeIdentifier( Ik )( Ii, 0 );
                }
                else
                {
                    tEntryOfTooFineDofs( tCounterTooFine++, 0 ) = Ii;
                }
            }

//            // Is this nessecary
//            //tEntryOfTooFineDofs.resize( tCounterTooFine, 1 );
//
//
//            // Loop over all refined dofs on this level
//            for ( moris::uint Ii = 0; Ii < tCounterTooFine; Ii++ )
//            {
//                // Ask for refined dofs on this level
//            }
//
//            mListAdofExtIndMap( Ik + 1 ).resize( tCounter, 1 );
//            mListAdofTypeTimeIdentifier( Ik + 1 ).resize( tCounter, 1 );
        }
    }
}
}
