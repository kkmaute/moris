/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Geometric_Multigrid.cpp
 *
 */

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_DLA_Geometric_Multigrid.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_HMR_Database.hpp"

namespace moris
{
namespace dla
{
    Geometric_Multigrid::Geometric_Multigrid( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ),
                                                                                      mMesh( mSolverInterface->get_mesh_pointer_for_multigrid() )
    {
        // Get the maximal mesh level
        moris::uint tNumBsplineMeshes = mMesh->get_num_interpolations();

        moris::sint tMaxMeshLevel = -1;

        for ( moris::uint Ia = 0; Ia < tNumBsplineMeshes; Ia++ )
        {
            moris::sint tMaxMeshLevelForMeshIndex = mMesh->get_max_level( Ia );

            tMaxMeshLevel = std::max( tMaxMeshLevel, tMaxMeshLevelForMeshIndex );
        }

        // Get the number of dofs per level which equal the current multigrid level or are coarser.
        moris::Matrix< DDUMat > tRemainingOldDofsOnLevel = mSolverInterface->get_number_remaining_dofs();

        // Get maps from MSI
        mListAdofExtIndMap          = mSolverInterface->get_lists_of_ext_index_multigrid();
        mListAdofTypeTimeIdentifier = mSolverInterface->get_lists_of_multigrid_identifiers() ;
        mMultigridMap               = mSolverInterface->get_multigrid_map();

        // Build matrix vector factory to build prolongation operators
        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );

        // Set size of List containing prolongation operators
        mProlongationList.resize( mListAdofExtIndMap.size() - 1 );

        // Loop over all coarse levels.
        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            // Create prolongation matrix
            mProlongationList( Ik - 1 ) = tMatFactory.create_matrix( mListAdofExtIndMap( Ik-1 ).numel(), mListAdofExtIndMap( Ik ).numel() );
        }

        // Loop over all coarse levels
        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            // Loop over coarse dofs
            for ( moris::uint Ii = 0; Ii < mListAdofExtIndMap( Ik ).numel(); Ii++ )
            {
                // Get this dof type
                moris::sint tDofType = mSolverInterface->get_type_time_identifier_to_type_map()( mListAdofTypeTimeIdentifier( Ik )( Ii ) );

                // Get order of this dof typegit a
                moris::sint tMeshIndex = mSolverInterface->get_adof_index_for_type( tDofType );

                // Get external dof id and identifier of this dof
                moris::uint tExtDofInd     = mListAdofExtIndMap( Ik )( Ii, 0 );
                moris::uint tDofIdentifier = mListAdofTypeTimeIdentifier( Ik )( Ii, 0 );

                // Ask mesh for the level of this dof index
                moris::uint tDofLevel = mMesh->get_basis_level( tMeshIndex, tExtDofInd );

                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
                if( ( tDofLevel <= tMaxMeshLevel - Ik ) && ( Ii < tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {
                   // These dofs stay unchanged. Thus assemble a zero to the diagonal
                   moris::Matrix< DDRMat > tIdentityMat( 1, 1, 1.0 );
                   moris::Matrix< DDSMat > tRowMat( 1, 1, Ii );

                   // Get row index of this dof on the fine mesh
                   moris::sint tColLevelPos = mMultigridMap( Ik-1 )( tDofIdentifier )( tExtDofInd, 0 );
                   moris::Matrix< DDSMat > tColMat( 1, 1, tColLevelPos );

                   mProlongationList( Ik-1 )->insert_values( tRowMat, tColMat, tIdentityMat );
                }
                // If coarse dof on this level is interpolated through fine dofs on thislevel + 1
                else if ( ( tDofLevel == tMaxMeshLevel - Ik ) && ( Ii >= tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {
                    // Create vector for row index
                    moris::Matrix< DDSMat > tRowMat( 1, 1, Ii );

                    // Get vector with external fine indices
                    moris::Matrix< DDSMat > tIndices = mMesh->get_fine_basis_inds_of_basis( tMeshIndex, tExtDofInd );

                    // Initialize vector with col indices
                    moris::Matrix< DDSMat > tColMat( tIndices.numel(), 1, -1 );

                    // Map col external fine indices to col index
                    for ( moris::uint Ia = 0; Ia < tIndices.numel(); Ia++ )
                    {
                        tColMat( Ia, 0 ) = mMultigridMap(Ik-1)(0)( tIndices( Ia, 0 ), 0 );           //FIXME add identifier
                    }

                    // Get weights
                    moris::Matrix< DDRMat > tWeights = mMesh->get_fine_basis_weights_of_basis( tMeshIndex, tExtDofInd  );

                    // Fill weights in operator
                    mProlongationList( Ik-1 )->insert_values( tRowMat, tColMat, tWeights );
                }
                else
                {
                    MORIS_ERROR(false, "Geometric_Multigrid::Geometric_Multigrid: Problem with Geometric multigrid. Dof either on a level which is too fine, or coarse but not refined  ");
                }
            }
            mProlongationList( Ik - 1 )->matrix_global_assembly();
//            mProlongationList( Ik - 1 )->print();
        }
    }

    Geometric_Multigrid::~Geometric_Multigrid()
    {
        for ( moris::uint Ik = 0; Ik < mProlongationList.size(); Ik++ )
        {
            delete( mProlongationList( Ik ) );
        }
    }
}
}

