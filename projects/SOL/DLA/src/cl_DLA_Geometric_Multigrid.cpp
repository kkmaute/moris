/*
 * cl_DLA_Geometric_Multigrid.cpp
 *
 *  Created on: Dez 12, 2018
 *      Author: schmidt */

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector.hpp"

#include "cl_DLA_Geometric_Multigrid.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Enums.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_HMR_Database.hpp"

namespace moris
{
namespace dla
{
    Geometric_Multigrid::Geometric_Multigrid( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ),
                                                                                      mMesh( mSolverInterface->get_mesh_pointer_for_multigrid() )
    {
        //FIXME Inser AdofOrderHack
        moris::uint tAdofOrderHack = 0;

        // Get the number of dofs per level which equal the current multigrid level or are coarser.
        moris::Matrix< DDUMat > tRemainingOldDofsOnLevel = mSolverInterface->get_number_remaining_dofs();

        // Get maps from MSI
        mListAdofExtIndMap          = mSolverInterface->get_lists_of_ext_index_multigrid();
        mListAdofTypeTimeIdentifier = mSolverInterface->get_lists_of_multigrid_identifiers() ;
        mMultigridMap               = mSolverInterface->get_multigrid_map();

        // Build matrix vector factory to build prolongation operators
        Matrix_Vector_Factory tMatFactory( MapType::Petsc );

        // Set size of List containing prolongation operators
        mProlongationList.resize( mListAdofExtIndMap.size() - 1 );

        // Loop over all coarse levels.
        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            // Create prolongation matrix
            mProlongationList( Ik - 1 ) = tMatFactory.create_matrix( mListAdofExtIndMap( Ik-1 ).length(), mListAdofExtIndMap( Ik ).length() );
        }

        // Get the maximal mesh level
        moris::uint tMaxMeshLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )->get_max_level();

        // Loop over all coarse levels
        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            // Loop over coarse dofs
            for ( moris::uint Ii = 0; Ii < mListAdofExtIndMap( Ik ).length(); Ii++ )
            {
                // Get external dof id and identifier of this dof
                moris::uint tExtDofInd     = mListAdofExtIndMap( Ik )( Ii, 0 );
                moris::uint tDofIdentifier = mListAdofTypeTimeIdentifier( Ik )( Ii, 0 );

                // Ask mesh for the level of this dof index
                moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )
                                                                 ->get_basis_by_index( tExtDofInd )
                                                                 ->get_level();

                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
                if( ( tDofLevel <= tMaxMeshLevel - Ik ) && ( Ii < tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {
                   // These dofs stay unchanged. Thus assemble a zero to the diagonal
                   moris::Matrix< DDRMat > tIdentityMat( 1, 1, 1.0 );
                   moris::Matrix< DDSMat > tRowMat( 1, 1, Ii );

                   // Get row index of this dof on the fine mesh
                   moris::sint tColLevelPos = mMultigridMap( Ik-1 )( tDofIdentifier )( tExtDofInd, 0 );
                   moris::Matrix< DDSMat > tColMat( 1, 1, tColLevelPos );

                   mProlongationList( Ik-1 )->fill_matrix_row( tIdentityMat, tRowMat, tColMat );
                }
                // If coarse dof on this level is interpolated through fine dofs on thislevel + 1
                else if ( ( tDofLevel == tMaxMeshLevel - Ik ) && ( Ii >= tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {
                    // Create vector for row index
                    moris::Matrix< DDSMat > tRowMat( 1, 1, Ii );

                    // Get vector with external fine indices
                    moris::Matrix< DDSMat > tIndices = mMesh->get_HMR_database()
                                                            ->get_bspline_mesh_by_index( tAdofOrderHack )
                                                            ->get_children_ind_for_basis( tExtDofInd );

                    // Initialize vector with col indices
                    moris::Matrix< DDSMat > tColMat( tIndices.length(), 1, -1 );

                    // Map col external fine indices to col index
                    for ( moris::uint Ia = 0; Ia < tIndices.length(); Ia++ )
                    {
                        tColMat( Ia, 0 ) = mMultigridMap(Ik-1)(0)( tIndices( Ia, 0 ), 0 );
                    }

                    // Get weights
                    moris::Matrix< DDRMat > tWeights = mMesh->get_HMR_database()
                                                            ->get_bspline_mesh_by_index( tAdofOrderHack )
                                                            ->get_children_weights_for_parent( tExtDofInd );

                    // Fill weights in operator
                    mProlongationList( Ik-1 )->fill_matrix_row( tWeights, tRowMat, tColMat );
                }
                else
                {
                    MORIS_ERROR(false, "Geometric_Multigrid::Geometric_Multigrid: Problem with Geometric multigrid. Dof either on a level which is too fine, or coarse but not refined  ");
                }
            }
            mProlongationList( Ik - 1 )->matrix_global_assembly();
            //mProlongationList( Ik - 1 )->print();
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
