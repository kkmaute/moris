/*
 * cl_DLA_Geometric_Multigrid.cpp
 *
 *  Created on: Nov 18, 2018
 *      Author: schmidt */

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_print.hpp"

#include "cl_DLA_Geometric_Multigrid.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Enums.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"

//#include "MSI_Adof_Order_Hack.hpp"

//using namespace MSI;

namespace moris
{
namespace dla
{
    Geometric_Multigrid::Geometric_Multigrid( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ),
                                                                                      mMesh( mSolverInterface->get_mesh_pointer_for_multigrid() )
    {
        PetscInitializeNoArguments();
        //FIXME Inser AdofOrderHack
        moris::uint tAdofOrderHack = 1;

        moris::Matrix< DDUMat > tRemainingOldDofsOnLevel = mSolverInterface->get_number_remaining_dofs();

        mListAdofExtIndMap = mSolverInterface->get_lists_of_ext_index_multigrid();
        mMultigridMap = mSolverInterface->get_multigrid_map();


        Matrix_Vector_Factory tMatFactory( MapType::Petsc );

        mProlongationList.resize( mListAdofExtIndMap.size() - 1 );

        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            // Build matrix
            mProlongationList( Ik - 1 ) = tMatFactory.create_matrix( mListAdofExtIndMap( Ik-1 ).length(), mListAdofExtIndMap( Ik ).length() );
//            mProlongationList( Ik - 1 )->matrix_global_asembly();
//            mProlongationList( Ik - 1 )->print_matrix_to_screen();
        }

        // Gets the maximal mesh level
        moris::uint tMaxMeshLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )->get_max_level();

        std::cout<<tMaxMeshLevel<<" Max level "<<std::endl;
        std::cout<<mListAdofExtIndMap.size()<<" List Size "<<std::endl;

        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            std::cout<<mListAdofExtIndMap( Ik ).length()<<" List Length "<<std::endl;

            for ( moris::uint Ii = 0; Ii < mListAdofExtIndMap( Ik ).length(); Ii++ )
            {
                moris::uint tExtDofInd = mListAdofExtIndMap( Ik )( Ii, 0 );


                // Ask mesh for the level of this mesh index
                moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( tAdofOrderHack )
                                                                 ->get_basis_by_index( tExtDofInd )
                                                                 ->get_level();

                std::cout<<tExtDofInd<<" Ext Dof Ind "<<tDofLevel<<" "<<Ii<<std::endl;


                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
                if( ( tDofLevel <= tMaxMeshLevel - Ik ) && ( Ii < tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {

                }
                else if ( ( tDofLevel == tMaxMeshLevel - Ik ) && ( Ii >= tRemainingOldDofsOnLevel( Ik-1, 0 ) ) )
                {
                    moris::Matrix< DDRMat > tWeights = mMesh->get_HMR_database()
                                                            ->get_bspline_mesh_by_index( tAdofOrderHack )
                                                            ->get_children_weights_for_parent( tExtDofInd );
                    print(tWeights, " weights");

                    moris::Matrix< DDSMat > tIndices = mMesh->get_HMR_database()
                                                            ->get_bspline_mesh_by_index( tAdofOrderHack )
                                                            ->get_children_ind_for_basis( tExtDofInd );
                    print(tIndices, " index");
                }
                else
                {
                	MORIS_ERROR(false, "problem in ");
                }

//                //mSolverInterface->read_multigrid_maps()
            }
        }
    }

}
}
