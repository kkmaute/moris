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

#include "cl_MTK_Mesh.hpp"
#include "cl_HMR_Database.hpp"

namespace moris
{
namespace dla
{
    Geometric_Multigrid::Geometric_Multigrid( Solver_Interface * aSolverInterface ) : mSolverInterface( aSolverInterface ),
                                                                                      mMesh( mSolverInterface->get_mesh_pointer_for_multigrid() )
    {
        mListAdofExtIndMap = mSolverInterface->get_lists_of_ext_index_multigrid();

        // Gets the maximal mesh level
        //moris::uint tMaxMeshLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )->get_max_level();


        for ( moris::uint Ik = 1; Ik < mListAdofExtIndMap.size(); Ik++ )
        {
            for ( moris::uint Ii = 1; Ii < mListAdofExtIndMap( Ik ).length(); Ii++ )
            {
//                moris::uint tExtDofInd = mListAdofExtIndMap( Ik )( Ii, 0 );
//
//                // Ask mesh for the level of this mesh index
//                moris::uint tDofLevel = mMesh->get_HMR_database()->get_bspline_mesh_by_index( gAdofOrderHack )
//                                                                 ->get_basis_by_index( mListAdofExtIndMap( Ik )( Ii, 0 ) )
//                                                                 ->get_level();
//
//                // If Index is inside of the set of dofs on this multigrid level, than add it to list.
//                if( tDofLevel < tMaxMeshLevel - Ik )
//                {
//
//                }
//                else
//                {
//
//                }
//
//                //mSolverInterface->read_multigrid_maps()
            }
        }
    }

}
}
