
#include "cl_NLA_Nonlinear_Solver_Manager.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_NLA_Nonlinear_Database.hpp"

using namespace moris;
using namespace NLA;

void Nonlinear_Database::create_solver_manager_dependencies()
{
    mListSolverManagerDepenencies.resize( mListNonlinerSolverManagers.size() );

    for ( uint Ik = 0 ; Ik < mListNonlinerSolverManagers.size(); Ik++ )
    {
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tNonlinerSolverManagerDofTypeList = mListNonlinerSolverManagers( Ik )->get_dof_type_list();

        mListSolverManagerDepenencies( Ik ).set_size( tNonlinerSolverManagerDofTypeList.size(), 1, -1 );

        for ( uint Ii = 0 ; Ii < tNonlinerSolverManagerDofTypeList.size(); Ii++ )
        {
            for ( uint Ij = 0 ; Ij < mListNonlinerSolverManagers.size(); Ij++ )
            {
                if( Ij != Ik )
                {
                    moris::Cell< enum MSI::Dof_Type > tUnionDofTypeList = mListNonlinerSolverManagers( Ij )->get_dof_type_union();

                    if( tNonlinerSolverManagerDofTypeList( Ii ).data() == tUnionDofTypeList.data() )
                    {
                        MORIS_ERROR( mListSolverManagerDepenencies( Ik )( Ii, 0 ) == -1,
                                "Nonlinear_Database::create_solver_manager_dependencies(): Two solvers are operating on the same dof types" );

                        mListSolverManagerDepenencies( Ik )( Ii, 0 ) = Ij;
                    }
                }
            }
        }
    }
}

//---------------------------------------------------------------------------------------------------------------------

moris::sint Nonlinear_Database::get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
                                                                    const moris::sint aDofTypeListIndex )
{
    MORIS_ERROR( mListSolverManagerDepenencies( aSolverManagerIndex )( aDofTypeListIndex, 0 ) != -1,
            "Nonlinear_Database::get_nonliner_solver_manager_index(): Return nonlinear solver manager -1. Check create_solver_manager_dependencies()" );

    return mListSolverManagerDepenencies( aSolverManagerIndex )( aDofTypeListIndex , 0 );
}

//---------------------------------------------------------------------------------------------------------------------



