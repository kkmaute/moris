#include "cl_NLA_Nonlinear_Solver_Manager.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_NLA_Nonlinear_Database.hpp"

#include "cl_Vector.hpp"
#include "cl_Map_Class.hpp"

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

void Nonlinear_Database::create_maps()
{

    Matrix_Vector_Factory    tMatFactory( MapType::Epetra );

    moris::uint tNumberNonlinSolverManager = mListNonlinerSolverManagers.size();

    mListOfFreeMaps.resize( tNumberNonlinSolverManager + 1, nullptr );

    // create map object
    mListOfFreeMaps( tNumberNonlinSolverManager ) = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map() );

    for ( uint Ik = 0 ; Ik < tNumberNonlinSolverManager; Ik++ )
    {
        moris::Cell< enum MSI::Dof_Type > tUnionDofTypeList = mListNonlinerSolverManagers( Ik )->get_dof_type_union();

        mSolverInterface->set_requested_dof_types( tUnionDofTypeList );

        mListOfFreeMaps( Ik ) = tMatFactory.create_map( mSolverInterface->get_my_local_global_map() );
    }
}

//---------------------------------------------------------------------------------------------------------------------

Map_Class * Nonlinear_Database::get_list_of_maps( const moris::sint aSolverManagerIndex )
{
    return mListOfFreeMaps( aSolverManagerIndex );
}

//---------------------------------------------------------------------------------------------------------------------

void Nonlinear_Database::finalize()
{
    this->create_solver_manager_dependencies();

    this->create_maps();

    Matrix_Vector_Factory    tMatFactory( MapType::Epetra );

    mFullVector = tMatFactory.create_vector( mSolverInterface, mListOfFreeMaps( mListNonlinerSolverManagers.size()), VectorType::FREE );
    mFullVector->vec_put_scalar( 0.0 );

    mListNonlinerSolverManagers( 0 )->finalize();

    for ( uint Ik = 0; Ik < mListNonlinerSolverManagers.size(); Ik++ )
    {
        mListNonlinerSolverManagers( Ik )->set_nonlinear_manager( this );
    }
}

