/*
 * cl_NLA_Nonlinear_Database.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Manager.hpp"

namespace moris
{
namespace NLA
{
    class Nonlinear_Problem;
    class Nonlinear_Solver;
    class Nonlinear_Database
    {
    private:
        moris::Cell< Nonlinear_Solver_Manager * > mListNonlinerSolverManagers;

        moris::Cell< moris::Matrix< DDSMat > > mListSolverManagerDepenencies;

        Solver_Interface * mSolverInput;

        void create_solver_manager_dependencies();

    protected:

    public:
        Nonlinear_Database( Solver_Interface * aSolverInput ) : mSolverInput( aSolverInput )
        {};

        ~Nonlinear_Database()
        {};

        void set_nonliner_solver_managers( Nonlinear_Solver_Manager * aNonlinerSolverManager )
        {
            mListNonlinerSolverManagers.push_back( aNonlinerSolverManager );

            mListNonlinerSolverManagers( mListNonlinerSolverManagers.size() -1 )->set_sonlinear_solver_manager_index( mListNonlinerSolverManagers.size() -1 );
        };

        void solve()
        {
            this->finalize();

            mListNonlinerSolverManagers( 0 )->solve( );
        };

        void finalize()
        {
            this->create_solver_manager_dependencies();

            mListNonlinerSolverManagers( 0 )->finalize();

            for ( uint Ik = 0; Ik < mListNonlinerSolverManagers.size(); Ik++ )
            {
                mListNonlinerSolverManagers( Ik )->set_nonlinear_manager( this );
            }
        };

        Solver_Interface * get_solver_interface()
        {
            return mSolverInput;
        };

        moris::Cell< Nonlinear_Solver_Manager * > & get_nonliner_solver_manager_list()
        {
            return mListNonlinerSolverManagers;
        };

        moris::sint get_nonlinear_solver_manager_index( const moris::sint aSolverManagerIndex,
                                                        const moris::sint aDofTypeListIndex );

    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_DATABASE_HPP_ */

