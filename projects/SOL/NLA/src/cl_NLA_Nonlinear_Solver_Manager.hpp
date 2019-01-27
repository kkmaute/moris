/*
 * cl_NLA_Nonlinear_Solver_Manager.hpp
 *
 *  Created on: Okt 6, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_MANAGER_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_MANAGER_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "typedefs.hpp" // CON/src
#include "cl_Cell.hpp"
#include <memory>
#include "cl_Param_List.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

namespace moris
{
namespace NLA
{
    class Nonlinear_Problem;
    class Nonlinear_Solver;
    class Nonlinear_Database;
    class Nonlinear_Solver_Manager
    {
    private:
        moris::Cell< moris::Cell< enum MSI::Dof_Type > >  mStaggeredDofTypeList;

        moris::sint mLevel = 0;

        moris::Cell< std::shared_ptr< Nonlinear_Solver > > mNonLinearSolverList;

        moris::uint mCallCounter = 0;

        bool mSolveMonolithically = true;

        Param_List< boost::variant< bool, sint, real > > mParameterListNonLinearSolver;

        enum NonlinearSolverType mNonLinSolverType = NonlinearSolverType::END_ENUM;

        Nonlinear_Database * mNonlinearManager;

        Nonlinear_Problem * mNonlinerProblem = nullptr;

        Solver_Interface * mSolverInput = nullptr;

        dla::Linear_Solver_Manager * mLinSolManager;



    protected:

    public:
        Nonlinear_Solver_Manager( const enum NonlinearSolverType aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

        ~Nonlinear_Solver_Manager();

        void set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
                                const moris::sint                       aLevel =  0)
        { mStaggeredDofTypeList.push_back( aStaggeredDofTypeList ); };

        void set_solver_level( const moris::sint aLevel )    { mLevel = aLevel;};

        Nonlinear_Database * get_nonlinear_database(  )    { return mNonlinearManager;};

        moris::sint get_solver_level()      {  return mLevel;};

        moris::Cell< moris::Cell< enum MSI::Dof_Type > > & get_dof_type_list()    { return mStaggeredDofTypeList; };

        void set_nonlinear_solver( std::shared_ptr< Nonlinear_Solver > aNonLinSolver );

        void set_nonlinear_solver( const moris::uint                         aListEntry,
                                         std::shared_ptr< Nonlinear_Solver > aLinSolver );

        void set_nonlinear_manager( Nonlinear_Database * aNonlinearManager );

        void solve();

        void solve( Nonlinear_Problem * aNonlinearProblem );

        void finalize()
        {

        };


        void set_nonlinear_solver_manager_parameters();

        boost::variant< bool, sint, real > & set_param( char const* aKey )
        {
            return mParameterListNonLinearSolver( aKey );
        }

        moris::Cell< enum MSI::Dof_Type > get_dof_type_union()
        {
            moris::sint tCounter = 0;

            for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
            {
                tCounter = tCounter + mStaggeredDofTypeList( Ik ).size();
            }
;
            moris::Cell< enum MSI::Dof_Type > tUnionEnumList( tCounter );
            tCounter = 0;

            for ( moris::uint Ik = 0; Ik < mStaggeredDofTypeList.size(); ++Ik )
            {
                for ( moris::uint Ii = 0; Ii < mStaggeredDofTypeList( Ik ).size(); ++Ii )
                {
                    tUnionEnumList( tCounter++ ) = mStaggeredDofTypeList( Ik )( Ii );
                }
            }

            return tUnionEnumList;
        };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_MANAGER_HPP_ */

