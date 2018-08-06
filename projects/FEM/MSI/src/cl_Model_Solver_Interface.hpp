/*
 * cl_Equation_Manager.hpp
 *
 *  Created on: Jul 22, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_MANAGER_HPP_
#define SRC_FEM_CL_EQUATION_MANAGER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mat.hpp"

#include "cl_Dof_Manager.hpp"

namespace moris
{

    namespace MSI
    {
    class MSI_Solver_Interface;
    class Model_Solver_Interface
    {
    private:
        moris::uint mNumEquationObjects;
        moris::Cell< Equation_Object* > mEquationObjectList;
        Dof_Manager mDofMgn;

        moris::Cell< MSI_Solver_Interface* > mAAA;

    public:
        Model_Solver_Interface( const moris::uint aNumEquationObj,
                                moris::Cell < Equation_Object* > & aListEqnObj ) : mNumEquationObjects( aNumEquationObj ),
                                                                                   mEquationObjectList( aListEqnObj ),
                                                                                   mDofMgn( aNumEquationObj, aListEqnObj )
        {
            //Dof_Manager tDofMgn ( aNumEquationObj, aListEqnObj );
        };

        ~Model_Solver_Interface()
        {};

        moris::uint get_num_eqn_objs()
        {
            return mNumEquationObjects;
        };

        void get_equation_obj_jacobian( const moris::uint               & aEqnObjInd,
                                              moris::Mat< moris::real > & aEqnObjMatrix)
        {
            mEquationObjectList( aEqnObjInd )->get_egn_obj_jacobian( aEqnObjMatrix );
        };

        void get_equation_obj_residual ( const moris::uint               & aEqnObjInd,
                                               moris::Mat< moris::real > & aEqnObjRHS)
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_residual( aEqnObjRHS );
        };

        void get_equation_obj_dof_ids( const moris::uint       & aEqnObjInd,
                                             moris::Mat< int > & aElementTopology )
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_dof_ids( aElementTopology );
        };


        void solve_system();

        void solve_system( moris::Mat< moris::real > & aSolution );

    };
    }
}



#endif /* SRC_FEM_CL_EQUATION_MANAGER_HPP_ */
