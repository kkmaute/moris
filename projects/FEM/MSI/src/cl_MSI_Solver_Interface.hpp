/*
 * cl_MSI_Solver_Interface.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_

#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_Solver_Interface.hpp"

namespace moris
{
class Dist_Vector;
    namespace MSI
    {
    class MSI_Solver_Interface : public moris::Solver_Interface
    {
    private:
        moris::MSI::Model_Solver_Interface*   mMSI;
        moris::MSI::Dof_Manager*   mDofMgn;

        Dist_Vector * mSolutionVector;

    public:
        MSI_Solver_Interface( )
        {};

        MSI_Solver_Interface( moris::MSI::Model_Solver_Interface * aMSI,
                              moris::MSI::Dof_Manager            * aDofMgn ) : mMSI( aMSI ),
                                                                               mDofMgn( aDofMgn )
        {};

        ~MSI_Solver_Interface()
        {};

        void set_solution_vector( Dist_Vector * aSolutionVector )
        {
            mSolutionVector = aSolutionVector;
        }

        // ----------------------------------------------------------------------------------------------
         // local dimension of the problem
         moris::uint get_num_my_dofs()
         {
             moris::uint tNumDofs = mDofMgn->get_num_adofs();
             return tNumDofs;
         };

         // ----------------------------------------------------------------------------------------------
         // local-to-global map
         Matrix< DDSMat > get_my_local_global_map()
         {
             Matrix< DDSMat > tLocalAdofIds = mDofMgn->get_local_adof_ids();
             return tLocalAdofIds;
         };

         Matrix< DDSMat > get_my_local_global_overlapping_map( )
         {
             Matrix< DDSMat > tLocalOverlappingAdofIds = mDofMgn->get_local_overlapping_adof_ids();
             return tLocalOverlappingAdofIds;
         };

         // ----------------------------------------------------------------------------------------------
         // element dofs
         moris::uint get_num_element_dof()
         {
             return 0;
         };

         // ----------------------------------------------------------------------------------------------
         // number of elements on proc
         moris::uint get_num_my_elements()
         {
             moris::uint tNumEquationObj= mMSI->get_num_eqn_objs();
             return tNumEquationObj;
         };

         // ----------------------------------------------------------------------------------------------
         void get_element_matrix( const moris::uint      & aMyElementInd,
                                        Matrix< DDRMat > & aElementMatrix )
         {
             mMSI->get_equation_obj_jacobian( aMyElementInd, aElementMatrix );
         };

         // ----------------------------------------------------------------------------------------------
         void  get_element_topology( const moris::uint      & aMyElementInd,
                                           Matrix< DDSMat > & aElementTopology )
         {
            mMSI->get_equation_obj_dof_ids( aMyElementInd, aElementTopology );
         };

         // ----------------------------------------------------------------------------------------------
         Matrix< DDUMat > get_constr_dof()
         {
             Matrix< DDUMat > tLocalConstrIds;// = mDofMgn->get_full_to_free_constraints();
             return tLocalConstrIds;
         };

         // ----------------------------------------------------------------------------------------------
         void get_element_rhs( const moris::uint      & aMyElementInd,
                                     Matrix< DDRMat > & aElementRHS )
         {
             mMSI->get_equation_obj_residual( aMyElementInd, aElementRHS, mSolutionVector );
         };
    };
    }
}



#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
