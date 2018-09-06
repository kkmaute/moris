/*
 * cl_MSI_Solver_Interface.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_

#include "linalg.hpp"

#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_Solver_Input.hpp"

namespace moris
{
    namespace MSI
    {
    class MSI_Solver_Interface : public moris::Solver_Input
    {
    private:
        moris::MSI::Model_Solver_Interface*   mMSI;
        moris::MSI::Dof_Manager*   mDofMgn;

    public:
        MSI_Solver_Interface( )
        {};

        MSI_Solver_Interface( moris::MSI::Model_Solver_Interface * aMSI,
                              moris::MSI::Dof_Manager            * aDofMgn ) : mMSI( aMSI ),
                                                                               mDofMgn( aDofMgn )
        {};

        ~MSI_Solver_Interface()
        {};

        // ----------------------------------------------------------------------------------------------
         // local dimension of the problem
         moris::uint get_num_my_dofs()
         {
             moris::uint tNumDofs = mDofMgn->get_num_adofs();
             return tNumDofs;
         };

         // ----------------------------------------------------------------------------------------------
         // local-to-global map
         Mat < int > get_my_local_global_map()
         {
             moris::Mat< int> tLocalAdofIds = mDofMgn->get_local_adof_ids();
             return tLocalAdofIds;
         };

         moris::Mat < int > get_my_local_global_overlapping_map( )
         {
             moris::Mat< int> tLocalOverlappingAdofIds = mDofMgn->get_local_overlapping_adof_ids();
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
         void get_element_matrix( const moris::uint               & aMyElementInd,
                                        moris::Mat< moris::real > & aElementMatrix )
         {
             mMSI->get_equation_obj_jacobian( aMyElementInd, aElementMatrix );
         };

         // ----------------------------------------------------------------------------------------------
         void  get_element_topology( const moris::uint       & aMyElementInd,
                                           moris::Mat< int > & aElementTopology )
         {
            mMSI->get_equation_obj_dof_ids( aMyElementInd, aElementTopology );
         };

         // ----------------------------------------------------------------------------------------------
         moris::Mat< moris::uint > get_constr_dof()
         {
             moris::Mat< moris::uint> tLocalConstrIds;// = mDofMgn->get_full_to_free_constraints();
              return tLocalConstrIds;
         };

         // ----------------------------------------------------------------------------------------------
         void get_element_rhs( const moris::uint               & aMyElementInd,
                                     moris::Mat< moris::real > & aElementRHS )
         {
             mMSI->get_equation_obj_residual( aMyElementInd, aElementRHS );
         };

         // ----------------------------------------------------------------------------------------------
         const moris::Mat< moris::uint > get_adof_ind_map()
         {
              return mDofMgn->get_adof_ind_map();
         };
    };
    }
}



#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
