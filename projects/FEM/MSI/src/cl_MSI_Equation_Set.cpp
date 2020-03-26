/*
 * cl_Equation_Object.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"

#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace MSI
    {

//------------------------------------------------------------------------------
        moris::Cell < enum MSI::Dof_Type > Equation_Set::get_requested_dof_types()
        {
            return mModelSolverInterface->get_solver_interface()->get_requested_dof_types();
        }

//------------------------------------------------------------------------------
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > Equation_Set::get_secundary_dof_types()
        {
            return mModelSolverInterface->get_solver_interface()->get_secundary_dof_types();
        }

//------------------------------------------------------------------------------
        moris::Cell < enum GEN_DV > Equation_Set::get_requested_dv_types()
        {
            moris::Cell< enum GEN_DV > tDvTypes;
            mDesignVariableInterface->get_ip_requested_dv_types( tDvTypes );
            return tDvTypes;
        }

//------------------------------------------------------------------------------
        moris::Cell< moris::Cell < moris_index > > & Equation_Set::get_QI_assembly_map()
        {
            return mDesignVariableInterface->get_QI_assembly_map();
        }

//-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */
