/*
 * cl_MSI_Design_Variable_Interface.cpp
 *
 *  Created on: Jan 10, 2020
 *      Author: schmidt
 */
#include "cl_MSI_Design_Variable_Interface_Proxy.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_MSI_Equation_Model.hpp" 
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MDL_Model.hpp"  

namespace moris
{
    namespace MSI
    {

    //---------------------------------------------------------------------------------------------------
    void Design_Variable_Interface_Proxy::set_requested_IQIs( const moris::Cell< std::string> & aRequestedIQINames )
    {
        mModel->get_fem_model()->set_requested_IQI_names(aRequestedIQINames);
    }

    //---------------------------------------------------------------------------------------------------

    }
}
