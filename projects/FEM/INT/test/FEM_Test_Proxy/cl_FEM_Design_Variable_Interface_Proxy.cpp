/*
 * cl_MSI_Design_Variable_Interface.cpp
 *
 *  Created on: Jan 10, 20120
 *      Author: schmidt
 */
#include "cl_FEM_Design_Variable_Interface_Proxy.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_MSI_Equation_Model.hpp" 
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MDL_Model.hpp"  

namespace moris
{
    namespace fem
    {

//-------------------------------------------------------------------------------------------------------
    void FEM_Design_Variable_Interface_Proxy::set_requested_IQI_type( const moris::Cell< moris::Cell< enum fem::IQI_Type > > & aRequestedIQIType )
    {
	    uint tNumEquationSets = mModel->get_fem_model()->get_equation_sets().size();
		
		for( uint Ik = 0; Ik <tNumEquationSets; Ik++ )
		{
			mModel->get_fem_model()->get_equation_sets()( Ik )->set_requested_IQI_types( aRequestedIQIType );
		}
    }

//-------------------------------------------------------------------------------------------------------

    }
}
