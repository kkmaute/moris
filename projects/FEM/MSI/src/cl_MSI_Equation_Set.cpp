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
        void Equation_Set::set_requested_IQI_types( const moris::Cell< moris::Cell< enum fem::IQI_Type > > & aRequestedIQITypes )
        {
            mRequestedIQITypes = aRequestedIQITypes;
        }

//------------------------------------------------------------------------------
        const moris::Cell< moris::Cell< enum fem::IQI_Type > > & Equation_Set::get_requested_IQI_types()
        {
            return mRequestedIQITypes;
        }

//------------------------------------------------------------------------------
        void Equation_Set::create_requested_IQI_type_map()
        {
            mRequestedIQITypeAssemblyMap.resize( mRequestedIQITypes.size() );

            for( uint Ik = 0; Ik < mRequestedIQITypes.size(); Ik++ )
            {
                mRequestedIQITypeAssemblyMap( Ik ).resize( static_cast< sint >( fem::IQI_Type::END_IQI_TYPE ), gNoIndex );
            }

            uint tCounter =0;
            for( uint Ik = 0; Ik < mRequestedIQITypes.size(); Ik++ )
            {
                for( uint Ii = 0; Ii < mRequestedIQITypes( Ik ).size(); Ii++ )
                {
                    mRequestedIQITypeAssemblyMap( Ik )( static_cast< sint >( mRequestedIQITypes( Ik )( Ii ) ) ) = tCounter++;
                }
            }
        }

//------------------------------------------------------------------------------
        moris::Cell < enum GEN_DV > Equation_Set::get_requested_dv_types()
        {
            moris::Cell< enum GEN_DV > tDvTypes;
            mDesignVariableInterface->get_ip_requested_dv_types( tDvTypes );
            return tDvTypes;
        }

//------------------------------------------------------------------------------
        // FIXME delete this one and build from requested IQI types
        moris_index Equation_Set::get_QI_assembly_index( const enum Phase_Type    aPhaseType,
                                                         const enum fem::IQI_Type aIQIType )
        {
            return mRequestedIQITypeAssemblyMap( static_cast< uint >( aPhaseType ) )( static_cast< uint >( aIQIType ) );
        }

//-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */
