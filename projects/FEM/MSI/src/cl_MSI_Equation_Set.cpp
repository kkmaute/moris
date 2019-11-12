/*
 * cl_Equation_Object.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"

#include "cl_Vector.hpp"

namespace moris
{
    namespace MSI
    {

    void Equation_Set::create_dof_assembly_map_EQ()
    {
        // get List of requested dof types
        Cell < enum MSI::Dof_Type >  tRequestedDofTypes =  mModelSolverInterface->get_solver_interface()
                                                                                ->get_requested_dof_types();

        uint tCounter = 0;

        // loop over requested dof types and check for contribution to jacobian. There will be no contribution if its a second or higher dof in an field interpolator
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = mMasterDofTypeMap( static_cast< int >( tRequestedDofTypes( Ik ) ) );

            if( tDofIndex != -1 )
            {
                mRequestedTypeToIndexMap( 0 )[ tRequestedDofTypes( Ik ) ] = Ik;

                tCounter++;
            }

            tDofIndex = mSlaveDofTypeMap( static_cast< int >( tRequestedDofTypes( Ik ) ) );

            if( tDofIndex != -1 )
            {
                mRequestedTypeToIndexMap( 1 )[ tRequestedDofTypes( Ik ) ] = Ik;

                tCounter++;
            }
        }

        mDofAssemblyMap_2.set_size( tCounter, 2, -1 );

        tCounter = 0;

        uint tCounter_2 = 0;

        // loop over requested dof types
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = mMasterDofTypeMap( static_cast< int >( tRequestedDofTypes( Ik ) ) );

            if( tDofIndex != -1 )
            {
                mDofAssemblyMap_2( Ik, 0 ) = tCounter;
                mDofAssemblyMap_2( Ik, 1 ) = tCounter + mNumDofMap( tDofIndex, 0 ) - 1;

                tCounter += mNumDofMap( tDofIndex, 0 );

                tCounter_2++;
            }
        }

        // loop over requested dof types
        for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
        {
            sint tDofIndex = mSlaveDofTypeMap( static_cast< int >( tRequestedDofTypes( Ik ) ) );

            if( tDofIndex != -1 )
            {
                mDofAssemblyMap_2( tCounter_2 + Ik, 0 ) = tCounter;
                mDofAssemblyMap_2( tCounter_2 + Ik, 1 ) = tCounter + mNumDofMap( tDofIndex, 1 ) - 1;

                tCounter += mNumDofMap( tDofIndex, 1 );
            }
        }
    }

//-------------------------------------------------------------------------------------------------

}
}
