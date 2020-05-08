/*
 * cl_Equation_Model.cpp
 *
 *  Created on: May 08, 2020
 *      Author: schmidt
 */

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace MSI
    {

//------------------------------------------------------------------------------
        void Equation_Model::compute_dQIdp()
        {
            // create map object
            moris::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );
            mdQiduMap = tMatFactory.create_map( mDesignVariableInterface->get_my_local_global_map() );

            uint tNumRHMS = this->get_num_rhs();

            // full vector and prev full vector
            mdQidu = tMatFactory.create_vector( mdQiduMap, tNumRHMS );

            mdQidu->vec_put_scalar( 0.0 );

            // Get local number of elements
            moris::uint tNumSets = mFemSets.size();

            // Loop over all local elements to build matrix graph
            for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
            {
                moris::uint tNumEquationObjectOnSet = mFemSets( Ii )->get_num_equation_objects();

                //this->initialize_set( Ii, true );

                for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
                {
                    mFemSets( Ii )->get_equation_object_list()( Ik )->compute_dQIdp();
                }

                //this->free_block_memory( Ii );
            }
        }

//-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */
