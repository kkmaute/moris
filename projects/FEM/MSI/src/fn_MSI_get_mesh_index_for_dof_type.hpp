/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_MSI_get_mesh_index_for_dof_type.hpp
 *
 */
#ifndef SRC_fn_MSI_get_mesh_index_for_dof_type
#define SRC_fn_MSI_get_mesh_index_for_dof_type

#include "moris_typedefs.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace MSI
    {
        //------------------------------------------------------------------------------

        inline moris_index
        get_mesh_index_for_dof_type(
                MSI::Dof_Type  aDofType,
                Parameter_List& aMSIParameterList )
        {
            // Note: Make sure to add for each DOF type a default interpolation index to the MSI Parameter list

            if ( aDofType == Dof_Type::TEMP )
            {
                return aMSIParameterList.get< moris::sint >( "TEMP" );
            }
            else if ( aDofType == Dof_Type::P )
            {
                return aMSIParameterList.get< moris::sint >( "P" );
            }
            else if ( aDofType == Dof_Type::RHO )
            {
                return aMSIParameterList.get< moris::sint >( "RHO" );
            }
            else if ( aDofType == Dof_Type::E )
            {
                return aMSIParameterList.get< moris::sint >( "E" );
            }
            else if ( aDofType == Dof_Type::EVP )
            {
                return aMSIParameterList.get< moris::sint >( "EVP" );
            }
            else if ( aDofType == Dof_Type::EVT )
            {
                return aMSIParameterList.get< moris::sint >( "EVT" );
            }

            else if ( aDofType == Dof_Type::UX )
            {
                return aMSIParameterList.get< moris::sint >( "UX" );
            }
            else if ( aDofType == Dof_Type::UY )
            {
                return aMSIParameterList.get< moris::sint >( "UY" );
            }
            else if ( aDofType == Dof_Type::UZ )
            {
                return aMSIParameterList.get< moris::sint >( "UZ" );
            }
            else if ( aDofType == Dof_Type::VX )
            {
                return aMSIParameterList.get< moris::sint >( "VX" );
            }
            else if ( aDofType == Dof_Type::VY )
            {
                return aMSIParameterList.get< moris::sint >( "VY" );
            }
            else if ( aDofType == Dof_Type::VZ )
            {
                return aMSIParameterList.get< moris::sint >( "VZ" );
            }
            else if ( aDofType == Dof_Type::MX )
            {
                return aMSIParameterList.get< moris::sint >( "MX" );
            }
            else if ( aDofType == Dof_Type::MY )
            {
                return aMSIParameterList.get< moris::sint >( "MY" );
            }
            else if ( aDofType == Dof_Type::MZ )
            {
                return aMSIParameterList.get< moris::sint >( "MZ" );
            }
            else if ( aDofType == Dof_Type::EVX )
            {
                return aMSIParameterList.get< moris::sint >( "EVX" );
            }
            else if ( aDofType == Dof_Type::EVY )
            {
                return aMSIParameterList.get< moris::sint >( "EVY" );
            }
            else if ( aDofType == Dof_Type::EVZ )
            {
                return aMSIParameterList.get< moris::sint >( "EVZ" );
            }
            else if ( aDofType == Dof_Type::NLSX )
            {
                return aMSIParameterList.get< moris::sint >( "NLSX" );
            }
            else if ( aDofType == Dof_Type::NLSY )
            {
                return aMSIParameterList.get< moris::sint >( "NLSY" );
            }
            else if ( aDofType == Dof_Type::NLSZ )
            {
                return aMSIParameterList.get< moris::sint >( "NLSZ" );
            }

            else if ( aDofType == Dof_Type::L2 )
            {
                return aMSIParameterList.get< moris::sint >( "L2" );
            }
            else if ( aDofType == Dof_Type::MAPPING_DOF )
            {
                return aMSIParameterList.get< moris::sint >( "MAPPING_DOF" );
            }
            else if ( aDofType == Dof_Type::LS1 )
            {
                return aMSIParameterList.get< moris::sint >( "LS1" );
            }
            else if ( aDofType == Dof_Type::LS2 )
            {
                return aMSIParameterList.get< moris::sint >( "LS2" );
            }

            else if ( aDofType == Dof_Type::THETA )
            {
                return aMSIParameterList.get< moris::sint >( "THETA" );
            }
            else if ( aDofType == Dof_Type::PHID )
            {
                return aMSIParameterList.get< moris::sint >( "PHID" );
            }
            else if ( aDofType == Dof_Type::PHISD )
            {
                return aMSIParameterList.get< moris::sint >( "PHISD" );
            }

            else if ( aDofType == Dof_Type::VISCOSITY )
            {
                return aMSIParameterList.get< moris::sint >( "VISCOSITY" );
            }
            else if ( aDofType == Dof_Type::STRESS_DOF )
            {
                return aMSIParameterList.get< moris::sint >( "STRESS_DOF" );
            }

            else if ( aDofType == Dof_Type::NLEQSTRAIN )
            {
                return aMSIParameterList.get< moris::sint >( "NLEQSTRAIN" );
            }
            else if ( aDofType == Dof_Type::HISTORY )
            {
                return aMSIParameterList.get< moris::sint >( "HISTORY" );
            }

            else if ( aDofType == Dof_Type::UNDEFINED )
            {
                return MORIS_INDEX_MAX;
            }

            else
            {
                MORIS_ERROR( false,
                        "get_mesh_index_for_dof_type(): Dof type does not exist in parameter list provided. Check dof type enums" );
                return 0;
            }
        }
    }    // namespace MSI
}    // namespace moris

#endif /* fn_MSI_get_mesh_index_for_dof_type.hpp */
