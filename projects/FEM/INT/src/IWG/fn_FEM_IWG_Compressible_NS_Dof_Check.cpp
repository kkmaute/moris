/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Compressible_NS_Dof_Check.cpp
 *
 */

#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        bool check_residual_dof_types(
                const Vector< Vector< MSI::Dof_Type > > & aResidualDofTypes )
        {
            // initialize
            bool tCheck = true;

            // FIXME: only density and pressure primitive variables supported for now

            // check that there are exactly 3 residual DoF types
            tCheck = tCheck && ( aResidualDofTypes.size() == 3 );
            if( aResidualDofTypes.size() != 3 )
            {
                MORIS_LOG_ERROR( "Compressible_NS_Dof_Check::check_residual_dof_types() - List of Residual DoF types must be of length 3 for pressure primitive variables." );
            }

            // check that the correct DoF types are present
            tCheck = tCheck && ( ( aResidualDofTypes( 0 )( 0 ) == MSI::Dof_Type::RHO ) || ( aResidualDofTypes( 0 )( 0 ) == MSI::Dof_Type::P ) );
            if ( !( ( aResidualDofTypes( 0 )( 0 ) == MSI::Dof_Type::RHO ) || ( aResidualDofTypes( 0 )( 0 ) == MSI::Dof_Type::P ) ) )
            {
                MORIS_LOG_ERROR( "Compressible_NS_Dof_Check::check_residual_dof_types() - First DoF type must be density or pressure." );
            }

            tCheck = tCheck && ( aResidualDofTypes( 1 )( 0 ) == MSI::Dof_Type::VX );
            tCheck = tCheck && ( aResidualDofTypes( 2 )( 0 ) == MSI::Dof_Type::TEMP );
            if ( !( ( aResidualDofTypes( 1 )( 0 ) == MSI::Dof_Type::VX ) and ( aResidualDofTypes( 2 )( 0 ) == MSI::Dof_Type::TEMP ) ) )
            {
                MORIS_LOG_ERROR( "Compressible_NS_Dof_Check::check_residual_dof_types() - Second and third DoF types must be velocity and temperature." );
            }

            // return whether check was successful or not
            return tCheck;
        }

        //------------------------------------------------------------------------------

        bool check_dof_dependencies(
                fem::Set                                          * aFemSet,
                const Vector< Vector< MSI::Dof_Type > > & aResidualDofTypes,
                const Vector< Vector< MSI::Dof_Type > > & aRequestedDofTypeList )
        {
            // get and check the number of dof dependencies
            uint tNumDofDependencies = aRequestedDofTypeList.size();
            if ( tNumDofDependencies != 3 )
            {
                MORIS_LOG_ERROR( "Compressible_NS_Dof_Check::check_dof_dependencies() - More than three DoF dependencies." );
            }

            // check dof dependencies and ordering
            bool tCheck = ( aRequestedDofTypeList( 0 )( 0 ) == MSI::Dof_Type::RHO ) || ( aRequestedDofTypeList( 0 )( 0 ) == MSI::Dof_Type::P );
            tCheck = tCheck && ( aRequestedDofTypeList( 1 )( 0 ) == MSI::Dof_Type::VX );
            tCheck = tCheck && ( aRequestedDofTypeList( 2 )( 0 ) == MSI::Dof_Type::TEMP );

            if ( !tCheck )
            {
                MORIS_LOG_ERROR(
                        "Compressible_NS_Dof_Check::check_dof_dependencies() - Only Pressure and Density Primitive variables supported for now."
                        "DoF ordering must be (rho,u,T) or (p,u,T)." );
            }

            // check that the assembly map is a connected block
            uint tLeaderDof1Index      = aFemSet->get_dof_index_for_type( aResidualDofTypes( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDof2Index      = aFemSet->get_dof_index_for_type( aResidualDofTypes( 1 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDof3Index      = aFemSet->get_dof_index_for_type( aResidualDofTypes( 2 )( 0 ), mtk::Leader_Follower::LEADER );
            sint tDofDep1Index         = aFemSet->get_dof_index_for_type( aRequestedDofTypeList( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            sint tDofDep2Index         = aFemSet->get_dof_index_for_type( aRequestedDofTypeList( 1 )( 0 ), mtk::Leader_Follower::LEADER );
            sint tDofDep3Index         = aFemSet->get_dof_index_for_type( aRequestedDofTypeList( 2 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDep1StopIndex  = aFemSet->get_jac_dof_assembly_map()( tLeaderDof1Index )( tDofDep1Index, 1 );
            uint tLeaderDep2StartIndex = aFemSet->get_jac_dof_assembly_map()( tLeaderDof2Index )( tDofDep2Index, 0 );
            uint tLeaderDep2StopIndex  = aFemSet->get_jac_dof_assembly_map()( tLeaderDof2Index )( tDofDep2Index, 1 );
            uint tLeaderDep3StartIndex = aFemSet->get_jac_dof_assembly_map()( tLeaderDof3Index )( tDofDep3Index, 0 );

            if ( !( ( tLeaderDep1StopIndex + 1 == tLeaderDep2StartIndex ) && ( tLeaderDep2StopIndex + 1 == tLeaderDep3StartIndex ) ) )
            {
                MORIS_LOG_ERROR( "Compressible_NS_Dof_Check::check_dof_dependencies() - Assembly map is not connected." );
            }

            return tCheck;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

