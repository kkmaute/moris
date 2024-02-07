/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Lagrange_Multiplier_L2.cpp
 *
 */

#include "cl_FEM_SP_Lagrange_Multiplier_L2.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

    	SP_Lagrange_Multiplier_L2::SP_Lagrange_Multiplier_L2()
        {
            // set size for the property pointer cell
            mLeaderProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //------------------------------------------------------------------------------

        Vector< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Leader_Follower > > SP_Lagrange_Multiplier_L2::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Lagrange_Multiplier_L2::eval_SP()
        {
            // get the material property
            const std::shared_ptr< Property > & tPropMaterial =
                    mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) / tPropMaterial->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void SP_Lagrange_Multiplier_L2::eval_dSPdLeaderDOF(
                const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mLeaderGlobalDofTypeMap( tDofType );

            // get FI for derivative dof type
            Field_Interpolator * tFIDer =
                    mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdLeaderDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the material property
            const std::shared_ptr< Property > & tPropMaterial =
                    mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // if material property depends on the dof type
            if ( tPropMaterial->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                MORIS_ERROR( false, "SP_Lagrange_Multiplier_L2::eval_dSPdLeaderDOF() - "
                		"This stabilization parameter doesn't have a derivative implemented. Implement this if you want to use a non-constant material parameter." );
            }
            else
            {
                mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

