/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Ghost_Displacement.cpp
 *
 */

#include "cl_FEM_SP_Ghost_Displacement.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    SP_Ghost_Displacement::SP_Ghost_Displacement()
    {
        mHasFollower = true;

        // set size for the property pointer cells
        mLeaderProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
        mFollowerProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
    }

    //------------------------------------------------------------------------------

    Vector< std::tuple<
            fem::Measure_Type,
            mtk::Primary_Void,
            mtk::Leader_Follower > >
    SP_Ghost_Displacement::get_cluster_measure_tuple_list()
    {
        return { mElementSizeTuple };
    }

    //------------------------------------------------------------------------------

    void
    SP_Ghost_Displacement::eval_SP()
    {
        // get element size cluster measure value
        real tElementSize =    //
                mCluster->get_cluster_measure(
                                std::get< 0 >( mElementSizeTuple ),
                                std::get< 1 >( mElementSizeTuple ),
                                std::get< 2 >( mElementSizeTuple ) )
                        ->val()( 0 );

        // compute stabilization parameter value
        mPPVal = mParameters( 0 )                                                                     //
               * std::pow( tElementSize, 2.0 * ( mOrder - ( this->*mGetWeakFormOrder )() ) + 1.0 )    //
               * mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );
    }

    //------------------------------------------------------------------------------

    void
    SP_Ghost_Displacement::eval_dSPdLeaderDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get element size cluster measure value
        real tElementSize =    //
                mCluster->get_cluster_measure(
                                std::get< 0 >( mElementSizeTuple ),
                                std::get< 1 >( mElementSizeTuple ),
                                std::get< 2 >( mElementSizeTuple ) )
                        ->val()( 0 );

        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mLeaderGlobalDofTypeMap( tDofType );

        // get FI for derivative dof type
        Field_Interpolator* tFIDer =
                mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // reset the matrix
        mdPPdLeaderDof( tDofIndex ).set_size(    //
                1,                               //
                tFIDer->get_number_of_space_time_coefficients() );

        // get the material property
        const std::shared_ptr< Property >& tPropMaterial =
                mLeaderProp( static_cast< uint >( SP_Property_Type::MATERIAL ) );

        if ( tPropMaterial->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mdPPdLeaderDof( tDofIndex ) =                                                                                   //
                    mParameters( 0 ) * std::pow( tElementSize, 2.0 * ( mOrder - ( this->*mGetWeakFormOrder )() ) + 1.0 )    //
                    * tPropMaterial->dPropdDOF( aDofTypes );
        }
        else
        {
            mdPPdLeaderDof( tDofIndex ).fill( 0.0 );
        }
    }

    //------------------------------------------------------------------------------

    real
    SP_Ghost_Displacement::get_weak_form_order_init()
    {

        if ( mParameters.size() == 2 )
        {
            // if the second parameter has been set, use it
            mWeakFormOrder = mParameters( 1 )( 0 );
        }
        else if ( mParameters.size() == 1 )
        {
            // default is 1.
            mWeakFormOrder = 1;
        }
        else
        {
            MORIS_ERROR( false,
                    "SP_Ghost_Displacement::get_weak_form_order_init() - incorrect number of sp function_parameters" );
        }

        // redirect the function pointer for all other points
        mGetWeakFormOrder = &SP_Ghost_Displacement::get_weak_form_order;

        return mWeakFormOrder;
    }

    //------------------------------------------------------------------------------

    real
    SP_Ghost_Displacement::get_weak_form_order()
    {
        return mWeakFormOrder;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
