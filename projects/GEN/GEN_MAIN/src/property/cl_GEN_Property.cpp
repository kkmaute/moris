/*
 * cl_GEN_Property.cpp
 *
 *  Created on: Dec 18, 2019
 *      Author: sonne
 */

#include "cl_GEN_Property.hpp"

namespace moris
{
namespace ge
{
//------------------------------------------------------------------------------
void GEN_Property::build_dv_type_map()
{
    // get number of dv types
    uint tNumDvTypes = mDvTypes.size();

    // determine the max Dv_Type enum
    sint tMaxEnum = 0;
    for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
    {
        tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDvTypes( iDV )( 0 ) ) );
    }
    tMaxEnum++;

    // set the Dv_Type map size
    mDvTypeMap.set_size( tMaxEnum, 1, -1 );

    // fill the Dv_Type map
    for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
    {
        // fill the property map
        mDvTypeMap( static_cast< int >( mDvTypes( iDV )( 0 ) ), 0 ) = iDV;
    }
}

//------------------------------------------------------------------------------
bool GEN_Property::check_dv_dependency( const moris::Cell< PDV_Type > aDvType )
{
    // set bool for dependency
    bool tDvDependency = false;

    // get index for dv type
    uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

    // if aDvType is an active dv type for the property
    if( tDvIndex < mDvTypeMap.numel() && mDvTypeMap( tDvIndex ) != -1 )
    {
        // bool is set to true
        tDvDependency = true;
    }
    // return bool for dependency
    return tDvDependency;
}

//------------------------------------------------------------------------------
void GEN_Property::set_dv_field_interpolators( moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolators )
{
    // check size
    MORIS_ASSERT( aFieldInterpolators.size() == mDvTypes.size(),
            "Property::set_dv_field_interpolators - wrong input size. " );

    // check field interpolator type
    bool tCheckFI = true;
    for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
    {
        tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == mDvTypes( iFI )( 0 ) );
    }
    MORIS_ASSERT( tCheckFI, "Property::set_dv_field_interpolators - wrong field interpolator dv type. ");

    // set field interpolators
    mDvFI = aFieldInterpolators;
}

//------------------------------------------------------------------------------
void GEN_Property::check_dv_field_interpolators()
{
    // get number of dv types
    uint tNumDvTypes = mDvTypes.size();

    // check field interpolators cell size
    MORIS_ASSERT( mDvFI.size() == tNumDvTypes, "Property::check_dv_field_interpolators - wrong FI size. " );

    // loop over the field interpolator pointers
    for( uint iFI = 0; iFI < tNumDvTypes; iFI++ )
    {
        // check that the field interpolator was set
        MORIS_ASSERT( mDvFI( iFI ) != nullptr, "Property::check_dv_field_interpolators - FI missing. " );
    }
}

//------------------------------------------------------------------------------
const Matrix< DDRMat > & GEN_Property::val()
{
    // if the property was not evaluated
    if( mPropEval )
    {
        // evaluate the property
        this->eval_Prop();

        // set bool for evaluation
        mPropEval = false;
    }
    // return the property value
    return mProp;
}
//------------------------------------------------------------------------------
void GEN_Property::eval_Prop()
{
    // check that mValFunc was assigned
    MORIS_ASSERT( mValFunction != nullptr, "Property::eval_Prop - mValFunction not assigned. " );

    this->check_dv_field_interpolators();

    // use mValFunction to evaluate the property
//    mProp = mValFunction( mParameters, mDvFI );
    mProp = mValFunction( mParameters );
}
//------------------------------------------------------------------------------
const Matrix< DDRMat > & GEN_Property::dPropdDV( const moris::Cell< PDV_Type > aDvType )
{
    // if aDvType is not an active dv type for the property
    MORIS_ERROR( this->check_dv_dependency( aDvType ), "Property::dPropdDV - no dependency in this dv type." );

    // get dv type index
    uint tDvIndex = mDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

    // if the derivative has not been evaluated yet
    if( mPropDvDerEval( tDvIndex ) )
    {
        // evaluate the derivative
        this->eval_dPropdDV( aDvType );

        // set bool for evaluation
        mPropDvDerEval( tDvIndex ) = false;
    }

    // return the derivative
    return mPropDvDer( tDvIndex );
}
//------------------------------------------------------------------------------
void GEN_Property::eval_dPropdDV( const moris::Cell< PDV_Type > aDvType )
{
    // get the dv index
    uint tDvIndex = mDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

    // check that mDofDerFunctions was assigned
    MORIS_ASSERT( mDvDerFunctions( tDvIndex ) != nullptr, "Property::eval_dPropdDV - mDvDerFunctions not assigned. " );

    // check that mGeometryInterpolator was assigned
    MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_dPropdDV - mGeometryInterpolator not assigned. " );

    // check that the field interpolators are set
    this->check_dv_field_interpolators();

    // if so use mDerivativeFunction to compute the derivative
//    mPropDvDer( tDvIndex ) = mDvDerFunctions( tDvIndex )( mParameters, mDvFI );
    mPropDvDer( tDvIndex ) = mDvDerFunctions( tDvIndex )( mParameters );
}

//------------------------------------------------------------------------------

}   // end ge namespace
}   // end moris namespace


