/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Property.cpp
 *
 */

#include "cl_FEM_Property.hpp"    //FEM/MSI/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    Property::Property() {}

    //------------------------------------------------------------------------------

    void
    Property::reset_eval_flags()
    {
        // reset the property value
        mPropEval = true;

        // reset the property derivatives wrt x
        mPropSpaceDerEval.fill( true );

        // reset the property derivatives wrt dof
        mPropDofDerEval.fill( true );

        // reset the property derivatives wrt x and dof
        mPropSpaceDofDerEval.fill( true );

        // reset the property derivatives wrt dv
        mPropDvDerEval.fill( true );
    }

    //------------------------------------------------------------------------------

    void
    Property::set_field_interpolator_manager(
            Field_Interpolator_Manager* aFieldInterpolatorManager )
    {
        // set field interpolator manager
        mFIManager = aFieldInterpolatorManager;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes )
    {
        // set dof type list
        mDofTypes = aDofTypes;

        // build dof type map
        this->build_dof_type_map();

        // number of dof types
        uint tNumDofTypes = mDofTypes.size();

        // set mDofDerFunctions size
        mDofDerFunctions.assign( tNumDofTypes, nullptr );

        // set mPropDofDerEval size
        mPropDofDerEval.set_size( tNumDofTypes, 1, true );

        // set mPropDofDer size
        mPropDofDer.resize( tNumDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    Property::build_dof_type_map()
    {
        // determine the max Dof_Type enum
        sint tMaxEnum = 0;
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDof )( 0 ) ) );
        }

        tMaxEnum++;

        // set the Dof_Type map size
        mDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Dof_Type map
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            // fill the property map
            mDofTypeMap( static_cast< int >( mDofTypes( iDof )( 0 ) ), 0 ) = iDof;
        }
    }

    //------------------------------------------------------------------------------

    bool
    Property::check_dof_dependency(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // get the dof type index
        uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

        // if aDofType is an active dof type for the property
        if ( tDofIndex >= mDofTypeMap.numel() )
        {
            return false;
        }

        if ( mDofTypeMap( tDofIndex ) != -1 )
        {
            return true;
        }

        return false;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_dv_type_list(
            const Vector< Vector< gen::PDV_Type > >& aDvTypes )
    {
        // set dv type list
        mDvTypes = aDvTypes;

        // build a dv type map
        this->build_dv_type_map();

        // number of dv types
        uint tNumDvTypes = mDvTypes.size();

        // set mDvDerFunctions size
        mDvDerFunctions.resize( tNumDvTypes, nullptr );

        // set mPropDvDerEval size
        mPropDvDerEval.set_size( tNumDvTypes, 1, true );

        // set mPropdvDer size
        mPropDvDer.resize( tNumDvTypes );
    }

    //------------------------------------------------------------------------------

    void
    Property::build_dv_type_map()
    {
        // get number of dv types
        uint tNumDvTypes = mDvTypes.size();

        // determine the max Dv_Type enum
        sint tMaxEnum = 0;

        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDvTypes( iDV )( 0 ) ) );
        }

        tMaxEnum++;

        // set the Dv_Type map size
        mDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Dof_Type map
        for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
        {
            // fill the property map
            mDvTypeMap( static_cast< int >( mDvTypes( iDV )( 0 ) ), 0 ) = iDV;
        }
    }

    //------------------------------------------------------------------------------

    bool
    Property::check_dv_dependency(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // set bool for dependency
        bool tDvDependency = false;

        // get index for dv type
        uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

        // if aDvType is an active dv type for the property
        if ( tDvIndex < mDvTypeMap.numel() && mDvTypeMap( tDvIndex ) != -1 )
        {
            // bool is set to true
            tDvDependency = true;
        }

        // return bool for dependency
        return tDvDependency;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_field_type_list(
            const Vector< Vector< mtk::Field_Type > >& aFieldTypes )
    {
        // set field type list
        mFieldTypes = aFieldTypes;

        // build a field type map
        this->build_field_type_map();
    }

    //------------------------------------------------------------------------------

    void
    Property::build_field_type_map()
    {
        // get number of field types
        uint tNumFieldTypes = mFieldTypes.size();

        // determine the max Field_Type enum
        sint tMaxEnum = 0;

        for ( uint iFi = 0; iFi < tNumFieldTypes; iFi++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFieldTypes( iFi )( 0 ) ) );
        }

        tMaxEnum++;

        // set the Field_Type map size
        mFieldTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Field_Type map
        for ( uint iFi = 0; iFi < tNumFieldTypes; iFi++ )
        {
            // fill the property map
            mFieldTypeMap( static_cast< int >( mFieldTypes( iFi )( 0 ) ), 0 ) = iFi;
        }
    }

    //------------------------------------------------------------------------------

    bool
    Property::check_field_dependency(
            const Vector< mtk::Field_Type >& aFieldType )
    {
        // set bool for dependency
        bool tFieldDependency = false;

        // get index for field type
        uint tFieldIndex = static_cast< uint >( aFieldType( 0 ) );

        // if aFieldType is an active field type for the property
        if ( tFieldIndex < mFieldTypeMap.numel() && mFieldTypeMap( tFieldIndex ) != -1 )
        {
            // bool is set to true
            tFieldDependency = true;
        }

        // return bool for dependency
        return tFieldDependency;
    }

    //------------------------------------------------------------------------------

    void
    Property::get_non_unique_dof_types(
            Vector< MSI::Dof_Type >& aDofTypes )
    {
        // init counter
        uint tCounter = 0;

        // loop over dof types
        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
        {
            // update counter
            tCounter += mDofTypes( iDOF ).size();
        }

        // reserve memory for dof type list
        aDofTypes.reserve( tCounter );

        // loop over dof types
        for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
        {
            // populate the dof type list
            aDofTypes.append( mDofTypes( iDOF ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    Property::get_non_unique_dof_dv_and_field_types(
            Vector< MSI::Dof_Type >&   aDofTypes,
            Vector< gen::PDV_Type >&   aDvTypes,
            Vector< mtk::Field_Type >& aFieldTypes )
    {
        // init counter
        uint tDofCounter   = 0;
        uint tDvCounter    = 0;
        uint tFieldCounter = 0;

        // loop over dof types
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            // update counter
            tDofCounter += mDofTypes( iDof ).size();
        }

        // loop over dv types
        for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
        {
            // update counter
            tDvCounter += mDvTypes( iDv ).size();
        }

        // loop over field types
        for ( uint iFi = 0; iFi < mFieldTypes.size(); iFi++ )
        {
            // update counter
            tFieldCounter += mFieldTypes( iFi ).size();
        }

        // reserve memory for dof, dv and field type lists
        aDofTypes.reserve( tDofCounter );
        aDvTypes.reserve( tDvCounter );
        aFieldTypes.reserve( tFieldCounter );

        // loop over dof types
        for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
        {
            // populate the dof type list
            aDofTypes.append( mDofTypes( iDof ) );
        }

        // loop over dv types
        for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
        {
            // populate the dv type list
            aDvTypes.append( mDvTypes( iDv ) );
        }

        // loop over field types
        for ( uint iFi = 0; iFi < mFieldTypes.size(); iFi++ )
        {
            // populate the dv type list
            aFieldTypes.append( mFieldTypes( iFi ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    Property::set_val_function( PropertyFunc aValFunction )
    {
        // set the value function
        mValFunction = std::move( aValFunction );

        // set setting flag
        mSetValFunction = true;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_space_derivative_functions(
            const Vector< PropertyFunc >& aSpaceDerFunctions )
    {
        // get the max order of space derivatives
        uint tMaxOrder = aSpaceDerFunctions.size();

        // for now allowed max order is 1
        MORIS_ASSERT( tMaxOrder <= mMaxSpaceDerOrder,           //
                "Property::set_space_derivative_functions()"    //
                " Higher order property space derivatives not supported." );

        // if there is a space dependency
        if ( tMaxOrder > 0 )
        {
            // set space dependency to true
            mSpaceDependency = true;
        }

        // set functions for derivatives wrt space
        mSpaceDerFunctions = aSpaceDerFunctions;

        // set storage
        mPropSpaceDer.resize( tMaxOrder );

        // set evaluation flag
        mPropSpaceDerEval.set_size( tMaxOrder, 1, true );

        // set setting flag
        mSetSpaceDerFunctions = true;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_dof_derivative_functions(
            const Vector< PropertyFunc >& aDofDerFunctions )
    {
        // set functions for derivatives wrt dof
        mDofDerFunctions = aDofDerFunctions;

        // set setting flag
        mSetDofDerFunctions = true;
    }

    //------------------------------------------------------------------------------

    void
    Property::set_space_dof_derivative_functions(
            const Vector< PropertyFunc >& aSpaceDofDerFunctions )
    {
        // get the max order of space derivatives
        // FIXME aSpaceDofDerFunctions should be a vector of vector
        // with space order as row and dof as column
        // right now limited to order 1
        uint tMaxOrder = 0;
        if ( aSpaceDofDerFunctions.size() > 0 )
        {
            tMaxOrder = 1;
        }

        // for now allowed max order is 1
        MORIS_ASSERT( tMaxOrder <= mMaxSpaceDerOrder,             //
                "Property - set_space_derivative_functions():"    //
                " Higher order property space derivatives not supported." );

        // number of dof types
        uint tNumDofTypes = mDofTypes.size();

        // set mDofDerFunctions size
        mSpaceDofDerFunctions.assign( tNumDofTypes, nullptr );

        // set functions for derivatives wrt space and wrt dof
        mSpaceDofDerFunctions = aSpaceDofDerFunctions;
        
        // set mPropSpaceDofDerEval size
        mPropSpaceDofDerEval.set_size( tMaxOrder, tNumDofTypes, true );

        // set mPropSpaceDofDer size
        mPropSpaceDofDer.resize( tMaxOrder );
        for ( uint iOrder = 0; iOrder < tMaxOrder; iOrder++ )
        {
            for ( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
            {
                mPropSpaceDofDer( iOrder ).resize( tNumDofTypes );
            }
        }

        // set setting flag
        mSetSpaceDofDerFunctions = true;
    }

    //------------------------------------------------------------------------------

    bool
    Property::check_space_dependency()
    {
        // return bool for dependency on space
        return mSpaceDependency;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Property::val()
    {
        // if the property was not evaluated
        if ( mPropEval )
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

    void
    Property::eval_Prop()
    {
        // check that mFIManager was assigned
        MORIS_ASSERT( mFIManager != nullptr,
                "Property::eval_Prop - mFIManager not assigned. " );

        // use mValFunction to evaluate the property
        mValFunction( mProp, mParameters, mFIManager );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Property::dPropdx( const uint& aOrder )
    {
        // if the property was not evaluated
        if ( mPropSpaceDerEval( aOrder - 1 ) )
        {
            // evaluate the property
            this->eval_dPropdx( aOrder );

            // set bool for evaluation
            mPropSpaceDerEval( aOrder - 1 ) = false;
        }

        // return the property value
        return mPropSpaceDer( aOrder - 1 );
    }

    //------------------------------------------------------------------------------

    void
    Property::eval_dPropdx( const uint& aOrder )
    {
        // check that mSetSpaceDerFunc was assigned
        MORIS_ASSERT( mSetSpaceDerFunctions == true,
                "Property::eval_dPropdx - mSpaceDerFunction not assigned for order. " );
        MORIS_ASSERT( mSpaceDerFunctions( aOrder - 1 ) != nullptr,
                "Property::eval_dPropdx - mSpaceDerFunctions not assigned for order. " );

        // check that mFIManager was assigned
        MORIS_ASSERT( mFIManager != nullptr,
                "Property::eval_dPropdx - mFIManager not assigned. " );

        // use mSpaceDerFunctions to evaluate the property
        mSpaceDerFunctions( aOrder - 1 )( mPropSpaceDer( aOrder - 1 ), mParameters, mFIManager );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Property::dPropdDOF(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // if aDofType is not an active dof type for the property
        MORIS_ERROR( this->check_dof_dependency( aDofType ),
                "Property::dPropdDOF - no dependency in this dof type." );

        // get the dof index
        uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mPropDofDerEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dPropdDOF( aDofType );

            // set bool for evaluation
            mPropDofDerEval( tDofIndex ) = false;
        }

        // return the derivative
        return mPropDofDer( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    Property::eval_dPropdDOF(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // get the dof index
        uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // check that mDofDerFunctions was assigned
        MORIS_ASSERT( mSetDofDerFunctions == true,
                "Property::eval_dPropdDOF - mDofDerFunctions not assigned. " );
        MORIS_ASSERT( mDofDerFunctions( tDofIndex ) != nullptr,
                "Property::eval_dPropdDOF - mDerFunction not assigned. " );

        // check that mFIManager was assigned
        MORIS_ASSERT( mFIManager != nullptr,
                "Property::eval_Prop - mFIManager not assigned. " );

        // if so use mDerivativeFunction to compute the derivative
        mDofDerFunctions( tDofIndex )( mPropDofDer( tDofIndex ), mParameters, mFIManager );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Property::dPropdxdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const uint&                    aOrder )
    {
        // check space dependency
        MORIS_ERROR( this->check_space_dependency(),
                "Property::dPropdxdDOF - No dependency on space." );

        // check dof dependency
        MORIS_ERROR( this->check_dof_dependency( aDofType ),
                "Property::dPropdxdDOF - No dependency on dof." );

        // get the dof index
        uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative of the property wrt space and dof was not evaluated
        if ( mPropSpaceDofDerEval( aOrder - 1, tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dPropdxdDOF( aDofType, aOrder );

            // set bool for evaluation
            mPropSpaceDofDerEval( aOrder - 1, tDofIndex ) = false;
        }

        // return the value
        return mPropSpaceDofDer( aOrder - 1 )( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    Property::eval_dPropdxdDOF(
            const Vector< MSI::Dof_Type >& aDofType,
            const uint&                    aOrder )
    {
        // get the dof index
        uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // check that mSpaceDofDerFunctions was assigned
        MORIS_ASSERT( mSetSpaceDofDerFunctions == true,
                "Property::eval_dPropdxdDOF - mSetSpaceDofDerFunctions not assigned. " );
        MORIS_ASSERT( mSpaceDofDerFunctions( tDofIndex ) != nullptr,
                "Property::eval_dPropdxdDOF - mSpaceDofDerFunctions not assigned. " );

        // check that mFIManager was assigned
        MORIS_ASSERT( mFIManager != nullptr,
                "Property::eval_dPropdxdDOF - mFIManager not assigned. " );

        // if so use mSpaceDofDerFunctions to compute the derivative
        mSpaceDofDerFunctions( tDofIndex )( mPropSpaceDofDer( aOrder - 1 )( tDofIndex ), mParameters, mFIManager );
    }

    //------------------------------------------------------------------------------

    void
    Property::set_dv_derivative_functions(
            const Vector< PropertyFunc >& aDvDerFunctions )
    {
        // set functions for derivatives wrt dv
        mDvDerFunctions = aDvDerFunctions;

        // set setting flag
        mSetDvDerFunctions = true;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Property::dPropdDV(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // if aDvType is not an active dv type for the property
        MORIS_ERROR( this->check_dv_dependency( aDvType ),
                "Property::dPropdDV - no dependency in this dv type." );

        // get dv type index
        uint tDvIndex = mDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mPropDvDerEval( tDvIndex ) )
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

    void
    Property::eval_dPropdDV(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // get the dv index
        uint tDvIndex = mDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

        // check that mDofDerFunctions was assigned
        MORIS_ASSERT( mSetDvDerFunctions == true,
                "Property::eval_dPropdDV - mDvDerFunctions not assigned. " );
        MORIS_ASSERT( mDvDerFunctions( tDvIndex ) != nullptr,
                "Property::eval_dPropdDV - mDvDerFunctions not assigned. " );

        // check that mFIManager was assigned
        MORIS_ASSERT( mFIManager != nullptr,
                "Property::eval_Prop - mFIManager not assigned. " );

        // if so use mDerivativeFunction to compute the derivative
        mDvDerFunctions( tDvIndex )( mPropDvDer( tDvIndex ), mParameters, mFIManager );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
