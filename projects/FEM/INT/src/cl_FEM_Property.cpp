
#include "cl_FEM_Property.hpp" //FEM/MSI/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        void Property::reset_eval_flags()
        {
            // reset the property value
            mPropEval = true;

            // reset the property derivative wrt x
            mdPropdxEval = true;

            // reset the property derivatives wrt dof type
            mPropDofDerEval.assign( mDofTypes.size(), true );

            // reset the property derivatives wrt dv type
            mPropDvDerEval.assign( mDvTypes.size(), true );
        }

        //------------------------------------------------------------------------------

        void Property::set_field_interpolator_manager(
                Field_Interpolator_Manager * aFieldInterpolatorManager )
        {
            // set field interpolator manager
            mFIManager = aFieldInterpolatorManager;
        }

        //------------------------------------------------------------------------------

        void Property::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes )
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
            mPropDofDerEval.assign( tNumDofTypes, true );

            // set mPropDofDer size
            mPropDofDer.resize( tNumDofTypes );
        }

        //------------------------------------------------------------------------------

        void Property::build_dof_type_map()
        {
            // determine the max Dof_Type enum
            sint tMaxEnum = 0;
            for( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDof )( 0 ) ) );
            }
            tMaxEnum++;

            // set the Dof_Type map size
            mDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // fill the Dof_Type map
            for( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // fill the property map
                mDofTypeMap( static_cast< int >( mDofTypes( iDof )( 0 ) ), 0 ) = iDof;
            }
        }

        //------------------------------------------------------------------------------

        bool Property::check_dof_dependency(
                const moris::Cell< MSI::Dof_Type > aDofType )
        {
            // set bool for dependency
            bool tDofDependency = false;

            // get the dof type index
            uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

            // if aDofType is an active dof type for the property
            if( tDofIndex < mDofTypeMap.numel() && mDofTypeMap( tDofIndex ) != -1 )
            {
                // bool is set to true
                tDofDependency = true;
            }
            // return bool for dependency
            return tDofDependency;
        }

        //------------------------------------------------------------------------------

        void Property::set_dv_type_list(
                const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes )
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
            mPropDvDerEval.assign( tNumDvTypes, true );

            // set mPropdvDer size
            mPropDvDer.resize( tNumDvTypes );
        }

        //------------------------------------------------------------------------------

        void Property::build_dv_type_map()
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

            // fill the Dof_Type map
            for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
            {
                // fill the property map
                mDvTypeMap( static_cast< int >( mDvTypes( iDV )( 0 ) ), 0 ) = iDV;
            }
        }

        //------------------------------------------------------------------------------

        bool Property::check_dv_dependency( const moris::Cell< PDV_Type > aDvType )
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

        void Property::get_non_unique_dof_types(
                moris::Cell< MSI::Dof_Type > & aDofTypes )
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

        void Property::get_non_unique_dof_and_dv_types(
                moris::Cell< MSI::Dof_Type > & aDofTypes,
                moris::Cell< PDV_Type >  & aDvTypes )
        {
            // init counter
            uint tDofCounter = 0;
            uint tDvCounter  = 0;

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

            // reserve memory for dof and dv type lists
            aDofTypes.reserve( tDofCounter );
            aDvTypes.reserve( tDvCounter );

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
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Property::val()
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

        void Property::eval_Prop()
        {
            // check that mValFunc was assigned
            MORIS_ASSERT( mSetValFunction == true,
                    "Property::eval_Prop - mValFunction not assigned. " );

            // check that mFIManager was assigned
            MORIS_ASSERT( mFIManager != nullptr,
                    "Property::eval_Prop - mFIManager not assigned. " );

            // use mValFunction to evaluate the property
            mValFunction( mProp, mParameters, mFIManager );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Property::dPropdx()
        {
            // if the property was not evaluated
            if( mdPropdxEval )
            {
                // evaluate the property
                this->eval_dPropdx();

                // set bool for evaluation
                mdPropdxEval = false;
            }
            // return the property value
            return mdPropdx;
        }

        //------------------------------------------------------------------------------

        void Property::eval_dPropdx()
        {
            // check that mSetSpaceDerFunc was assigned
            MORIS_ASSERT( mSetSpaceDerFunction == true,
                    "Property::eval_dPropdx - mSpaceDerFunction not assigned. " );

            // check that mFIManager was assigned
            MORIS_ASSERT( mFIManager != nullptr,
                    "Property::eval_dPropdx - mFIManager not assigned. " );

            // use mSpaceDerFunction to evaluate the property
            mSpaceDerFunction( mdPropdx, mParameters, mFIManager );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Property::dPropdDOF(
                const moris::Cell< MSI::Dof_Type > aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR( this->check_dof_dependency( aDofType ),
                    "Property::dPropdDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mPropDofDerEval( tDofIndex ) )
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

        void Property::eval_dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
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

        const Matrix< DDRMat > & Property::dPropdDV(
                const moris::Cell< PDV_Type > aDvType )
        {
            // if aDvType is not an active dv type for the property
            MORIS_ERROR( this->check_dv_dependency( aDvType ),
                    "Property::dPropdDV - no dependency in this dv type." );

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

        void Property::eval_dPropdDV( const moris::Cell< PDV_Type > aDvType )
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
    } /* namespace fem */
} /* namespace moris */
