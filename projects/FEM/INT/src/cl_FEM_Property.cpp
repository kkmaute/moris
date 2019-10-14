
#include "cl_FEM_Property.hpp" //FEM/MSI/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        Property::Property( fem::Property_Type                          aPropertyType,
                            moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                            moris::Cell< moris::Matrix< DDRMat > >      aParameters,
                            PropertyFunc                                aValFunction,
                            moris::Cell< PropertyFunc >                 aDerFunctions,
                            Geometry_Interpolator*                      aGeometryInterpolator )
        : mPropertyType( aPropertyType ),
          mDofTypes( aDofTypes ),
          mParameters( aParameters ),
          mValFunction( aValFunction ),
          mDerFunctions( aDerFunctions ),
          mGeometryInterpolator( aGeometryInterpolator )
        {
            // build a dof type map
            this->build_dof_type_map();

            // set mDerFunctions size
            mDerFunctions.resize( mDofTypes.size(), nullptr );

            // set mPropDerEval size
            mPropDerEval.resize( mDofTypes.size(), true );

            // set mPropDer size
            mPropDer.resize( mDofTypes.size() );
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
        bool Property::check_dof_dependency( const moris::Cell< MSI::Dof_Type > aDofType )
        {
            // set bool for dependency
            bool tDofDependency = false;

            // if aDofType is an active dof type for the property
            if( static_cast< uint >( aDofType( 0 ) ) < mDofTypeMap.numel() && mDofTypeMap(static_cast< uint >( aDofType( 0 ) ) ) != -1 )
            {
                // bool is set to true
                tDofDependency = true;
            }
            // return bool for dependency
            return tDofDependency;
        }

//------------------------------------------------------------------------------
        void Property::set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check size
            MORIS_ASSERT( aFieldInterpolators.size() == mDofTypes.size(),
                          "Property::set_field_interpolators - wrong input size. " );

            // check field interpolator type
            bool tCheckFI = true;
            for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
            {
                tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == mDofTypes( iFI )( 0 ) );
            }
            MORIS_ASSERT( tCheckFI, "Property::set_field_interpolators - wrong field interpolator dof type. ");

            // set field interpolators
            mFieldInterpolators = aFieldInterpolators;
        }

//------------------------------------------------------------------------------
        void Property::check_field_interpolators()
        {
            // check field interpolators cell size
            MORIS_ASSERT( mFieldInterpolators.size() == mDofTypes.size(),
                          "Property::check_field_interpolators - wrong FI size. " );

           // loop over the field interpolator pointers
           for( uint iFI = 0; iFI < mDofTypes.size(); iFI++ )
           {
               // check that the field interpolator was set
               MORIS_ASSERT( mFieldInterpolators( iFI ) != nullptr,
                             "Property::check_field_interpolators - FI missing. " );
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
            }
            // return the property value
            return mProp;
        }

//------------------------------------------------------------------------------
        void Property::eval_Prop()
        {
            // check that mValFunc was assigned
            MORIS_ASSERT( mValFunction != nullptr, "Property::eval_val - mValFunction not assigned. " );

            // check that mGeometryInterpolator was assigned
            MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_val - mGeometryInterpolator not assigned. " );

            // check that the field interpolators are set
            this->check_field_interpolators();

            // use mValFunction to evaluate the property
            mProp = mValFunction( mParameters, mFieldInterpolators, mGeometryInterpolator );

            // set bool for evaluation
            mPropEval = false;
        }

//------------------------------------------------------------------------------
        const Matrix< DDRMat > & Property::dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
        {
           // if aDofType is not an active dof type for the property
           MORIS_ERROR( this->check_dof_dependency( aDofType ), "Property::dPropdDOF - no dependency in this dof type." );

           // if the derivative has not been evaluated yet
           if( mPropDerEval( mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) ) ) )
           {
               // evaluate the derivative
               this->eval_dPropdDOF( aDofType );
           }

           // retuirn the derivative
           return mPropDer( mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) ) );
        }

//------------------------------------------------------------------------------
        void Property::eval_dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
        {
            // check that mDerFunctions was assigned
            MORIS_ASSERT( mDerFunctions( mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) ) ) != nullptr, "Property::dPdDOF - mDerFunction not assigned. " );

            // check that mGeometryInterpolator was assigned
            MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_dPropdDOF - mGeometryInterpolator not assigned. " );

            // check that the field interpolators are set
            this->check_field_interpolators();

            // if so use mDerivativeFunction to compute the derivative
            mPropDer( mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) ) )
                = mDerFunctions( mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) ) )( mParameters, mFieldInterpolators, mGeometryInterpolator );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
