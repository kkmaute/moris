
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
                            moris::Cell< PropertyFunc >                 aDofDerFunctions,
                            Geometry_Interpolator*                      aGeometryInterpolator )
        : mPropertyType( aPropertyType ),
          mDofTypes( aDofTypes ),
          mParameters( aParameters ),
          mValFunction( aValFunction ),
          mDofDerFunctions( aDofDerFunctions ),
          mGeometryInterpolator( aGeometryInterpolator )
        {
            // build a dof type map
            this->build_dof_type_map();

            // set mDofDerFunctions size
            mDofDerFunctions.resize( mDofTypes.size(), nullptr );

            // set mPropDofDerEval size
            mPropDofDerEval.resize( mDofTypes.size(), true );

            // set mPropDofDer size
            mPropDofDer.resize( mDofTypes.size() );
        }

        Property::Property( fem::Property_Type                          aPropertyType,
                            moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                            moris::Cell< moris::Cell< MSI::Dv_Type > >  aDvTypes,
                            moris::Cell< moris::Matrix< DDRMat > >      aParameters,
                            PropertyFunc                                aValFunction,
                            moris::Cell< PropertyFunc >                 aDofDerFunctions,
                            moris::Cell< PropertyFunc >                 aDvDerFunctions,
                            Geometry_Interpolator*                      aGeometryInterpolator )
        : mPropertyType( aPropertyType ),
          mDofTypes( aDofTypes ),
          mDvTypes( aDvTypes ),
          mParameters( aParameters ),
          mValFunction( aValFunction ),
          mDofDerFunctions( aDofDerFunctions ),
          mDvDerFunctions( aDvDerFunctions ),
          mGeometryInterpolator( aGeometryInterpolator )
        {
            // build a dof type map
            this->build_dof_type_map();

            // build a dv type map
            this->build_dv_type_map();

            // number of dof types
            uint tNumDofTypes = mDofTypes.size();

            // number of dv types
            uint tNumDvTypes = mDvTypes.size();

            // set mDofDerFunctions size
            mDofDerFunctions.resize( tNumDofTypes, nullptr );
            mDvDerFunctions.resize( tNumDvTypes, nullptr );

            // set mPropDofDerEval size
            mPropDofDerEval.resize( tNumDofTypes, true );
            mPropDvDerEval.resize( tNumDvTypes, true );

            // set mPropDofDer size
            mPropDofDer.resize( tNumDofTypes );
            mPropDvDer.resize( tNumDvTypes );
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
        void Property::set_dof_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
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
            mDofFI = aFieldInterpolators;
        }

//------------------------------------------------------------------------------
        void Property::check_dof_field_interpolators()
        {
            // check field interpolators cell size
            MORIS_ASSERT( mDofFI.size() == mDofTypes.size(),
                          "Property::check_field_interpolators - wrong FI size. " );

           // loop over the field interpolator pointers
           for( uint iFI = 0; iFI < mDofTypes.size(); iFI++ )
           {
               // check that the field interpolator was set
               MORIS_ASSERT( mDofFI( iFI ) != nullptr,
                             "Property::check_field_interpolators - FI missing. " );
           }
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
        bool Property::check_dv_dependency( const moris::Cell< MSI::Dv_Type > aDvType )
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
        void Property::set_dv_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
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
        void Property::check_dv_field_interpolators()
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
            this->check_dof_field_interpolators();

            // use mValFunction to evaluate the property
            mProp = mValFunction( mParameters, mDofFI, mDvFI, mGeometryInterpolator );

            // set bool for evaluation
            mPropEval = false;
        }

//------------------------------------------------------------------------------
        const Matrix< DDRMat > & Property::dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
        {
           // if aDofType is not an active dof type for the property
           MORIS_ERROR( this->check_dof_dependency( aDofType ), "Property::dPropdDOF - no dependency in this dof type." );

           // get the dof index
           uint tDofIndex = mDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

           // if the derivative has not been evaluated yet
           if( mPropDofDerEval( tDofIndex ) )
           {
               // evaluate the derivative
               this->eval_dPropdDOF( aDofType );
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
            MORIS_ASSERT( mDofDerFunctions( tDofIndex ) != nullptr, "Property::dPdDOF - mDerFunction not assigned. " );

            // check that mGeometryInterpolator was assigned
            MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_dPropdDOF - mGeometryInterpolator not assigned. " );

            // check that the field interpolators are set
            this->check_dof_field_interpolators();

            // if so use mDerivativeFunction to compute the derivative

            mPropDofDer( tDofIndex ) = mDofDerFunctions( tDofIndex )( mParameters,
                                                                   mDofFI,
                                                                   mDvFI,
                                                                   mGeometryInterpolator );
        }

//------------------------------------------------------------------------------
        const Matrix< DDRMat > & Property::dPropdDV( const moris::Cell< MSI::Dv_Type > aDvType )
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
           }

           // return the derivative
           return mPropDvDer( tDvIndex );
        }

//------------------------------------------------------------------------------
        void Property::eval_dPropdDV( const moris::Cell< MSI::Dv_Type > aDvType )
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

            mPropDvDer( tDvIndex ) = mDvDerFunctions( tDvIndex )( mParameters,
                                                                  mDofFI,
                                                                  mDvFI,
                                                                  mGeometryInterpolator );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
