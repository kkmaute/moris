
#include "cl_FEM_Property.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        /**
         * constructor
         */
        Property::Property( fem::Property_Type                          aPropertyType,
                            moris::Cell< moris::Cell< MSI::Dof_Type > > aActiveDofTypes,
                            PropertyFunc                                aValFunction,
                            moris::Cell< PropertyFunc >                 aDerFunctions )
        : mPropertyType( aPropertyType ),
          mDofTypes( aActiveDofTypes ),
          mValFunction( aValFunction ),
          mDerFunctions( aDerFunctions )
        {
            // build a dof type map
            this->build_dof_type_map();

            // set mDerFunctions size
            mDerFunctions.resize( mDofTypes.size(), nullptr );

            // set mPropDerEval size
            mPropDerEval.resize( mDofTypes.size(), true );

            // set mPropDer size
            mPropDer.resize( mDofTypes.size() );

            // init mPropDerZero size
            mPropDerZero.set_size( 1, 1, 0.0 );
        }


        Property::Property( fem::Property_Type                          aPropertyType,
                            moris::Cell< moris::Cell< MSI::Dof_Type > > aActiveDofTypes,
                            PropertyFunc                                aValFunction,
                            moris::Cell< PropertyFunc >                 aDerFunctions,
                            Geometry_Interpolator*                      aGeometryInterpolator )
        : mPropertyType( aPropertyType ),
          mDofTypes( aActiveDofTypes ),
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

            // init mPropDerZero size
            mPropDerZero.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------
        /**
         * build a dof type map
         */
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
        /**
         * set field interpolators
         */
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

        /**
         * check that field interpolators were assigned
         */
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
        /**
         * set coefficients
         */
        void Property::set_coefficients( const moris::Cell< Matrix< DDRMat > > & aCoefficients )
        {
             // fixme check right number of coeff?

            // set coefficients
            mCoeff = aCoefficients;
        }

//------------------------------------------------------------------------------
        /**
         * get the property value
         */
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
        /**
         * evaluate property in terms of coefficients and variables
         */
        void Property::eval_Prop()
        {
            // check that mValFunc was assigned
            MORIS_ASSERT( mValFunction != nullptr, "Property::eval_val - mValFunction not assigned. " );

            // check that mGeometryInterpolator was assigned
            MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_val - mGeometryInterpolator not assigned. " );

            // check that the field interpolators are set
            this->check_field_interpolators();

            // use mValFunction to evaluate the property
            mProp = mValFunction( mCoeff, mFieldInterpolators, mGeometryInterpolator );

            // set bool for evaluation
            mPropEval = false;
        }

//------------------------------------------------------------------------------
        /**
         * evaluate property derivatives in terms of coefficients and variables
         */
        const Matrix< DDRMat > & Property::dPropdDOF( MSI::Dof_Type aDofType )
        {
           // if aDofType is an active dof type for the property
           if( static_cast< uint >( aDofType ) < mDofTypeMap.numel() && mDofTypeMap(static_cast< uint >( aDofType ) ) != -1 )
           {
               if( mPropDerEval( mDofTypeMap(static_cast< uint >( aDofType ) ) ) )
               {
                   this->eval_dPropdDOF( aDofType );
               }
               return mPropDer( mDofTypeMap( static_cast< uint >( aDofType ) ) );
           }
           // if adofType is not an active dof type for the property
           else
           {
               return mPropDerZero;
           }
        }

//------------------------------------------------------------------------------
        /**
         * evaluate property derivatives in terms of coefficients and variables
         */
        void Property::eval_dPropdDOF( MSI::Dof_Type aDofType )
        {
            // check that mDerFunctions was assigned
            MORIS_ASSERT( mDerFunctions( mDofTypeMap( static_cast< uint >( aDofType ) ) ) != nullptr, "Property::dPdDOF - mDerFunction not assigned. " );

            // check that mGeometryInterpolator was assigned
            MORIS_ASSERT( mGeometryInterpolator != nullptr, "Property::eval_dPropdDOF - mGeometryInterpolator not assigned. " );

            // check that the field interpolators are set
            this->check_field_interpolators();

            // if so use mDerivativeFunction to compute the derivative
            mPropDer( mDofTypeMap( static_cast< uint >( aDofType ) ) )
                = mDerFunctions( mDofTypeMap( static_cast< uint >( aDofType ) ) )( mCoeff, mFieldInterpolators, mGeometryInterpolator );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
