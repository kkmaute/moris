/*
 * cl_FEM_Property.hpp
 *
 *  Created on: Jul 12, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PROPERTY_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include "fn_equal_to.hpp"
#include <functional>

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Property
         */
        class Property
        {
        protected :

            // property type
            fem::Property_Type mPropertyType;

            // active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > mActiveDofTypes;

            // field interpolators
            moris::Cell< Field_Interpolator* > mFieldInterpolators;

            // coefficients
            moris::Cell< Matrix< DDRMat> > mCoeff;

            // val function
//            Matrix< DDRMat > ( *mValFunction )( moris::Cell< Matrix< DDRMat > >    & aCoeff,
//                                                moris::Cell< Field_Interpolator* > & aFieldInterpolator ) = nullptr;
            std::function< Matrix< DDRMat >( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                             moris::Cell< Field_Interpolator* > & aFieldInterpolator ) > mValFunction = nullptr;
            // derivative functions
            moris::Cell< std::function< Matrix< DDRMat >( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                                          moris::Cell< Field_Interpolator* > & aFieldInterpolator ) > > mDerFunctions;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property( fem::Property_Type aPropertyType,
                      moris::Cell< moris::Cell< MSI::Dof_Type > > aActiveDofTypes ) : mPropertyType( aPropertyType ),
                                                                                      mActiveDofTypes( aActiveDofTypes )
            {};

            Property( fem::Property_Type aPropertyType,
                      moris::Cell< moris::Cell< MSI::Dof_Type > > aActiveDofTypes,
                      std::function< Matrix< DDRMat >( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                                       moris::Cell< Field_Interpolator* > & aFieldInterpolator ) > aValFunction,
                      moris::Cell< std::function< Matrix< DDRMat >( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                                                    moris::Cell< Field_Interpolator* > & aFieldInterpolator ) > > aDerFunctions )
            : mPropertyType( aPropertyType ),
              mActiveDofTypes( aActiveDofTypes ),
              mValFunction( aValFunction ),
              mDerFunctions( aDerFunctions )
            {
                // set mDerFunctions size
                mDerFunctions.resize( mActiveDofTypes.size(), nullptr );
            };

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Property(){};

//------------------------------------------------------------------------------
            /**
             * returns property type
             */
            fem::Property_Type get_property_type() const
            {
                return mPropertyType;
            };

//------------------------------------------------------------------------------
            /**
             * returns a list of active dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_active_dof_types() const
            {
                return mActiveDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * set field interpolators
             */
            void set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                // check size
                MORIS_ASSERT( aFieldInterpolators.size() == mActiveDofTypes.size(), "Property::set_field_interpolators - wrong input size. " );

                // set field interpolators
                mFieldInterpolators = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * set coefficients
             */
            void set_coefficients( moris::Cell< Matrix< DDRMat > > & aCoefficients )
            {
                // set coefficients
                mCoeff = aCoefficients;
            }
//------------------------------------------------------------------------------
            /**
             * evaluate property in terms of coefficients and variables
             */
            void val( Matrix< DDRMat > & aPropertyValue )
            {
                // check that mValFunc was assigned
                MORIS_ASSERT( mValFunction != nullptr, "Property::val - mValFunction not assigned. " );

                // check that the field interpolators are set
                MORIS_ASSERT( mFieldInterpolators.size() > 0, "Property::val - mFieldInterpolators not assigned. " );

                // use mValFunction to evaluate the property
                aPropertyValue = mValFunction( mCoeff, mFieldInterpolators );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate property derivatives in terms of coefficients and variables
             */
            void dPdDOF( MSI::Dof_Type                        aDofType,
                         Matrix< DDRMat >                   & aPropertyDerivative )
            {
                // check that the field interpolators are set
                MORIS_ASSERT( mFieldInterpolators.size() > 0, "Property::val - mFieldInterpolators not assigned. " );

                // loop over the property dof dependencies
                bool tDerEval = false;
                for( uint iDOF = 0; iDOF < mActiveDofTypes.size(); iDOF++ )
                {
                    // check aDofType = active DofType
                    if ( equal_to( static_cast< uint >( aDofType ), static_cast< uint >( mActiveDofTypes( iDOF )( 0 ) ) ) )
                    {
                        // check that mDerFunctions was assigned
                        MORIS_ASSERT( mDerFunctions( iDOF ) != nullptr, "Property::dPdDOF - mDerFunction not assigned. " );

                        // if so use mDerivativeFunction to compute the derivative
                        aPropertyDerivative = mDerFunctions( iDOF )( mCoeff, mFieldInterpolators );

                        // set bool to true as derivative was evaluated
                        tDerEval = true;
                        break;
                    }
                }

                // if property does not depend on aDofType
                if( !tDerEval )
                {
                    // set the derivative to zero
                    aPropertyDerivative.set_size( 1, 1, 0.0 );
                }
            }

        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_HPP_ */
