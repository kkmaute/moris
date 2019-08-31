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
    typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                              moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator) > PropertyFunc;
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
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // active dof type map
            Matrix< DDSMat > mDofTypeMap;

            // field interpolators
            moris::Cell< Field_Interpolator* > mFieldInterpolators;

            // coefficients
            moris::Cell< Matrix< DDRMat> > mCoeff;

            // value function
            PropertyFunc mValFunction = nullptr;

            // derivative functions
            moris::Cell< PropertyFunc > mDerFunctions;

            // flag for evaluation
            bool mPropEval = true;
            moris::Cell< bool > mPropDerEval;

            // storage
            Matrix< DDRMat > mProp;
            moris::Cell< Matrix< DDRMat > > mPropDer;
            Matrix< DDRMat > mPropDerZero;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property( fem::Property_Type aPropertyType,
                      moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes ) : mPropertyType( aPropertyType ),
                                                                                mDofTypes( aDofTypes )
            {};

            Property( fem::Property_Type                          aPropertyType,
                      moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                      PropertyFunc                                aValFunction,
                      moris::Cell< PropertyFunc >                 aDerFunctions );

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
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list() const
            {
                return mDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * build a dof type map
             */
            void build_dof_type_map();

//------------------------------------------------------------------------------
            /**
             * set field interpolators
             */
            void set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators );

            moris::Cell< fem::Field_Interpolator* > & get_field_interpolators()
            {
                return mFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * check that field interpolators are assigned
             */
            void check_field_interpolators();

//------------------------------------------------------------------------------
            /**
             * set coefficients
             */
            void set_coefficients( const moris::Cell< Matrix< DDRMat > > & aCoefficients );

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags()
            {
                mPropEval = true;

                mPropDerEval.resize( mDofTypes.size(), true );
            }

//------------------------------------------------------------------------------
            /**
             * get the property value
             */
            const Matrix< DDRMat > & val();

//------------------------------------------------------------------------------
            /**
             * evaluate property in terms of coefficients and variables
             */
            void eval_Prop();

//------------------------------------------------------------------------------
            /**
             * evaluate property derivatives in terms of coefficients and variables
             */
            const Matrix< DDRMat > & dPropdDOF( MSI::Dof_Type aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate property derivatives in terms of coefficients and variables
             */
            void eval_dPropdDOF( MSI::Dof_Type aDofType );

        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_HPP_ */
