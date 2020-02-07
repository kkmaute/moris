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

#include "cl_GEN_Dv_Enums.hpp"

#include <functional>

namespace moris
{
    namespace fem
    {
    class Set;
    class Field_Interpolator_Manager;

    typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >        & aCoeff,
                                              moris::fem::Field_Interpolator_Manager * aFIManager ) > PropertyFunc;

//------------------------------------------------------------------------------
        /**
         * Property
         */
        class Property
        {
        protected :

            // field interpolator manager pointer
            Field_Interpolator_Manager * mFIManager = nullptr;

            // FEM set pointer
            Set * mSet = nullptr;

            // active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // active dof type map
            Matrix< DDSMat > mDofTypeMap;

            // active dv types
            moris::Cell< moris::Cell< GEN_DV > > mDvTypes;

            // active dv type map
            Matrix< DDSMat > mDvTypeMap;

            // parameters
            moris::Cell< Matrix< DDRMat> > mParameters;

            // value function
            PropertyFunc mValFunction = nullptr;

            // dof and dv derivative functions
            moris::Cell< PropertyFunc > mDofDerFunctions;
            moris::Cell< PropertyFunc > mDvDerFunctions;

            // storage
            Matrix< DDRMat > mProp;
            moris::Cell< Matrix< DDRMat > > mPropDofDer;
            moris::Cell< Matrix< DDRMat > > mPropDvDer;

        private:

            // flag for evaluation
            bool mPropEval = true;
            moris::Cell< bool > mPropDofDerEval;
            moris::Cell< bool > mPropDvDerEval;

            // flag for setting of mValFunction, mDofDerFunctions, mDvDerFunctions
            bool mSetValFunction     = false;
            bool mSetDofDerFunctions = false;
            bool mSetDvDerFunctions  = false;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Property(){};

//------------------------------------------------------------------------------
            /**
             * set field interpolator manager
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager );

//------------------------------------------------------------------------------
            /*
             * set member set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_set_pointer( Set * aSetPointer )
            {
                mSet = aSetPointer;
            }

//------------------------------------------------------------------------------
            /**
             * set parameters
             * @param[ in ] aParameters list of parameters for property evaluation
             */
            void set_parameters( moris::Cell< moris::Matrix< DDRMat > > aParameters )
            {
                mParameters = aParameters;
            };

//------------------------------------------------------------------------------
            /**
             * get parameters
             * @param[ out ] mParameters list of parameters for property evaluation
             */
            const moris::Cell< moris::Matrix< DDRMat > > & get_parameters() const
            {
                return mParameters;
            };

//------------------------------------------------------------------------------
            /**
             * set val function
             * @param[ in ] aValFunction function for property evaluation
             */
            void set_val_function( PropertyFunc aValFunction )
            {
                // set the value function
                mValFunction = aValFunction;

                // set setting flag
                mSetValFunction = true;
            };

//------------------------------------------------------------------------------
            /**
             * get val function
             * @param[ out ] mValFunction function for property evaluation
             */
            const PropertyFunc & get_val_function() const
            {
                return mValFunction;
            };

//------------------------------------------------------------------------------
            /**
             * set dof derivative functions
             * @param[ in ] aDofDerFunctions list function for property derivatives wrt dof
             */
            void set_dof_derivative_functions( moris::Cell< PropertyFunc > aDofDerFunctions )
            {
                // set functions for derivatives wrt dof
                mDofDerFunctions = aDofDerFunctions;

                // set setting flag
                mSetDofDerFunctions = true;
            };

//------------------------------------------------------------------------------
            /**
             * get dof derivative functions
             * @param[ out ] mDofDerFunctions list function for property derivatives wrt dof
             */
            const moris::Cell< PropertyFunc > & get_dof_derivative_functions() const
            {
                return mDofDerFunctions;
            };

//------------------------------------------------------------------------------
            /**
             * set dv derivative functions
             * @param[ in ] aDofDerFunctions list function for property derivatives wrt dv
             */
            void set_dv_derivative_functions( moris::Cell< PropertyFunc > aDvDerFunctions )
            {
                // set functions for derivatives wrt dv
                mDvDerFunctions = aDvDerFunctions;

                // set setting flag
                mSetDvDerFunctions = true;
            };

//------------------------------------------------------------------------------
            /**
             * get dv derivative functions
             * @param[ out ] mDofDerFunctions list function for property derivatives wrt dv
             */
            const moris::Cell< PropertyFunc > & get_dv_derivative_functions() const
            {
                return mDvDerFunctions;
            };

//------------------------------------------------------------------------------
            /**
             * set parameters, val, dof derivative and dv derivative functions
             * @param[ in ] aParameters      list of parameters for property evaluation
             * @param[ in ] aValFunction     function for property value
             * @param[ in ] aDvDerFunctions  list function for property derivatives wrt dof
             * @param[ in ] aDofDerFunctions list function for property derivatives wrt dv
             */
            void set_parameters_and_functions( moris::Cell< moris::Matrix< DDRMat > > aParameters,
                                               PropertyFunc aValFunction,
                                               moris::Cell< PropertyFunc > aDofDerFunctions,
                                               moris::Cell< PropertyFunc > aDvDerFunctions )
            {
                mParameters      = aParameters;
                mValFunction     = aValFunction;
                mDofDerFunctions = aDofDerFunctions;
                mDvDerFunctions  = aDvDerFunctions;

                // set setting flag
                mSetValFunction     = true;
                mSetDofDerFunctions = true;
                mSetDvDerFunctions  = true;
            };

//------------------------------------------------------------------------------
            /**
             * set a list of active dof types
             * @param[ in ] aDofTypes list of dof types
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes )
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
            };

//------------------------------------------------------------------------------
            /**
             * return a list of active dof types
             * @param[ out ] mDofTypes list of dof types
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
             * get a dof type map
             * @param[ out ] mDofTypeMap map of the dof types
             */
            const moris::Matrix< DDSMat > & get_dof_type_map() const
            {
                return mDofTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * check if the property depends on a particular group of dof type
             * @param[ in ]  aDofType cell of dof type
             * @param[ out ] aBool    boolean, true if dependency on the dof type
             */
            bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > aDofType );

//------------------------------------------------------------------------------
            /**
             * set a list of dv types
             * @param[ in ] aDvTypes list of dv type
             */
            void set_dv_type_list( const moris::Cell< moris::Cell< GEN_DV > > & aDvTypes )
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
            };

//------------------------------------------------------------------------------
            /**
             * return a list of dv types
             * @param[ out ] mDvTypes list of dv type
             */
            const moris::Cell< moris::Cell< GEN_DV > > & get_dv_type_list() const
            {
                return mDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * build a dv type map
             */
            void build_dv_type_map();

//------------------------------------------------------------------------------
            /**
             * get a dv type map
             * @param[ out ] mDvTypeMap map of the dv types
             */
            const moris::Matrix< DDSMat > & get_dv_type_map() const
            {
                return mDvTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * check if the property depends on a particular group of dv types
             * @param[ in ]  aDvType cell of dv type
             * @param[ out ] aBool   boolean, true if dependency on the dv type
             */
            bool check_dv_dependency( const moris::Cell< GEN_DV > aDvType );

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags()
            {
                // reset the property value
                mPropEval = true;

                // reset the property derivatives wrt dof type
                mPropDofDerEval.assign( mDofTypes.size(), true );

                // reset the property derivatives wrt dv type
                mPropDvDerEval.assign( mDvTypes.size(), true );
            }

//------------------------------------------------------------------------------
            /**
             * get the property value
             * @param[ out ] aVal matrix with property value
             */
            const Matrix< DDRMat > & val();

//------------------------------------------------------------------------------
            /**
             * evaluate property in terms of coefficients and variables
             */
            void eval_Prop();

//------------------------------------------------------------------------------
            /**
             * get property derivatives wrt a dof type
             * @param[ in ]  aDofType   cell of dof type
             * @param[ out ] adPropdDOF matrix with derivative wrt to the dof type
             */
            const Matrix< DDRMat > & dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate property derivatives wrt a dof type
             * @param[ in ] aDofType cell of dof type
             */
            void eval_dPropdDOF( const moris::Cell< MSI::Dof_Type > aDofType );

//------------------------------------------------------------------------------
            /**
             * get the property derivatives wrt a design variable
             * @param[ in ]  aDvType   cell of dv type
             * @param[ out ] adPropdDV matrix with derivative wrt to the dv type
             */
            const Matrix< DDRMat > & dPropdDV( const moris::Cell< GEN_DV > aDvType );

//------------------------------------------------------------------------------
            /**
             * evaluate property derivatives wrt a design variable
             * @param[ in ] aDvType cell of dv type
             */
            void eval_dPropdDV( const moris::Cell< GEN_DV > aDvType );

//------------------------------------------------------------------------------
            /**
             * get non unique dof type list
             * @param[ in ] aDofType cell of dof type
             */
            void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
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
            /**
             * get non unique dof type list
             * @param[ in ] aDofType cell of dof type
             */
            void get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                  moris::Cell< GEN_DV >  & aDvTypes )
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
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_HPP_ */
