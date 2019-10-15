/*
 * cl_FEM_Property_User_Defined_Info.hpp
 *
 *  Created on: Aug 23, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include <functional>

namespace moris
{
    namespace fem
    {
    class Field_Interpolator;
    class Geometry_Interpolator;

    typedef std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                              moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                              moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                              fem::Geometry_Interpolator              * aGeometryInterpolator ) > PropertyFunc;

//------------------------------------------------------------------------------
        /**
         * Property_User_Defined_Info
         */
        class Property_User_Defined_Info
        {
        protected :

        // property type
        fem::Property_Type mPropertyType;

        // property dof type dependency list
        moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

        // property dv type dependency list
        moris::Cell< moris::Cell< MSI::Dv_Type > > mDvTypes;

        // property coeff list
        moris::Cell< Matrix< DDRMat > > mParamList;

        // property value function list
        PropertyFunc mValFunc;

        // property derivative functions list wrt dof
        moris::Cell< PropertyFunc > mDofDerFuncList;

        // property derivative functions list wrt dv
        moris::Cell< PropertyFunc > mDvDerFuncList;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property_User_Defined_Info(){};

            /**
             * constructor without dv dependency
             * @param[ in ] aPropertyType   property type
             * @param[ in ] aDofTypes       list of dof types
             * @param[ in ] aParamList      list of parameters
             * @param[ in ] aValFunc        property function for property value
             * @param[ in ] aDofDerFuncList list of property functions for derivatives wrt dof
             */
            Property_User_Defined_Info( fem::Property_Type                          aPropertyType,
                                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                        moris::Cell< Matrix< DDRMat > >             aParamList,
                                        PropertyFunc                                aValFunc,
                                        moris::Cell< PropertyFunc >                 aDofDerFuncList )
            : mPropertyType( aPropertyType ),
              mDofTypes( aDofTypes ),
              mParamList( aParamList ),
              mValFunc( aValFunc ),
              mDofDerFuncList( aDofDerFuncList )
            {};

            /**
             * constructor with dv dependency
             * @param[ in ] aPropertyType   property type
             * @param[ in ] aDofTypes       list of dof types
             * @param[ in ] aDvTypes        list of dv types
             * @param[ in ] aParamList      list of parameters
             * @param[ in ] aValFunc        property function for property value
             * @param[ in ] aDofDerFuncList list of property functions for derivatives wrt dof
             * @param[ in ] aDvDerFuncList  list of property functions for derivatives wrt dv
             */
            Property_User_Defined_Info( fem::Property_Type                          aPropertyType,
                                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                        moris::Cell< moris::Cell< MSI::Dv_Type > >  aDvTypes,
                                        moris::Cell< Matrix< DDRMat > >             aParamList,
                                        PropertyFunc                                aValFunc,
                                        moris::Cell< PropertyFunc >                 aDofDerFuncList,
                                        moris::Cell< PropertyFunc >                 aDvDerFuncList )
            : mPropertyType( aPropertyType ),
              mDofTypes( aDofTypes ),
              mDvTypes( aDvTypes ),
              mParamList( aParamList ),
              mValFunc( aValFunc ),
              mDofDerFuncList( aDofDerFuncList ),
              mDvDerFuncList( aDvDerFuncList )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
             ~Property_User_Defined_Info(){};

//------------------------------------------------------------------------------
            /**
             * set property type
             * @param[ in ] aPropertyType property type
             */
            void set_property_type( fem::Property_Type aPropertyType )
            {
                mPropertyType = aPropertyType;
            };

//------------------------------------------------------------------------------
            /**
             * return property type
             * @param[ out ] mPropertyType property type
             */
            const fem::Property_Type & get_property_type() const
            {
                return mPropertyType;
            };

//------------------------------------------------------------------------------
            /**
             * set property dof type dependency
             * @param[ in ] aDofTypes list of dof types
             */
            void set_property_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
            {
                mDofTypes = aDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * return property dof type dependency
             * @param[ out ] mDofTypes list of dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_property_dof_type_list() const
            {
                return mDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * set property dv type dependency
             * @param[ in ] aDvTypes list of dv types
             */
            void set_property_dv_type_list( moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes )
            {
                mDvTypes = aDvTypes;
            }

//------------------------------------------------------------------------------
            /**
             * return property dv type dependency
             * @param[ out ] aDvTypes list of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_property_dv_type_list() const
            {
                return mDvTypes;
            }

//------------------------------------------------------------------------------
            /**
             * set property parameter list
             * @param[ in ] aParamList list of parameters
             */
            void set_property_param_list( moris::Cell< Matrix< DDRMat > > & aParamList )
            {
                mParamList = aParamList;
            };

//------------------------------------------------------------------------------
            /**
             * return property parameter list
             * @param[ out ] aParamList list of parameters
             */
            const moris::Cell< Matrix< DDRMat > > & get_property_param_list() const
            {
                return mParamList;
            };

//------------------------------------------------------------------------------
            /**
             * set property value function
             * @param[ in ] aValFunc property function
             */
            void set_property_valFunc( PropertyFunc & aValFunc )
            {
                mValFunc = aValFunc;
            };

//------------------------------------------------------------------------------
            /**
             * return property value function
             * @param[ out ] aValFunc property function
             */
            const PropertyFunc & get_property_valFunc() const
            {
                return mValFunc;
            };

//------------------------------------------------------------------------------
            /**
             * set property derivative function list wrt dof
             * @param[ in ] aDofDerFuncList list of property functions for derivative wrt dof
             */
            void set_property_dof_derFunc_list( moris::Cell< PropertyFunc > & aDofDerFuncList )
            {
                mDofDerFuncList = aDofDerFuncList;
            };

//------------------------------------------------------------------------------
            /**
             * return property derivative function list wrt dof
             * @param[ in ] mDofDerFuncList list of property functions for derivative wrt dof
             */
            const moris::Cell< PropertyFunc > & get_property_dof_derFunc_list() const
            {
                return mDofDerFuncList;
            };

//------------------------------------------------------------------------------
            /**
             * set property derivative function list wrt dv
             * @param[ in ] aDvDerFuncList list of property functions for derivative wrt dv
             */
            void set_property_dv_derFunc_list( moris::Cell< PropertyFunc > & aDvDerFuncList )
            {
                mDvDerFuncList = aDvDerFuncList;
            };

//------------------------------------------------------------------------------
            /**
             * return property derivative function list wrt dv
             * @param[ in ] mDvDerFuncList list of property functions for derivative wrt dv
             */
            const moris::Cell< PropertyFunc > & get_property_dv_derFunc_list() const
            {
                return mDvDerFuncList;
            };

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_ */
