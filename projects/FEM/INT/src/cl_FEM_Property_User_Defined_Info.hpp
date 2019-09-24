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
                                              moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
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

        // property coeff list
        moris::Cell< Matrix< DDRMat > > mParamList;

        // property value function list
        PropertyFunc mValFunc;

        // property derivative functions list
        moris::Cell< PropertyFunc > mDerFuncList;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Property_User_Defined_Info(){};

            Property_User_Defined_Info( fem::Property_Type                          aPropertyType,
                                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                        moris::Cell< Matrix< DDRMat > >             aParamList,
                                        PropertyFunc                                aValFunc,
                                        moris::Cell< PropertyFunc >                 aDerFuncList )
            : mPropertyType( aPropertyType ),
              mDofTypes( aDofTypes ),
              mParamList( aParamList ),
              mValFunc( aValFunc ),
              mDerFuncList( aDerFuncList )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
             ~Property_User_Defined_Info(){};

////------------------------------------------------------------------------------
            /**
             * returns property type
             */
            const fem::Property_Type & get_property_type() const
            {
                return mPropertyType;
            };

//------------------------------------------------------------------------------
            /**
             * returns property dof type dependency
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_property_dof_type_list() const
            {
                return mDofTypes;
            };
//------------------------------------------------------------------------------
            /**
             * returns property coeff list
             */
            const moris::Cell< Matrix< DDRMat > > & get_property_param_list() const
            {
                return mParamList;
            };
//------------------------------------------------------------------------------
            /**
             * returns property value function
             */
            const PropertyFunc & get_property_valFunc() const
            {
                return mValFunc;
            };
//------------------------------------------------------------------------------
            /**
             * returns property derivative function list
             */
            const moris::Cell< PropertyFunc > & get_property_derFunc_list() const
            {
                return mDerFuncList;
            };

        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_USER_DEFINED_INFO_HPP_ */
