/*
 * cl_FEM_Constitutive_User_Defined_Info.hpp
 *
 *  Created on: Sep 18, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_USER_DEFINED_INFO_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_USER_DEFINED_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
    class Field_Interpolator;
    class Property;

//------------------------------------------------------------------------------
        /**
         * Constitutive_User_Defined_Info
         */
        class Constitutive_User_Defined_Info
        {
        protected :

        // constitutive type
        fem::Constitutive_Type mConstitutiveType;

        // constitutive model dof type list
        moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

        // constitutive model property type list
        moris::Cell< fem::Property_Type > mPropTypes;

        // spatial dimension
        uint mSpaceDim;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Constitutive_User_Defined_Info(){};

            Constitutive_User_Defined_Info( fem::Constitutive_Type                      aConstitutiveType,
                                            moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                            moris::Cell< fem::Property_Type >           aPropTypes,
                                            uint                                        aSpaceDim )
            : mConstitutiveType( aConstitutiveType ),
              mDofTypes( aDofTypes ),
              mPropTypes( aPropTypes ),
              mSpaceDim( aSpaceDim )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
             ~Constitutive_User_Defined_Info(){};

//------------------------------------------------------------------------------
            /**
             * sets constitutive type
             */
            void set_constitutive_type( fem::Constitutive_Type aConstitutiveType )
            {
                mConstitutiveType = aConstitutiveType;
            };

//------------------------------------------------------------------------------
            /**
             * returns constitutive type
             */
            const fem::Constitutive_Type & get_constitutive_type() const
            {
                return mConstitutiveType;
            };

//------------------------------------------------------------------------------
            /**
             * sets constitutive model dof types
             */
            void set_constitutive_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
            {
                mDofTypes = aDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * returns constitutive model dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_constitutive_dof_type_list() const
            {
                return mDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * sets constitutive model property types
             */
            void set_constitutive_property_type_list( moris::Cell< fem::Property_Type > & aPropTypes )
            {
                mPropTypes = aPropTypes;
            };

//------------------------------------------------------------------------------
            /**
             * returns constitutive model property types
             */
            const moris::Cell< fem::Property_Type > & get_constitutive_property_type_list() const
            {
                return mPropTypes;
            };

//------------------------------------------------------------------------------
            /**
             * sets spatial dimension
             */
            void set_constitutive_space_dim( uint aSpaceDim )
            {
                mSpaceDim = aSpaceDim;
            };

//------------------------------------------------------------------------------
            /**
             * returns spatial dimension
             */
            uint get_constitutive_space_dim() const
            {
                return mSpaceDim;
            };

        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_USER_DEFINED_INFO_HPP_ */
