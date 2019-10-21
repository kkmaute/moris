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

        // constitutive model dv type list
        moris::Cell< moris::Cell< MSI::Dv_Type > > mDvTypes;

        // constitutive model property type list
        moris::Cell< fem::Property_Type > mPropTypes;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Constitutive_User_Defined_Info(){};

            /**
             * constructor
             * @param[ in ] aConstitutiveType constitutive type
             * @param[ in ] aDofTypes         list of group of dof types
             * @param[ in ] aPropTypes        list of property types
             */
            Constitutive_User_Defined_Info( fem::Constitutive_Type                      aConstitutiveType,
                                            moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                            moris::Cell< fem::Property_Type >           aPropTypes )
            : mConstitutiveType( aConstitutiveType ),
              mDofTypes( aDofTypes ),
              mPropTypes( aPropTypes )
            {};

            /**
             * constructor
             * @param[ in ] aConstitutiveType constitutive type
             * @param[ in ] aDofTypes         list of group of dof types
             * @param[ in ] aDvTypes          list of group of dv types
             * @param[ in ] aPropTypes        list of property types
             */
            Constitutive_User_Defined_Info( fem::Constitutive_Type                      aConstitutiveType,
                                            moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                            moris::Cell< moris::Cell< MSI::Dv_Type > >  aDvTypes,
                                            moris::Cell< fem::Property_Type >           aPropTypes )
            : mConstitutiveType( aConstitutiveType ),
              mDofTypes( aDofTypes ),
              mDvTypes( aDvTypes ),
              mPropTypes( aPropTypes )
            {};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
             ~Constitutive_User_Defined_Info(){};

//------------------------------------------------------------------------------
            /**
             * set constitutive type
             * @param[ in ] aConstitutiveType constitutive type
             */
            void set_constitutive_type( fem::Constitutive_Type aConstitutiveType )
            {
                mConstitutiveType = aConstitutiveType;
            };

//------------------------------------------------------------------------------
            /**
             * return constitutive type
             * @param[ out ] mConstitutiveType constitutive type
             */
            const fem::Constitutive_Type & get_constitutive_type() const
            {
                return mConstitutiveType;
            };

//------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes list of group of dof types
             */
            void set_constitutive_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
            {
                mDofTypes = aDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * return constitutive model dof types
             * @param[ out ] aDofTypes list of group of dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_constitutive_dof_type_list() const
            {
                return mDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * set constitutive model dv types
             * @param[ in ] aDvTypes list of group of dv types
             */
            void set_constitutive_dv_type_list( moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes )
            {
                mDvTypes = aDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * return constitutive model dv types
             * @param[ out ] aDvTypes list of group of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_constitutive_dv_type_list() const
            {
                return mDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * set constitutive model property types
             * @param[ in ] aPropTypes list of property types
             */
            void set_constitutive_property_type_list( moris::Cell< fem::Property_Type > & aPropTypes )
            {
                mPropTypes = aPropTypes;
            };

//------------------------------------------------------------------------------
            /**
             * return constitutive model property types
             * @param[ out ] mPropTypes list of property types
             */
            const moris::Cell< fem::Property_Type > & get_constitutive_property_type_list() const
            {
                return mPropTypes;
            };

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_USER_DEFINED_INFO_HPP_ */
