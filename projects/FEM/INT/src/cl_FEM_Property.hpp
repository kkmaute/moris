/*
 * cl_FEM_Property.hpp
 *
 *  Created on: Jul 12, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PROPERTY_HPP_
#define SRC_FEM_CL_FEM_PROPERTY_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"        //FEM/MSI/src


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

            // active property types
            moris::Cell< fem::Property_Type > mActivePropertyTypes;

            // boolean
            // true  if associated field interpolator uses the coefficients
            // false if associated field interpolator does not use the coefficients
            bool mUseCoeff;

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
            moris::Cell< moris::Cell< MSI::Dof_Type > > get_active_dof_types() const
            {
                return mActiveDofTypes;
            };
//------------------------------------------------------------------------------
            /**
             * returns a list of active property types
             */
            moris::Cell< fem::Property_Type > get_active_property_types() const
            {
                return mActivePropertyTypes;
            };

//------------------------------------------------------------------------------
            /**
             * evaluate coefficients
             */
            virtual void val_coeff( Matrix< DDRMat > & aCoeff )
            {
                MORIS_ERROR( false, "Property::val_coeff - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate property in terms of x and t
             */
            virtual void val( Matrix< DDRMat > & aCoeff,
                              Matrix< DDRMat > & aSpacePhysCoords,
                              Matrix< DDRMat > & aTimePhysCoords)
            {
                MORIS_ERROR( false, "Property::val - This function does nothing. " );
            }

        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PROPERTY_HPP_ */
