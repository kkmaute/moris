/*
 * cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp
 *
 *  Created on: Jul 20, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINMATIC_VISCOSITY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINMATIC_VISCOSITY_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Turbulent_Kinematic_Viscosity : public IQI
        {
            private:

                // viscosity dof type (default)
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;

                // property type for IQI
                enum class Property_Type
                {
                    DYNAMIC_VISCOSITY, // fluid dynamic viscosity
                    DENSITY,           // fluid density
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, Property_Type > mPropertyMap;

                // Spalart-Allmaras turbulence model constants
                real mCv1 = 7.1;

                //------------------------------------------------------------------------------

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Turbulent_Kinematic_Viscosity();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Turbulent_Kinematic_Viscosity(){};

                //------------------------------------------------------------------------------
                /**
                 * set dof type list for IQI
                 * @param[ in ] aDofTypes   list of group of dof types
                 * @param[ in ] aDofStrings list of names for group of dof types
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aQI quantity of interest matrix to fill
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_QI( moris::real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest
                 * wrt requested dof types
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dQIdu( moris::real aWStar );

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------
                /**
                 * compute fv1 = chi³ / ( chi³ + cv1³)
                 * @param[ out ] fv1
                 */
                real compute_fv1();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv1 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
                 */
                void compute_dfv1du(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adfv1du );

                //------------------------------------------------------------------------------
                /**
                 * compute chi = viscosityDof / viscosityProp
                 * @param[ out ] chi
                 */
                real compute_chi();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidu    a matrix to fill with dchidu
                 */
                void compute_dchidu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adchidu );

        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINEMATIC_VISCOSITY_HPP_ */
