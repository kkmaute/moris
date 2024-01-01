/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINMATIC_VISCOSITY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINMATIC_VISCOSITY_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_Turbulent_Kinematic_Viscosity : public IQI
        {
            private:

                // viscosity dof type (default)
                MSI::Dof_Type mLeaderDofViscosity = MSI::Dof_Type::VISCOSITY;

                // property type for IQI
                enum class Property_Type
                {
                        KINEMATIC_VISCOSITY, // fluid dynamic viscosity
                        DENSITY,           // fluid density
                        MAX_ENUM
                };

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
                        Vector< Vector< MSI::Dof_Type > > & aDofTypes,
                        Vector< std::string >                  & aDofStrings,
                        mtk::Leader_Follower                             aIsLeader = mtk::Leader_Follower::LEADER );

            private:

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_QI( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * Evaluate the quantity of interest and fill aQI with value
                 * @param[ in ] aQI IQI value at evaluation point
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dQIdu( real aWStar )
                {
                    MORIS_ERROR( false, "IQI_Turbulent_Kinematic_Viscosity::compute_dQIdu - not implemented." );
                }

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
                 * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
                 */
                void compute_dQIdu(
                        Vector< MSI::Dof_Type > & aDofType,
                        Matrix< DDRMat >             & adQIdu )
                {
                    MORIS_ERROR( false, "IQI_Turbulent_Kinematic_Viscosity::compute_dQIdu() - not implemented for a drag/lift coefficient IQI.");
                }
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TURBULENT_KINEMATIC_VISCOSITY_HPP_ */

