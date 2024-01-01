/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_MM_Perfect_Gas.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_MM_PERFECT_GAS_HPP_
#define SRC_FEM_CL_FEM_MM_PERFECT_GAS_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Material_Model.hpp"        //FEM/INT/src

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------

        class MM_Perfect_Gas : public Material_Model
        {

                //------------------------------------------------------------------------------
            protected:
                // default local properties
                std::shared_ptr< Property > mPropIsochoricHeatCapacity  = nullptr;
                std::shared_ptr< Property > mPropSpecificGasConstant    = nullptr;

            private:

                // property type for MM
                enum class MM_Property_Type
                {
                        ISOCHORIC_HEAT_CAPACITY,  // heat capacity at constant density
                        SPECIFIC_GAS_CONSTANT,    // specific gas constant for fluid (R = R_0/M)
                        MAX_ENUM
                };

                // thermodynamic variables as protected DoF types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                MM_Perfect_Gas();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~MM_Perfect_Gas(){};

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dof types
                 * @param[ in ] aDofTypes a list of group of dof types
                 * @param[ in ] aDofStrings a list of strings to describe the dof types
                 */
                void set_dof_type_list(
                        Vector< Vector< MSI::Dof_Type > > aDofTypes,
                        Vector< std::string >                  aDofStrings );

                //------------------------------------------------------------------------------
                /**
                 * set local properties
                 */
                void set_local_properties();

                //------------------------------------------------------------------------------
                // SPECIFIC INTERNAL ENERGY (FIRST EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the specific internal energy and its derivatives
                 */
                void eval_Eint();
                void eval_EintDot();
                void eval_dEintdx();
                void eval_d2Eintdx2();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the specific internal energy derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_EintDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_EintDotDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_dEintdxDOF( const Vector< MSI::Dof_Type > & aDofTypes);
                void eval_d2Eintdx2DOF( const Vector< MSI::Dof_Type > & aDofTypes);

                //------------------------------------------------------------------------------
                // DENSITY (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the density and its derivatives
                 */
                void eval_density();
                void eval_DensityDot();
                void eval_dDensitydx();
                void eval_d2Densitydx2();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic density derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_DensityDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_DensityDotDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_dDensitydxDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_d2Densitydx2DOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                // PRESSURE (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the pressure and its derivatives
                 */
                void eval_pressure();
                void eval_PressureDot();
                void eval_dPressuredx();
                void eval_d2Pressuredx2();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic pressure derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_PressureDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_PressureDotDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_dPressuredxDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_d2Pressuredx2DOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                // TEMPERATURE (SECOND EQUATION OF STATE)
                //------------------------------------------------------------------------------
                /**
                 * evaluate the temperature and its derivatives
                 */
                void eval_temperature();
                void eval_TemperatureDot();
                void eval_dTemperaturedx();
                void eval_d2Temperaturedx2();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic temperature derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_TemperatureDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_TemperatureDotDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_dTemperaturedxDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_d2Temperaturedx2DOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                // THERMODYNAMIC QUANTITIES
                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic quantities
                 */
                void eval_VolumeExpansivity();
                void eval_IsothermalCompressibility();
                void eval_Cv();
                void eval_Cp();
                void eval_Gamma();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermodynamic quantity derivatives wrt to the dof types
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_VolumeExpansivityDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_IsothermalCompressibilityDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_CvDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_CpDOF( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_GammaDOF( const Vector< MSI::Dof_Type > & aDofTypes );
        };

        //--------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_MM_PERFECT_GAS_HPP_ */

