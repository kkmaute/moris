/*
 * cl_FEM_IQI_Strain_Energy.hpp
 *
 *  Created on: Feb, 2020
 *      Author: doble
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_

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

        class IQI_Stress : public IQI
        {

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------

                // stress type to evaluate
                enum Stress_Type mStressType;

                enum class IQI_Property_Type
                {
                        MAX_ENUM
                };

                // Local string to property enum map
                std::map< std::string, IQI_Property_Type > mPropertyMap;

                enum class IQI_Constitutive_Type
                {
                        ELAST_LIN_ISO,
                        MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IQI_Constitutive_Type > mConstitutiveMap;

                enum class IQI_Stabilization_Type
                {
                        MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IQI_Stabilization_Type > mStabilizationMap;

                //------------------------------------------------------------------------------

            protected:

                //------------------------------------------------------------------------------
                /*
                 * default constructor
                 */
                IQI_Stress(){
                    MORIS_ERROR( false, " IQI_Stress: Default constructor unavailable. Use enum to construct IQI. ");
                };

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Stress( enum Stress_Type aStressType );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Stress(){};

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property( std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
                {
                    // can only be master
                    MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                            "IQI::set_property - can only be master." );

                    // FIXME check that property type makes sense?

                    // set the property in the property cell
                    this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
                }

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model
                 * @param[ in ] aConstitutiveModel  a constitutive model pointer
                 * @param[ in ] aConstitutiveString a string defining the constitutive model
                 * @param[ in ] aIsMaster           an enum for master or slave
                 */
                void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                        std::string                           aConstitutiveString,
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
                {
                    // can only be master
                    MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                            "IQI::set_constitutive model - can only be master." );

                    // FIXME check that constitutive string makes sense?

                    // set the constitutive model in the constitutive model cell
                    this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
                }

                //------------------------------------------------------------------------------
                /**
                 * set stabilization parameter
                 * @param[ in ] aStabilizationParameter a stabilization parameter pointer
                 * @param[ in ] aStabilizationString    a string defining the stabilization parameter
                 */
                void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                        std::string                                aStabilizationString )
                {
                    // FIXME check that stabilization string makes sense?

                    // set the stabilization parameter in the stabilization parameter cell
                    this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
                }

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aQI quantity of interest matrix to fill
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] adQIdDof derivative of quantity of interest matrix to fill
                 */
                void compute_dQIdDof( Matrix< DDRMat > & adQIdDof );

                //------------------------------------------------------------------------------
                /**
                 * gets the stress vector from the constitutive model and sorts it into a format standardized
                 * (independent of 2D-plane strain, 2D-plain stress, and 3D)
                 * @param[ out ] tStressVector vector with 6 entries containing the normal and shear stress values
                 */
                Matrix< DDRMat > get_stress_vector();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the Von-Mises stress using the stress vector provided by the constitutive model
                 * @param[ out ] tStressValue the value of the Von-Mises stress
                 */
                real eval_Von_Mises_stress();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th principal stress using the stress vector provided by the constitutive model
                 * @param[ in ] aPrincipalStressIndex number i of the principal stress to get,
                 *               i.e. 1,2,3 for the first, second and third principal stress, respectively
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_principal_stress( uint aPrincipalStressIndex );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th normal stress using the stress vector provided by the constitutive model
                 * @param[ in ] aStressIndex number of the principal stress to get,
                 *               i.e. 1,2,3 for the x-, y-, and z-normal stresses, respectively
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_normal_stress( uint aStressIndex );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th shear stress using the stress vector provided by the constitutive model
                 * @param[ in ] aStressIndex number of the principal stress to get,
                 *               i.e. 1,2,3 for the yz-, xz-, and xy-shear stresses, respectively
                 *               in 2D the input is ignored, as there's only one shear stress
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_shear_stress( uint aStressIndex );

                //------------------------------------------------------------------------------

        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_ */
