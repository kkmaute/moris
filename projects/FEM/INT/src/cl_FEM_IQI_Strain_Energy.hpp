/*
 * cl_FEM_IQI_Strain_Energy.hpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_

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

        class IQI_Strain_Energy : public IQI
        {
                //------------------------------------------------------------------------------

                enum class IQI_Constitutive_Type
                {
                    ELAST,
                    MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IQI_Constitutive_Type > mConstitutiveMap;

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Strain_Energy();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Strain_Energy(){};

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model
                 * @param[ in ] aConstitutiveModel  a constitutive model pointer
                 * @param[ in ] aConstitutiveString a string defining the constitutive model
                 * @param[ in ] aIsMaster           an enum for master or slave
                 */
                void set_constitutive_model(
                        std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                        std::string                           aConstitutiveString,
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
                {
                    // check that aConstitutiveString makes sense
                    MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                            "IQI_Strain_Energy::set_constitutive_model - Unknown aConstitutiveString." );

                    // check no slave allowed
                    MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                            "IQI_Strain_Energy::set_constitutive_model - No slave allowed." );

                    // set the constitutive model in the constitutive model cell
                    this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
                }

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
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_ */
