/*
 * cl_FEM_IQI_Volume.hpp
 *
 *  Created on: Nov 20, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_

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

        class IQI_Volume : public IQI
        {
            enum class IQI_Property_Type
            {
                DENSITY,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IQI_Property_Type > mPropertyMap;

            //------------------------------------------------------------------------------
        public:

            //------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IQI_Volume()
            {
                // set size for the property pointer cell
                mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

                // populate the property map
                mPropertyMap[ "Density" ] = IQI_Property_Type::DENSITY;
            }

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Volume(){};

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

        private:

            /**
             * compute the quantity of interest
             * @param[ in ] aQI quantity of interest matrix to fill
             */
            void compute_QI( Matrix< DDRMat > & aQI );

            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] adQIdu derivative of quantity of interest matrix to fill
             */
            void compute_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu );

            //------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_ */
