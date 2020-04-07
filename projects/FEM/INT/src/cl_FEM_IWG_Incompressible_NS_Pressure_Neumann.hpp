/*
 * cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.hpp
 *
 *  Created on: Mar 4, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_NEUMANN_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_NEUMANN_HPP_

#include <map>
#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Incompressible_NS_Pressure_Neumann : public IWG
        {

//------------------------------------------------------------------------------
        public:

            // local property enums
            enum class IWG_Property_Type
            {
                PRESSURE,
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Incompressible_NS_Pressure_Neumann();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Incompressible_NS_Pressure_Neumann(){};

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
                // check that aPropertyString makes sense
                MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                             "IWG_Incompressible_NS_Pressure_Neumann::set_property - Unknown aPropertyString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Incompressible_NS_Pressure_Neumann::set_property - No slave allowed." );

                // set the property in the property cell
                this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            }

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian_and_residual( real aWStar );

//------------------------------------------------------------------------------
            /**
             * compute the derivative of the residual wrt design variables
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_dRdp( real aWStar );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */




#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_BULK_HPP_ */
