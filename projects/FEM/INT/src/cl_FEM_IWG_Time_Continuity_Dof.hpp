/*
 * cl_FEM_IWG_Time_Continuity_Dof.hpp
 *
 *  Created on: Apr 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_
#define SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_

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

        class IWG_Time_Continuity_Dof : public IWG
        {

//------------------------------------------------------------------------------
        public:

            enum class IWG_Property_Type
            {
                WEIGHT_CURRENT,
                WEIGHT_PREVIOUS,
                INITIAL_CONDITION,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IWG_Property_Type > mPropertyMap;

//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Time_Continuity_Dof();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Time_Continuity_Dof(){};

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
                             "IWG_Time_Continuity_Dof::set_property - Unknown aPropertyString." );

                // check no slave allowed
                MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                             "IWG_Time_Continuity_Dof::set_property - No slave allowed." );

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

#endif /* SRC_FEM_CL_FEM_IWG_TIME_CONTINUITY_DOF_HPP_ */
