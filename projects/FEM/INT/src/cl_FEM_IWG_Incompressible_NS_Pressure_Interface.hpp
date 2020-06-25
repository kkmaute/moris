/*
 * cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp
 *
 *  Created on: Jun 02, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Incompressible_NS_Pressure_Interface : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // sign for symmetric/unsymmetric Nitsche
                sint mBeta = 1;

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                    FLUID_INCOMPRESSIBLE,
                    MAX_ENUM
                };

                // local string to constitutive enum map
                std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

                enum class IWG_Stabilization_Type
                {
                    MASTER_WEIGHT_INTERFACE,
                    SLAVE_WEIGHT_INTERFACE,
                    MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * @param[ in ] aBeta +1 or -1 for symmetric/skew symmetric Nitsche
                 */
                IWG_Incompressible_NS_Pressure_Interface( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Incompressible_NS_Pressure_Interface(){};

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
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER );


                //------------------------------------------------------------------------------
                /**
                 * set stabilization parameter
                 * @param[ in ] aStabilizationParameter a stabilization parameter pointer
                 * @param[ in ] aStabilizationString    a string defining the stabilization parameter
                 */
                void set_stabilization_parameter(
                        std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                        std::string                                aStabilizationString );

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

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_HPP_ */
