/*
 * cl_FEM_IWG_Incompressible_NS_Viscous_Velocity_Ghost.hpp
 *
 *  Created on: Mar 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VISCOUS_VELOCITY_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VISCOUS_VELOCITY_GHOST_HPP_

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

        class IWG_Incompressible_NS_Viscous_Velocity_Ghost : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // interpolation order for residual dof type
                uint mOrder;

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                        VISCOUS_GHOST,
                        TIME_VELOCITY_GHOST,
                        MAX_ENUM
                };

                // local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Incompressible_NS_Viscous_Velocity_Ghost();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Incompressible_NS_Viscous_Velocity_Ghost(){};

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
            private:

                //------------------------------------------------------------------------------
                /**
                 * compute the block matrix for directional derivative
                 * (dnNdxn . normal)
                 * @param[ in ] aFlatdnNdxn matrix to fill with derivatives
                 * @param[ in ] aOrder      interpolation order for residual dof type
                 * @param[ in ] aIsMaster   enum for master or slave
                 */
                void compute_flat_dnNdxn(
                        Matrix< DDRMat >  & aFlatdnNdxn,
                        uint                aOrder,
                        mtk::Master_Slave   aIsMaster );

                //------------------------------------------------------------------------------
                /**
                 * compute flattened normal matrix
                 * @param[ in ] aFlatNormal matrix to fill with flattened normal
                 * @param[ in ] aOrder      interpolation order for residual dof type
                 */
                void get_normal_matrix(
                        Matrix< DDRMat > & aFlatNormal,
                        uint               aOrder );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VISCOUS_VELOCITY_GHOST_HPP_ */
