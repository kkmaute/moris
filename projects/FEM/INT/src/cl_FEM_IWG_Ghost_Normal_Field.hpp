/*
 * cl_FEM_IWG_Ghost_Normal_Field.hpp
 *
 *  Created on: Jul 27, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_
#define SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_IWG.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Ghost_Normal_Field : public IWG
        {
                //------------------------------------------------------------------------------
            private:

                enum class IWG_Stabilization_Type
                {
                    GHOST_SP,
                    MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

            public:
                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Ghost_Normal_Field();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Ghost_Normal_Field(){};

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
                 * @param[ in ] aResidual cell of residual vectors to fill
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
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_Ghost_Normal_Field_HPP_ */
