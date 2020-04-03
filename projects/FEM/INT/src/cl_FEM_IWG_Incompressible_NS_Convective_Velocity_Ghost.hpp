/*
 * cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp
 *
 *  Created on: Mar 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_

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

        class IWG_Incompressible_NS_Convective_Velocity_Ghost : public IWG
        {

//------------------------------------------------------------------------------
        public:

            // local stabilization enums
            enum class IWG_Stabilization_Type
            {
                CONVECTIVE_GHOST,
                MAX_ENUM
            };

            // local string to constitutive enum map
            std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Incompressible_NS_Convective_Velocity_Ghost();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Incompressible_NS_Convective_Velocity_Ghost(){};

//------------------------------------------------------------------------------
            /**
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                              std::string                                aStabilizationString )
            {
                // check that aConstitutiveString makes sense
                MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                             "IWG_Incompressible_NS_Convective_Velocity_Ghost::set_stabilization_parameter - Unknown aStabilizationString." );

                // set the stabilization parameter in the stabilization parameter cell
                this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
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
        private:

//------------------------------------------------------------------------------
            /**
             * compute the block matrix for dnNdxn
             * @param[ in ] adNdx     matrix to fill with derivatives
             * @param[ in ] aIsMaster enum for master or slave
             */
            void compute_dnNdxn( Matrix< DDRMat >  & adNdx,
                                 mtk::Master_Slave   aIsMaster );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST_HPP_ */
