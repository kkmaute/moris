/*
 * cl_FEM_IWG_Heat_Method_Bulk.hpp
 *
 *  Created on: Nov 04, 2020
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_HEAT_METHOD_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_HEAT_METHOD_BULK_HPP_

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

        class IWG_Heat_Method_Bulk : public IWG
        {
                //------------------------------------------------------------------------------

            private:

                // default dof types
                MSI::Dof_Type mDofTheta  = MSI::Dof_Type::THETA;
                MSI::Dof_Type mDofPhiD   = MSI::Dof_Type::PHID;

                //------------------------------------------------------------------------------
            public:

                enum class IWG_Property_Type
                {
                    MAX_ENUM
                };

                // Local string to property enum map
                std::map< std::string, IWG_Property_Type > mPropertyMap;

                enum class IWG_Constitutive_Type
                {
                        DIFFUSION,
                        MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

                enum class IWG_Stabilization_Type
                {
                        GGLS_DIFFUSION,
                        MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Heat_Method_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Heat_Method_Bulk(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_residual( real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_jacobian( real tWStar );

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

#endif /* SRC_FEM_CL_FEM_IWG_HEAT_METHOD_BULK_HPP_ */
