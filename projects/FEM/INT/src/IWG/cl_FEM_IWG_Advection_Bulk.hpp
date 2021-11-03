/*
 * cl_FEM_IWG_Advection_Bulk.hpp
 *
 *  Created on: Apr 12, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_ADVECTION_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ADVECTION_BULK_HPP_

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

        class IWG_Advection_Bulk : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // local property enums
                enum class IWG_Property_Type
                {
                    BODY_LOAD,
                    MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                    DIFFUSION,
                    MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                    SUPG,
                    YZBETA,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Advection_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Advection_Bulk(){};

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

            private:
                //------------------------------------------------------------------------------
                /**
                 * compute the residual strong form
                 * @param[ in ] aRT a matrix to fill with RM
                 */
                void compute_residual_strong_form( Matrix< DDRMat > & aRT );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual strong form
                 * @param[ in ] aDofTypes a list of dof type wrt which
                 *                        the derivative is requested
                 * @param[ in ] aJT       a matrix to fill with dRMdDof
                 */
                void compute_jacobian_strong_form(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aJT );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ADVECTION_BULK_HPP_ */
