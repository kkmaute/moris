/*
 * cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.hpp
 *
 *  Created on: Mar 03, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_BULK_HPP_

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

        class IWG_Incompressible_NS_Velocity_Bulk : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // local property enums
                enum class IWG_Property_Type
                {
                    GRAVITY,
                    THERMAL_EXPANSION,
                    REF_TEMP,
                    INV_PERMEABILITY,
                    MASS_SOURCE,
                    BODY_LOAD,
                    HOMOTOPY,
                    MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                    INCOMPRESSIBLE_FLUID,
                    MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                    INCOMPRESSIBLE_FLOW,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Incompressible_NS_Velocity_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Incompressible_NS_Velocity_Bulk(){};

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
                 * @param[ in ] aRM a matrix to fill with RM
                 * @param[ in ] aRC a matrix to fill with RC
                 */
                void compute_residual_strong_form(
                        Matrix< DDRMat > & aRM,
                        real             & aRC );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual strong form
                 * @param[ in ] aDofTypes a list of dof type wrt which
                 *                        the derivative is requested
                 * @param[ in ] aJM       a matrix to fill with dRMdDof
                 * @param[ in ] aJC       a matrix to fill with dRCdDof
                 */
                void compute_jacobian_strong_form(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & aJM,
                        Matrix< DDRMat >                   & aJC );

                //------------------------------------------------------------------------------
                /**
                 * compute the term uj vij
                 * @param[ in ] aujvij a matrix to fill with uj vij
                 */
                void compute_ujvij( Matrix< DDRMat > & aujvij );

                //------------------------------------------------------------------------------
                /**
                 * compute the term uj vij rm
                 * @param[ in ] aujvijrm a matrix to fill with uj vij rm
                 * @param[ in ] arm      provided strong form residual
                 */
                void compute_ujvijrm(
                        Matrix< DDRMat > & aujvijrm,
                        Matrix< DDRMat > & arm );

                //------------------------------------------------------------------------------
                // FIXME provided directly by the field interpolator?
                /**
                 * compute the term dnNdtn
                 * @param[ in ] adnNdtn a matrix to fill with dnNdtn
                 */
                void compute_dnNdtn( Matrix< DDRMat > & adnNdtn );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_BULK_HPP_ */
