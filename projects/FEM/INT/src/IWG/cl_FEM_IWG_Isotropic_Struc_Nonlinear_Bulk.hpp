/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_BULK_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Nonlinear_Bulk : public IWG
        {
        	//------------------------------------------------------------------------------

        	private:
            // stress and strain type to evaluate the IWG
            enum CM_Function_Type mStressType;
            enum CM_Function_Type mStrainType;

                //------------------------------------------------------------------------------
            public:

                enum class IWG_Property_Type
                {
                    LOAD,
                    BEDDING,
                    THICKNESS,
                    MAX_ENUM
                };

                enum class IWG_Constitutive_Type
                {
                    ELAST_LIN_ISO,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IWG_Isotropic_Struc_Nonlinear_Bulk(
                		enum CM_Function_Type aStressType,
        				enum CM_Function_Type aStrainType );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Isotropic_Struc_Nonlinear_Bulk(){};

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

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_NONLINEAR_BULK_HPP_ */

