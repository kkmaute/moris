/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_FEM_IWG_Compressible_NS_Boundary.hpp  
 * 
 */
#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_

#include <map>

#include "moris_typedefs.hpp"                         //MRS/COR/src
#include "cl_Vector.hpp"                              //MRS/CNT/src

#include "cl_Matrix.hpp"                        //LINALG/src
#include "linalg_typedefs.hpp"                  //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"        //FEM/INT/src
#include "cl_FEM_IWG_Compressible_NS_Base.hpp"  //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Compressible_NS_Boundary : public IWG_Compressible_NS_Base
    {
        //------------------------------------------------------------------------------

      private:
        // reset flags for storage variables
        bool mTractionEval    = true;
        bool mTractionDofEval = true;

        // cells of matrices containing the traction terms (K*n)
        Matrix< DDRMat > mTraction;
        Matrix< DDRMat > mTractionDOF;

        //------------------------------------------------------------------------------

      public:
        enum class IWG_Property_Type
        {
            DYNAMIC_VISCOSITY,
            THERMAL_CONDUCTIVITY,
            BODY_FORCE,
            BODY_HEAT_LOAD,
            HEAT_FLUX,
            TRACTION,
            PRESSURE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Compressible_NS_Boundary();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Compressible_NS_Boundary() override{};

        //------------------------------------------------------------------------------
        /**
         * reset eval flags specific to this IWG
         */
        void reset_child_eval_flags() override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated with evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated with evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aWStar weight associated with evaluation point
         */
        void compute_jacobian_and_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the residual wrt design variables
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dRdp( real aWStar ) override;

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /**
         * evaluate and get the Traction term ( K_ij * Y_,j * n_i )
         */
        const Matrix< DDRMat > &Traction();

        //------------------------------------------------------------------------------
        /**
         * evaluate and get the Dof-derivative of the Traction term d( K_ij * Y_,j * n_i ) / dDof
         */
        const Matrix< DDRMat > &dTractiondDOF();

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_ */
