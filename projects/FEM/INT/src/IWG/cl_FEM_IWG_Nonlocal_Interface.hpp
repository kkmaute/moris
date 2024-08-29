/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Nonlocal_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_NONLOCAL_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_NONLOCAL_INTERFACE_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Nonlocal_Interface : public IWG
    {

      public:
        // sint for symmetric/unsymmetric Nitsche
        sint mBeta = 1.0;

        real             mCharacteristicLength = 1.0;
        uint             mOrder                = 1;
        Matrix< DDRMat > mOrderCoeff           = { { 2.0 }, { 8.0 }, { 48.0 } };

        enum class IWG_Property_Type
        {
            WEIGHT,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            MAX_ENUM
        };

        enum class IWG_Stabilization_Type
        {
            NITSCHE_INTERFACE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Nonlocal_Interface( sint aBeta );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Nonlocal_Interface() override{};

        //------------------------------------------------------------------------------
        /**
         * set parameters
         * @param[ in ] aParameters a list of parameters
         */
        void set_parameters( const Vector< Matrix< DDRMat > >& aParameters ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian_and_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the residual wrt design variables
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dRdp( real aWStar ) override;

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_NONLOCAL_INTERFACE_HPP_ */
