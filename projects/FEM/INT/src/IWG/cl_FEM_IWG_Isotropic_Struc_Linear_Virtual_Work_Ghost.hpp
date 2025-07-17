/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_

#include "moris_typedefs.hpp"                     //MRS/COR/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost : public IWG
    {

        //------------------------------------------------------------------------------

      private:
        enum class IWG_Property_Type
        {
            THICKNESS,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            ELAST_LIN_ISO,
            MAX_ENUM
        };

        enum class IWG_Stabilization_Type
        {
            GHOST_VW,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real tWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real tWStar ) override;

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

      private:
        //------------------------------------------------------------------------------
        /**
         * method to assemble "normal matrix" from normal vector needed for
         * 2nd and 3rd order Ghost formulations
         * @param[ in ] aFlatNormal flattened normal
         * @param[ in ] aOrder      order of derivatives and ghost formulation
         */
        void get_flat_normal_matrix(
                Matrix< DDRMat > &aFlatNormal,
                uint              aOrder );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_ */

