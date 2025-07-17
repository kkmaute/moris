/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_L2_Damage_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_L2_Damage_Bulk : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        // source type to evaluate
        uint mSourceType = 0;

        // FIXME for now only implemented for the specific case of equivalent strain
        // parameters
        real             mCharacteristicLength = 0.0;
        uint             mOrder                = 1;
        Matrix< DDRMat > mOrderCoeff           = { { 2.0 }, { 8.0 }, { 48.0 } };

        enum class IWG_Property_Type
        {
            WEIGHT,
            LUMP,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            ELASTIC_DAMAGE,
            MAX_ENUM
        };

        enum class IWG_Stabilization_Type
        {
            MAX_ENUM
        };

        // pointer to implementation for different source
        const Matrix< DDRMat >& ( IWG_L2_Damage_Bulk::*m_get_source )() = nullptr;
        const Matrix< DDRMat >& ( IWG_L2_Damage_Bulk::*m_get_dsourcedu )(
                const moris::Vector< MSI::Dof_Type >& aDofType ) = nullptr;

        //------------------------------------------------------------------------------
        /*
         *  constructor
         * @param[ in ] aSourceType uint for source type
         * 0 - equivalent strain
         * 1 - history
         * 2 - smooth damage
         */
        IWG_L2_Damage_Bulk( uint aSourceType );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_L2_Damage_Bulk() override{};

        //------------------------------------------------------------------------------
        /**
         * set parameters
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

      private:
        //------------------------------------------------------------------------------
        /**
         * get source contribution
         */
        const Matrix< DDRMat >& get_source();

        const Matrix< DDRMat >& get_source_eqStrain();
        const Matrix< DDRMat >& get_source_smoothDam();
        const Matrix< DDRMat >& get_source_history();

        //------------------------------------------------------------------------------
        /**
         * get source derivative
         */
        const Matrix< DDRMat >&
        get_dsourcedu( const moris::Vector< MSI::Dof_Type >& aDofType );

        const Matrix< DDRMat >&
        get_dsourcedu_eqStrain( const moris::Vector< MSI::Dof_Type >& aDofType );
        const Matrix< DDRMat >&
        get_dsourcedu_smoothDam( const moris::Vector< MSI::Dof_Type >& aDofType );
        const Matrix< DDRMat >&
        get_dsourcedu_history( const moris::Vector< MSI::Dof_Type >& aDofType );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_L2_DAMAGE_BULK_HPP_ */
