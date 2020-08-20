/*
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp
 *
 *  Created on: Jun 02, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_
//MRS/COR/src
#include <map>
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

        class IWG_Spalart_Allmaras_Turbulence_Interface : public IWG
        {

                //------------------------------------------------------------------------------
            private:

                // sint for symmetric/unsymmetric Nitsche formulation
                sint mBeta = 1.0;

                // local property enums
                enum class IWG_Property_Type
                {
                    VISCOSITY,
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, IWG_Property_Type > mPropertyMap;

                // local stabilization parameter enums
                enum class IWG_Stabilization_Type
                {
                        NITSCHE_INTERFACE,
                        MAX_ENUM
                };

                // local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                // Spalart Allmaras model constants
                real mCb1 = 0.1355;
                real mCb2 = 0.6220;
                real mSigma = 2.0/3.0;
                real mKappa = 0.41;
                real mCw1 = mCb1 / std::pow( mKappa, 2.0 ) + ( 1.0 + mCb2 ) / mSigma;
                real mCw2 = 0.3;
                real mCw3 = 2.0;
                real mCt3 = 1.2;
                real mCt4 = 0.5;
                real mCv1 = 7.1;
                real mCv2 = 0.7;
                real mCv3 = 0.9;
                real mRLim = 10.0;
                real mCn1 = 16.0;

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * @param[ in ] aBeta sint for symmetric/unsymmetric Nitsche formulation
                 */
                IWG_Spalart_Allmaras_Turbulence_Interface( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Spalart_Allmaras_Turbulence_Interface(){};

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

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
                 * compute diffusion coefficient
                 * Diff = kinViscosity + modViscosity
                 * if modViscosity >= 0
                 * Diff = kinViscosity + modViscosity * fn
                 * if modViscosity <  0
                 */
                real compute_diffusion_coefficient(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the diffusion coefficient
                 * @param[ in ] aDofTypes     a list of dof type wrt which
                 *                            the derivative is requested
                 * @param[ in ] addiffusiondu a matrix to fill with ddiffusiondu
                 */
                void compute_ddiffusiondu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & addiffusiondu,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute fn
                 * fn = ( cn1 + chi³ ) / ( cn1 - chi³)
                 */
                real compute_fn( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fn wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfndu     a matrix to fill with dfndu
                 */
                void compute_dfndu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adfndu,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute chi = viscosityDof / viscosityPtop
                 * @param[ out ] chi
                 */
                real compute_chi( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidu    a matrix to fill with dchidu
                 */
                void compute_dchidu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adchidu,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the traction = ( v + vtilde ) * grad vtilde / mSigma
                 * @param[ in ] aTraction a matrix to fill with traction
                 */
                void compute_traction(
                        Matrix< DDRMat >  & aTraction,
                        mtk::Master_Slave   aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the traction = ( v + vtilde ) * grad vtilde / mSigma
                 * wrt dof type aDofTypes
                 * @param[ in ] aDofTypes    group of dervative dof types
                 * @param[ in ] adtractiondu a matrix to fill with dtractiondu
                 */
                void compute_dtractiondu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adtractiondu,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the test traction
                 *  = delta ( ( v + vtilde ) * grad vtilde  / mSigma )
                 * @param[ in ] aTestDofTypes group of test dof types
                 * @param[ in ] aTestTraction a matrix to fill with test traction
                 */
                void compute_testtraction(
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        Matrix< DDRMat >                   & aTestTraction,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the test traction
                 * = delta ( ( v + vtilde ) * grad vtilde  / mSigma )
                 * @param[ in ] aDofTypes        group of derivative dof types
                 * @param[ in ] aTestDofTypes    group of test dof types
                 * @param[ in ] adtesttractiondu a matrix to fill with test traction
                 */
                void compute_dtesttractiondu(
                        const moris::Cell< MSI::Dof_Type> & aDofTypes,
                        const moris::Cell< MSI::Dof_Type> & aTestDofTypes,
                        Matrix< DDRMat >                  & adtesttractiondu,
                        mtk::Master_Slave                   aIsMaster = mtk::Master_Slave::MASTER );

                void compute_dtesttractiondu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        Matrix< DDRMat >                   & adtesttractiondu_FD,
                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER,
                        real                                 aPerturbation = 1e-6,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_ */
