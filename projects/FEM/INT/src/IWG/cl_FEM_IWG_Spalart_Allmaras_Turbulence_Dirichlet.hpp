/*
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp
 *
 *  Created on: May 27, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_
//MRS/COR/src
#include <map>
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG.hpp"
#include "fn_FEM_FD_Scheme.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Spalart_Allmaras_Turbulence_Dirichlet : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // sint for symmetric/unsymmetric Nitsche formulation
                sint mBeta = 1.0;

                // local property enums
                enum class IWG_Property_Type
                {
                        DIRICHLET,
                        VISCOSITY,
                        SELECT,
                        UPWIND,
                        MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                        NITSCHE,
                        MAX_ENUM
                };

            private:

//                // Spalart Allmaras model constants
//                real mCb1 = 0.1355;
                real mCb2 = 0.6220;
                real mSigma = 2.0/3.0;
//                real mKappa = 0.41;
//                real mCw1 = mCb1 / std::pow( mKappa, 2.0 ) + ( 1.0 + mCb2 ) / mSigma;
//                real mCw2 = 0.3;
//                real mCw3 = 2.0;
//                real mCt3 = 1.2;
//                real mCt4 = 0.5;
//                real mCv1 = 7.1;
//                real mCv2 = 0.7;
//                real mCv3 = 0.9;
//                real mRLim = 10.0;
//                real mCn1 = 16.0;

            public:
                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 *  aBeta signed int for symmetric/unsymmetric Nitsche
                 */
                IWG_Spalart_Allmaras_Turbulence_Dirichlet( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Spalart_Allmaras_Turbulence_Dirichlet(){};

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
//                //------------------------------------------------------------------------------
//                /**
//                 * compute diffusion coefficient
//                 * Diff = kinViscosity + modViscosity
//                 * if modViscosity >= 0
//                 * Diff = kinViscosity + modViscosity * fn
//                 * if modViscosity <  0
//                 */
//                real compute_diffusion_coefficient();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of the diffusion coefficient
//                 * @param[ in ] aDofTypes     a list of dof type wrt which
//                 *                            the derivative is requested
//                 * @param[ in ] addiffusiondu a matrix to fill with ddiffusiondu
//                 */
//                void compute_ddiffusiondu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & addiffusiondu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute fn
//                 * fn = ( cn1 + chi³ ) / ( cn1 - chi³)
//                 */
//                real compute_fn();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fn wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfndu     a matrix to fill with dfndu
//                 */
//                void compute_dfndu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfndu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute chi = viscosityDof / viscosityPtop
//                 * @param[ out ] chi
//                 */
//                real compute_chi();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of chi wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adchidu    a matrix to fill with dchidu
//                 */
//                void compute_dchidu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adchidu );

                //------------------------------------------------------------------------------
                /**
                 * compute the traction = ( v + vtilde ) * grad vtilde / mSigma
                 * @param[ in ] aTraction a matrix to fill with traction
                 */
                void compute_traction( Matrix< DDRMat > & aTraction );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the traction = ( v + vtilde ) * grad vtilde / mSigma
                 * wrt dof type aDofTypes
                 * @param[ in ] aDofTypes    group of dervative dof types
                 * @param[ in ] adtractiondu a matrix to fill with dtractiondu
                 */
                void compute_dtractiondu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adtractiondu );

                //------------------------------------------------------------------------------
                /**
                 * compute the test traction
                 *  = delta ( ( v + vtilde ) * grad vtilde  / mSigma )
                 * @param[ in ] aTestDofTypes group of test dof types
                 * @param[ in ] aTestTraction a matrix to fill with test traction
                 */
                void compute_testtraction(
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        Matrix< DDRMat >                   & aTestTraction );

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
                        Matrix< DDRMat >                  & adtesttractiondu );

                void compute_dtesttractiondu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        Matrix< DDRMat >                   & adtesttractiondu_FD,
                        real                                 aPerturbation = 1e-6,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_DIRICHLET_HPP_ */
