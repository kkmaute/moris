/*
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp
 *
 *  Created on: Mar 23, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_
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

        class IWG_Spalart_Allmaras_Turbulence_Bulk : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // local property enums
                enum class IWG_Property_Type
                {
                    WALL_DISTANCE,
                    VISCOSITY,
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, IWG_Property_Type > mPropertyMap;

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                    SUPG,
                    MAX_ENUM
                };

                // local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                // FIXME temp all the constants
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

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Spalart_Allmaras_Turbulence_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Spalart_Allmaras_Turbulence_Bulk(){};

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
                 * compute the residual strong form
                 * @param[ in ] aR real to fill with R
                 */
                void compute_residual_strong_form( Matrix< DDRMat > & aR );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian strong form
                 * @param[ in ] aDofTypes a list of dof type wrt which
                 *                        the derivative is requested
                 * @param[ in ] aJ        a matrix to fill with dRdDof
                 */
                void compute_jacobian_strong_form(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & aJ );

                //------------------------------------------------------------------------------
                /**
                 * compute Wij = 0.5 * ( dui/dxj - duj/dxi )
                 * @param[ out ] Wij
                 */
                void compute_divflux( Matrix< DDRMat > & aDivFlux );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of Wij wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adwijdu    a matrix to fill with dwijdu
                 */
                void compute_ddivfluxdu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adDivFluxdu );

                //------------------------------------------------------------------------------
                /**
                 * compute Wij = 0.5 * ( dui/dxj - duj/dxi )
                 * @param[ out ] Wij
                 */
                void compute_wij( Matrix< DDRMat > & aWij );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of Wij wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adwijdu    a matrix to fill with dwijdu
                 */
                void compute_dwijdu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adwijdu );

                //------------------------------------------------------------------------------
                /**
                 * compute chi = viscosityDof / viscosityPtop
                 * @param[ out ] chi
                 */
                real compute_chi();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidu    a matrix to fill with dchidu
                 */
                void compute_dchidu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adchidu );

                //------------------------------------------------------------------------------
                /**
                 * compute fv1 = chi³ / ( chi³ + cv1³)
                 * @param[ out ] fv1
                 */
                real compute_fv1();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv1 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
                 */
                void compute_dfv1du(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adfv1du );

                //------------------------------------------------------------------------------
                /**
                 * compute fv2 = 1 - chi / ( 1 + chi * fv1 )
                 * @param[ out ] fv2
                 */
                real compute_fv2();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv2 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv2du    a matrix to fill with dfv2du
                 */
                void compute_dfv2du(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adfv2du );

                //------------------------------------------------------------------------------
                /**
                 * compute SBar = sqrt( 2 * Wij * Wij )
                 * @param[ out ] SBar
                 */
                real compute_sbar();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of SBar wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adSBardu a matrix to fill with dSBardu
                 */
                void compute_dsbardu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adsbardu );

                //------------------------------------------------------------------------------
                /**
                 * compute STilde = SBar + viscosity * fv2 / ( kappa² * d² )
                 * @param[ out ] STilde
                 */
                real compute_stilde();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of STilde wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adSTildedu a matrix to fill with dSTildedu
                 */
                void compute_dstildedu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adstildedu );

                //------------------------------------------------------------------------------
                /**
                 * compute ft2 = ct3 * exp( -ct4 * chi²), chi = viscosityDof/viscosityProp
                 * @param[ out ] ft2
                 */
                real compute_ft2();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of ft2 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adft2du a matrix to fill with dft2du
                 */
                void compute_dft2du(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >          & adft2du );

                //------------------------------------------------------------------------------
                /**
                 * compute r = min( viscosity / ( STilde * kappa² * d² ), 10 )
                 * @param[ out ] r
                 */
                real compute_r();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of r wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adrdu a matrix to fill with drdu
                 */
                void compute_drdu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adrdu );

                //------------------------------------------------------------------------------
                /**
                 * compute g = r + cw2 * ( r⁶ - r )
                 * @param[ out ] g
                 */
                real compute_g();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of g wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adgdu a matrix to fill with dgdu
                 */
                void compute_dgdu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adgdu );

                //------------------------------------------------------------------------------
                /**
                 * compute fw = g * ( ( 1 + cw3⁶ ) / ( g⁶ + cw3⁶ ) )^(1/6)
                 * g = r + cw2 * ( r⁶ - r )
                 * r = min( viscosity / ( STilde * kappa² * d² ), 10 )
                 * @param[ out ] fw
                 */
                real compute_fw();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fw wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfwdu a matrix to fill with dfwdu
                 */
                void compute_dfwdu(
                        moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >             & adfwdu );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_BULK_HPP_ */
