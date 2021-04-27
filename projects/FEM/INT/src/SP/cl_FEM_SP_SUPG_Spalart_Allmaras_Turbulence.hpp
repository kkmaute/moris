/*
 * cl_FEM_SP_Spalart_Allmaras_Turbulence.hpp
 *
 *  Created on: May 19, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Cluster.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_SUPG_Spalart_Allmaras_Turbulence : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private :

                // default tuple for element size to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Master_Slave > mElementSizeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_LENGTH_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Master_Slave::MASTER );

                // default dof type
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;
                MSI::Dof_Type mMasterDofVelocity  = MSI::Dof_Type::VX;

                // Spalart-Allmaras model constants
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

                // internal threshold
                const real mEpsilon = 1e-12;//1e-18;

                // property type for the SP
                enum class SP_Property_Type
                {
                    VISCOSITY, // fluid viscosity
                    WALL_DISTANCE,
                    MAX_ENUM
                };

                /*
                 * Rem: mParameters( 0 ) -
                 */

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_SUPG_Spalart_Allmaras_Turbulence();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_SUPG_Spalart_Allmaras_Turbulence(){};

                //------------------------------------------------------------------------------
                /**
                 * set function pointers for evaluation
                 */
                void set_function_pointers();

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                        moris::Cell< std::string >             & aDvStrings,
                        mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsMaster );
                }

                //------------------------------------------------------------------------------
                /**
                 * get cluster measure tuples
                 * @param[ in ] aClusterMeasureTuples list of tuples describing the cluster measure types
                 */
                moris::Cell< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Master_Slave > > get_cluster_measure_tuple_list();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter value
                 */
                void eval_SP();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 */
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdMasterDV - not implemented." );
                }

                //------------------------------------------------------------------------------
            private:
//                //------------------------------------------------------------------------------
//                /**
//                 * compute production term
//                 * P = cb1 * ( 1 - ft2 ) * STilde * modViscosity
//                 * if modViscosity >= 0
//                 * P = cb1 * ( 1 - ct3 ) * S * modViscosity
//                 * if modViscosity <  0
//                 */
//                real compute_production_coefficient();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of the production term
//                 * @param[ in ] aDofTypes      a list of dof type wrt which
//                 *                             the derivative is requested
//                 * @param[ in ] adproductiondu a matrix to fill with dproductiondu
//                 */
//                void compute_dproductiondu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adproductiondu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute wall destruction term
//                 * D = ( cw1 * fw - cb1 * ft2 / kappa² ) * ( modViscosity / wallD )²
//                 * if modViscosity >= 0
//                 * D = - cw1 * ( modViscosity / wallD )²
//                 * if modViscosity <  0
//                 */
//                real compute_wall_destruction_coefficient();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of the wall destruction term
//                 * @param[ in ] aDofTypes           a list of dof type wrt which
//                 *                                  the derivative is requested
//                 * @param[ in ] adwalldestructiondu a matrix to fill with dwalldestructiondu
//                 */
//                void compute_dwalldestructiondu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adwalldestructiondu );
//
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
//                 * compute Wij = 0.5 * ( dui/dxj - duj/dxi )
//                 * @param[ out ] Wij
//                 */
//                void compute_wij( Matrix< DDRMat > & aWij );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of Wij wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adwijdu    a matrix to fill with dwijdu
//                 */
//                void compute_dwijdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adwijdu );
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
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute fv1 = chi³ / ( chi³ + cv1³)
//                 * @param[ out ] fv1
//                 */
//                real compute_fv1();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fv1 wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
//                 */
//                void compute_dfv1du(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfv1du );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute fv2 = 1 - chi / ( 1 + chi * fv1 )
//                 * @param[ out ] fv2
//                 */
//                real compute_fv2();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fv2 wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfv2du    a matrix to fill with dfv2du
//                 */
//                void compute_dfv2du(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfv2du );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute SBar = sqrt( 2 * Wij * Wij )
//                 * @param[ out ] SBar
//                 */
//                real compute_sbar();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of SBar wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adSBardu a matrix to fill with dSBardu
//                 */
//                void compute_dsbardu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adsbardu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute S = viscosity * fv2 / ( kappa² * d² )
//                 * @param[ out ] S
//                 */
//                real compute_s();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of S wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adSdu a matrix to fill with dSdu
//                 */
//                void compute_dsdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >             & adsdu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute Smod =
//                 * @param[ out ] Smod
//                 */
//                real compute_smod();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of Smod wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adSmoddu a matrix to fill with dSmoddu
//                 */
//                void compute_dsmoddu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >             & adsmoddu );
//
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute STilde = SBar + viscosity * fv2 / ( kappa² * d² )
//                 * @param[ out ] STilde
//                 */
//                real compute_stilde();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of STilde wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adSTildedu a matrix to fill with dSTildedu
//                 */
//                void compute_dstildedu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adstildedu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute ft2 = ct3 * exp( -ct4 * chi²), chi = viscosityDof/viscosityProp
//                 * @param[ out ] ft2
//                 */
//                real compute_ft2();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of ft2 wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adft2du a matrix to fill with dft2du
//                 */
//                void compute_dft2du(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adft2du );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute r = min( viscosity / ( STilde * kappa² * d² ), 10 )
//                 * @param[ out ] r
//                 */
//                real compute_r();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of r wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adrdu a matrix to fill with drdu
//                 */
//                void compute_drdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adrdu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute g = r + cw2 * ( r⁶ - r )
//                 * @param[ out ] g
//                 */
//                real compute_g();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of g wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adgdu a matrix to fill with dgdu
//                 */
//                void compute_dgdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adgdu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute fw = g * ( ( 1 + cw3⁶ ) / ( g⁶ + cw3⁶ ) )^(1/6)
//                 * g = r + cw2 * ( r⁶ - r )
//                 * r = min( viscosity / ( STilde * kappa² * d² ), 10 )
//                 * @param[ out ] fw
//                 */
//                real compute_fw();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fw wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfwdu a matrix to fill with dfwdu
//                 */
//                void compute_dfwdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfwdu );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_SPALART_ALLMARAS_TURBULENCE_HPP_ */
