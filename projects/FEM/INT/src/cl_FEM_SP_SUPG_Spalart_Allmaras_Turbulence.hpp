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

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_SUPG_Spalart_Allmaras_Turbulence : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private :

                // cluster measures
                real mElementSize = 1.0;

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

            public :

                // property type for the SP
                enum class Property_Type
                {
                    DENSITY,   // fluid density
                    VISCOSITY, // fluid viscosity
                    WALL_DISTANCE,
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, Property_Type > mPropertyMap;

                /*
                 * Rem: mParameters( 0 ) -
                 */

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
                 * reset the cluster measures required for this SP
                 */
                void reset_cluster_measures();

                //------------------------------------------------------------------------------
                /**
                 * set function pointers for evaluation
                 */
                void set_function_pointers();

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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
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
                        moris::Cell< MSI::Dof_Type >   aDofTypes,
                        Matrix< DDRMat >          & adfwdu );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_SPALART_ALLMARAS_TURBULENCE_HPP_ */
