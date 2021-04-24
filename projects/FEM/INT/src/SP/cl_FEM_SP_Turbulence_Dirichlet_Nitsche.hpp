/*
 * cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp
 *
 *  Created on: Jun 04, 2020
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_Turbulence_Dirichlet_Nitsche : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // default tuple for element size to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Master_Slave > mElementSizeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_LENGTH_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Master_Slave::MASTER );

                // populate the dof map (default)
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;

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

                // Property type for the SP
                enum class SP_Property_Type
                {
                    VISCOSITY,
                    MAX_ENUM
                };

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * Rem: mParameters( 0 ) - gamma penalty parameter
                 */
                SP_Turbulence_Dirichlet_Nitsche();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Turbulence_Dirichlet_Nitsche(){};

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
                 * dPPdMasterDOF ( 1 x numDerDof )
                 */
                void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dPPdMasterDV ( 1 x numDerDv )
                 */
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Turbulence_Dirichlet_Nitsche - eval_dSPdMasterDV: not implemented." );
                }

                //------------------------------------------------------------------------------
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
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_TURBULENCE_DIRICHLET_NITSCHE_HPP_ */
