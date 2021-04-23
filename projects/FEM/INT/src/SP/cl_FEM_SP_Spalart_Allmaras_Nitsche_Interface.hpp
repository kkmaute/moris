/*
 * cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp
 *
 *  Created on: Aug 12, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_

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

        class SP_Spalart_Allmaras_Nitsche_Interface : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // populate the dof map (default)
                MSI::Dof_Type mMasterDofViscosity = MSI::Dof_Type::VISCOSITY;
                MSI::Dof_Type mSlaveDofViscosity  = MSI::Dof_Type::VISCOSITY;

                // default tuple for master volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Master_Slave > mMasterVolumeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Master_Slave::MASTER );

                // default tuple for slave volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Master_Slave > mSlaveVolumeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Master_Slave::SLAVE );

                // default tuple for slave volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Master_Slave > mInterfaceSurfaceTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_SIDE_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Master_Slave::MASTER );

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

                enum class SP_Property_Type
                {
                    MATERIAL,
                    MAX_ENUM
                };

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Spalart_Allmaras_Nitsche_Interface();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Spalart_Allmaras_Nitsche_Interface(){};

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
                 * evaluate the stabilization parameter derivative wrt to a master dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdMasterDOF ( 1 x numDerDof )
                 */
                void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a slave dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdSlaveDOF ( 1 x numDerDof )
                 */
                void eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a master dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dPPdMasterDV ( 1 x numDerDv )
                 */
                void eval_dSPdMasterDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Nitsche_Interface::eval_dSPdMasterDV: not implemented." );
                }

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a slave dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dSPdSlaveDV ( 1 x numDerDv )
                 */
                void eval_dSPdSlaveDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Nitsche_Interface::eval_dSPdSlaveDV: not implemented." );
                }
                //------------------------------------------------------------------------------

            private:
//                //------------------------------------------------------------------------------
//
//                real compute_diffusion_coefficient(
//                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );
//
//                //------------------------------------------------------------------------------
//
//                void compute_ddiffusiondu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & addiffusiondu,
//                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );
//
//                //------------------------------------------------------------------------------
//
//                real compute_fn(
//                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );
//
//                //------------------------------------------------------------------------------
//
//                void compute_dfndu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfndu,
//                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );
//
//                //------------------------------------------------------------------------------
//
//                real compute_chi(
//                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );
//
//                //------------------------------------------------------------------------------
//
//                void compute_dchidu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adchidu,
//                        mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER );

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_ */
