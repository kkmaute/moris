/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_

#include <map>
//MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
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

        class SP_Spalart_Allmaras_Nitsche_Interface : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // populate the dof map (default)
                MSI::Dof_Type mLeaderDofViscosity = MSI::Dof_Type::VISCOSITY;
                MSI::Dof_Type mFollowerDofViscosity  = MSI::Dof_Type::VISCOSITY;

                // default tuple for leader volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Leader_Follower > mLeaderVolumeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::LEADER );

                // default tuple for follower volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Leader_Follower > mFollowerVolumeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::FOLLOWER );

                // default tuple for follower volume to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Leader_Follower > mInterfaceSurfaceTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_SIDE_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::LEADER );

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        SPALART_ALLMARAS_TURBULENCE,
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
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void set_dof_type_list(
                        Vector< Vector< MSI::Dof_Type > > & aDofTypes,
                        Vector< std::string >                  & aDofStrings,
                        mtk::Leader_Follower                             aIsLeader = mtk::Leader_Follower::LEADER );

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void set_dv_type_list(
                        Vector< Vector< gen::PDV_Type > > & aDvTypes,
                        Vector< std::string >             & aDvStrings,
                        mtk::Leader_Follower                        aIsLeader = mtk::Leader_Follower::LEADER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsLeader );
                }

                //------------------------------------------------------------------------------
                /**
                 * get cluster measure tuples
                 * @param[ in ] aClusterMeasureTuples list of tuples describing the cluster measure types
                 */
                Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > > get_cluster_measure_tuple_list();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter value
                 */
                void eval_SP();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a leader dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdLeaderDOF ( 1 x numDerDof )
                 */
                void eval_dSPdLeaderDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a follower dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dSPdFollowerDOF ( 1 x numDerDof )
                 */
                void eval_dSPdFollowerDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a leader dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dPPdLeaderDV ( 1 x numDerDv )
                 */
                void eval_dSPdLeaderDV( const Vector< gen::PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Nitsche_Interface::eval_dSPdLeaderDV: not implemented." );
                }

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a follower dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dSPdFollowerDV ( 1 x numDerDv )
                 */
                void eval_dSPdFollowerDV( const Vector< gen::PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Nitsche_Interface::eval_dSPdFollowerDV: not implemented." );
                }
                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_SP_SPALART_ALLMARAS_NITSCHE_INTERFACE_HPP_ */

