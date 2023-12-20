/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Lagrange_Multiplier_L2.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_LAGRANGE_MULTIPLIER_L2_HPP_
#define SRC_FEM_CL_FEM_SP_LAGRANGE_MULTIPLIER_L2_HPP_

//#ifndef SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_
//#define SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"
namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class SP_Lagrange_Multiplier_L2 : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // default tuple for element size to define cluster measure
                std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mElementSizeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_LENGTH_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::LEADER );

                // Property type for the SP
                enum class SP_Property_Type
                {
                    MATERIAL,
                    MAX_ENUM
                };

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * Rem: mParameters( 0 ) - gamma penalty parameter
                 */
                SP_Lagrange_Multiplier_L2();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Lagrange_Multiplier_L2(){};

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< std::string >                  & aDofStrings,
                        mtk::Leader_Follower                             aIsLeader = mtk::Leader_Follower::LEADER )
                {
                    Stabilization_Parameter::set_dof_type_list( aDofTypes, aIsLeader );
                }

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                        moris::Cell< std::string >             & aDvStrings,
                        mtk::Leader_Follower                        aIsLeader = mtk::Leader_Follower::LEADER )
                {
                    Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsLeader );
                }

                //------------------------------------------------------------------------------
                /**
                 * get cluster measure tuples
                 * @param[ in ] aClusterMeasureTuples list of tuples describing the cluster measure types
                 */
                moris::Cell< std::tuple<
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
                 * evaluate the penalty parameter derivative wrt to a leader dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dPPdLeaderDOF ( 1 x numDerDof )
                 */
                void eval_dSPdLeaderDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a leader dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 * dPPdLeaderDV ( 1 x numDerDv )
                 */
                void eval_dSPdLeaderDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Lagrange_Multiplier_L2 - eval_dSPdLeaderDV: not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_DIRICHLET_NITSCHE_HPP_ */

