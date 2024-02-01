/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Velocity_SlipBoundary_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_VELOCITY_SLIPBOUNDARY_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_VELOCITY_SLIPBOUNDARY_NITSCHE_HPP_

#include <map>
//MRS/CNT/src
#include "moris_typedefs.hpp"
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
        /*
         * Stabilization parameter for SlipBoundary BC on velocity for fluid problem
         * applied with Nitsche's formulation
         *
         * normal direction:
         *            gamma_n = alpha_N * ( viscosity / h
         *                     + density * norm_inf( u ) / 6.0
         *                     + density * h / (12 alpha_time Delta_t )
         *
         * from Schott et al. (2015)
         *
         * tangential direction:
         *             gamma_t1 = alpha_tgang/(alpha_tang*sliplength+h)
         *             gamma_t2 = h/( alpha_tang*sliplength+h )
         */
        class SP_Velocity_SlipBoundary_Nitsche : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // populate the dof map (default)
                MSI::Dof_Type mLeaderDofVelocity = MSI::Dof_Type::VX;

                // default tuple for element size to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Leader_Follower > mElementSizeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_LENGTH_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::LEADER );

                // property type for the SP
                enum class Property_Type
                {
                        VISCOSITY,  // fluid viscosity
                        DENSITY,    // fluid density
                        SLIPLENGTH, // slip length
                        MAX_ENUM
                };

                /*
                 * Rem: mParameters( 0 )(0) - alpha_n
                 *      mParameters( 1 )(0) - alpha_time
                 *      mParameters( 2 )(0) - alpha_tang
                 */

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Velocity_SlipBoundary_Nitsche();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Velocity_SlipBoundary_Nitsche(){};

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
                        mtk::Leader_Follower                             aIsLeader = mtk::Leader_Follower::LEADER );

                //------------------------------------------------------------------------------
                /**
                 * set dv types
                 * @param[ in ] aDvTypes   a cell of group of dv types
                 * @param[ in ] aDvStrings list of strings describing the dv types
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< gen::PDV_Type > > & aDvTypes,
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
                 * evaluate the stabilization parameter values
                 *
                 * mPPVal(0) : gamma_n  -  normal direction
                 * mPPVal(1) : gamma_t1 -  tangential direction
                 * mPPVal(2) : gamma_t2 -  tangential direction
                 */
                void eval_SP();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the stabilization parameter derivative wrt to a leader dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_dSPdLeaderDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the penalty parameter derivative wrt to a leader dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 */
                void eval_dSPdLeaderDV( const moris::Cell< gen::PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Velocity_SlipBoundary_Nitsche::eval_dSPdLeaderDV - not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_VELOCITY_SLIPBOUNDARY_NITSCHE_HPP_ */

