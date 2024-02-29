/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_VELOCITY_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_VELOCITY_DIRICHLET_NITSCHE_HPP_

#include <map>
// MRS/CNT/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// FEM/INT/src
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
         * Stabilization parameter for Dirichlet BC on velocity for fluid problem
         * applied with Nitsche's formulation
         *
         * gamma_N = alpha_N * ( viscosity / h
         *                     + density * norm_inf( u ) / 6.0
         *                     + density * h / (12 alpha_time Delta_time )
         *
         * from Schott et al. (2015) and contribution from impermeability
         *
         * gamma_N = gamma_N + alpha_N * beta * h / 12
         */
        class SP_Velocity_Dirichlet_Nitsche : public Stabilization_Parameter
        {

            //------------------------------------------------------------------------------

          private:
            // populate the dof map (default)
            MSI::Dof_Type mLeaderDofVelocity = MSI::Dof_Type::VX;

            // default tuple for element size to define cluster measure
            std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mElementSizeTuple =
                    std::make_tuple(
                            fem::Measure_Type::CELL_LENGTH_MEASURE,
                            mtk::Primary_Void::PRIMARY,
                            mtk::Leader_Follower::LEADER );

            // default tuple for leader volume to define cluster measure
            std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mLeaderVolumeTuple =
                    std::make_tuple(
                            fem::Measure_Type::CELL_MEASURE,
                            mtk::Primary_Void::PRIMARY,
                            mtk::Leader_Follower::LEADER );

            std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mInterfaceSurfaceTuple =
                    std::make_tuple(
                            fem::Measure_Type::CELL_SIDE_MEASURE,
                            mtk::Primary_Void::PRIMARY,
                            mtk::Leader_Follower::LEADER );

            // property type for the SP
            enum class Property_Type
            {
                VISCOSITY,           // fluid viscosity
                DENSITY,             // fluid density
                INV_PERMEABILITY,    // fluid impermeability (Brinkman term)
                MAX_ENUM
            };

            /*
             * Rem: mParameters( 0 )      - alpha_N
             *      mParameters( 1 )( 0 ) - alpha_time
             *      mParameters( 2 )( 0 ) - geometry formulation for cluster measure
             *                              0: cluster edge length
             *                              1: surface measure / primary volume
             */

            // parameters
            real mAlphaN    = 0.0;
            real mAlphaTime = 0.0;
            uint mGeometryFormulation = 0;

            // set parameter bool
            bool mSetAlphaN    = false;
            bool mSetAlphaTime = false;

            // pointer to function for geometry measure
            real ( SP_Velocity_Dirichlet_Nitsche::*mGeometryMeasureFunc )() = nullptr;

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Velocity_Dirichlet_Nitsche();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Velocity_Dirichlet_Nitsche(){};

            //------------------------------------------------------------------------------
            /**
             * set parameters
             * @param[ in ] aParameters a list of parameters
             */
            void set_parameters( Vector< Matrix< DDRMat > > aParameters );

            //------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aDofStrings list of strings describing the dof types
             * @param[ in ] aIsLeader enum for leader or follower
             */
            void set_dof_type_list(
                    Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                    Vector< std::string >&                  aDofStrings,
                    mtk::Leader_Follower                            aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes   a cell of group of dv types
             * @param[ in ] aDvStrings list of strings describing the dv types
             * @param[ in ] aIsLeader enum for leader or follower
             */
            void
            set_dv_type_list(
                    Vector< Vector< gen::PDV_Type > >& aDvTypes,
                    Vector< std::string >&             aDvStrings,
                    mtk::Leader_Follower                       aIsLeader = mtk::Leader_Follower::LEADER )
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
                    mtk::Leader_Follower > >
            get_cluster_measure_tuple_list();

            //------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter value
             */
            void eval_SP();

            //------------------------------------------------------------------------------
            /**
             * evaluate the stabilization parameter derivative wrt to a leader dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dSPdLeaderDOF( const Vector< MSI::Dof_Type >& aDofTypes );

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a leader dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            void
            eval_dSPdLeaderDV( const Vector< gen::PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, "SP_Velocity_Dirichlet_Nitsche::eval_dSPdLeaderDV - not implemented." );
            }

          private:

            //------------------------------------------------------------------------------
            /**
             * evaluate the specified geometry measure for the stabilization parameter
             * @param[ out ] real value of specified geometry measure
             */
            real eval_geometry_measure_edge_length();
            real eval_geometry_measure_relative_edge_length();

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_VELOCITY_DIRICHLET_NITSCHE_HPP_ */
