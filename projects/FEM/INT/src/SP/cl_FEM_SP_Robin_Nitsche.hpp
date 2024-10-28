/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Robin_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_ROBIN_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_SP_ROBIN_NITSCHE_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_FEM_Cluster.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------
    /*
     * Stabilization parameter for SlipBoundary BC on velocity for fluid problem
     * applied with Nitsche's formulation
     *
     * by Juntunen abd Stenberg 2009
     *
     * tangential direction:
     *             alpha_tang = 1/gamma ;
     *             gamma_t1 = alpha_tgang/(alpha_tang*sliplength+h)
     *             gamma_t2 = h/( alpha_tang*sliplength+h )
     */
    class SP_Robin_Nitsche : public Stabilization_Parameter
    {

        //------------------------------------------------------------------------------

      private:
        // populate the dof map (default)
        MSI::Dof_Type mLeaderDofTemp = MSI::Dof_Type::THETA;

        // default tuple for element size to define cluster measure
        std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mElementSizeTuple =
                std::make_tuple(
                        fem::Measure_Type::CELL_LENGTH_MEASURE,
                        mtk::Primary_Void::PRIMARY,
                        mtk::Leader_Follower::LEADER );

        // property type for the SP
        enum class Property_Type
        {
            NEUMANN_PENALTY,
            MAX_ENUM
        };

        /*
         *      mParameters( 0 )(0) - alpha_tang
         */

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        SP_Robin_Nitsche();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~SP_Robin_Nitsche() override {};

        //------------------------------------------------------------------------------
        /**
         * set dof types
         * @param[ in ] aDofTypes a cell of cell of dof types
         * @param[ in ] aDofStrings list of strings describing the dof types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                Vector< std::string >&             aDofStrings,
                mtk::Leader_Follower               aIsLeader = mtk::Leader_Follower::LEADER ) override;

        //------------------------------------------------------------------------------
        /**
         * set dv types
         * @param[ in ] aDvTypes   a cell of group of dv types
         * @param[ in ] aDvStrings list of strings describing the dv types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void
        set_dv_type_list(
                Vector< gen::PDV_Type >& aDvTypes,
                Vector< std::string >&   aDvStrings,
                mtk::Leader_Follower     aIsLeader = mtk::Leader_Follower::LEADER ) override
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
        get_cluster_measure_tuple_list() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the stabilization parameter values
         *
         * mPPVal(0) : gamma_t1 -  stability coeff
         * mPPVal(1) : gamma_t2 -  adjoint coeff
         */
        void eval_SP() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the stabilization parameter derivative wrt to a leader dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dSPdLeaderDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter derivative wrt to a leader dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         */
        void
        eval_dSPdLeaderDV( const Vector< gen::PDV_Type >& aDvTypes ) override
        {
            MORIS_ERROR( false, "SP_Robin_Nitsche::eval_dSPdLeaderDV - not implemented." );
        }

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_SP_Robin_Nitsche_HPP_ */
