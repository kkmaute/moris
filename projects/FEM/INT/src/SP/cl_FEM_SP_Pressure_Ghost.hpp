/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Pressure_Ghost.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_
#define SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_

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

namespace moris::fem
{
    //------------------------------------------------------------------------------
    /*
     * Stabilization parameter for face-oriented pressure ghost-penalty
     * for different flow regimes,
     * i.e., viscous, convective, and transient regime and contribution from inverse permeability
     * gamma_p = alpha_p * delta_p * h^(2*i)
     * where i       = interpolation order
     *       delta_p = ( ( viscosity / h )
     *                 + ( density * norm_inf( u ) / 6.0 )
     *                 + ( density * h / ( 12.0 * theta * deltat ) )
     *                 + ( impermeability * h / 12.0 ) )^(-1)
     * from Schott et al. (2015). Note this reference does not account for impermeability
     *
     */
    class SP_Pressure_Ghost : public Stabilization_Parameter
    {

        //------------------------------------------------------------------------------

      private:
        // default tuple for element size to define cluster measure
        std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mElementSizeTuple =
                std::make_tuple(
                        fem::Measure_Type::CELL_LENGTH_MEASURE,
                        mtk::Primary_Void::PRIMARY,
                        mtk::Leader_Follower::LEADER );

        // default dof type
        MSI::Dof_Type mLeaderDofVelocity = MSI::Dof_Type::VX;

        // property type for the SP
        enum class Property_Type
        {
            VISCOSITY,           // fluid viscosity
            DENSITY,             // fluid density
            INV_PERMEABILITY,    // fluid inverse permeability
            MAX_ENUM
        };

        // parameters
        real mAlphaP = 0.0;
        real mTheta  = 0.0;

        // set parameter bool
        bool mSetAlphaP = false;
        bool mSetTheta  = false;

      public:
        /*
         * Rem: mParameters( 0 ) - alpha_p
         *      mParameters( 1 ) - theta
         */

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        SP_Pressure_Ghost();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~SP_Pressure_Ghost() override {};

        //------------------------------------------------------------------------------

        void set_parameters( const Vector< Matrix< DDRMat > >& aParameters ) override;

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
         * evaluate the stabilization parameter value
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
            MORIS_ERROR( false, "SP_Pressure_Ghost::eval_dSPdLeaderDV - not implemented." );
        }

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_SP_PRESSURE_GHOST_HPP_ */
