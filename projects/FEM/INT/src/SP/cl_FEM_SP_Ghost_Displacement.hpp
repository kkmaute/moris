/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Ghost_Displacement.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_
#define SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"     //FEM/INT/src
#include "cl_FEM_Cluster.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class SP_Ghost_Displacement : public Stabilization_Parameter
    {

        //------------------------------------------------------------------------------

      private:
        // default tuple for element size to define cluster measure
        std::tuple< fem::Measure_Type, mtk::Primary_Void, mtk::Leader_Follower > mElementSizeTuple =
                std::make_tuple(
                        fem::Measure_Type::CELL_LENGTH_MEASURE,
                        mtk::Primary_Void::PRIMARY,
                        mtk::Leader_Follower::LEADER );

        // order of weak form of governing questions
        real mWeakFormOrder;

        // property type for the SP
        enum class SP_Property_Type
        {
            MATERIAL,
            MAX_ENUM
        };

        // function to determine and set the weak form order the first time
        real get_weak_form_order_init();

        // function to get the weak form order after being initialized
        real get_weak_form_order();

        // init function pointer to get the weak form order. Pointer used to avoid a conditional call for every iteration
        real ( SP_Ghost_Displacement::*mGetWeakFormOrder )() = &SP_Ghost_Displacement::get_weak_form_order_init;

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * Rem: mParameters( 0 ) - gamma penalty parameter
         */
        SP_Ghost_Displacement();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~SP_Ghost_Displacement() override{};

        //------------------------------------------------------------------------------
        /**
         * set dof types
         * @param[ in ] aDofTypes a cell of cell of dof types
         * @param[ in ] aDofStrings list of strings describing the dof types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > > &aDofTypes,
                Vector< std::string >             &aDofStrings,
                mtk::Leader_Follower               aIsLeader = mtk::Leader_Follower::LEADER ) override
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
                Vector< Vector< gen::PDV_Type > > &aDvTypes,
                Vector< std::string >             &aDvStrings,
                mtk::Leader_Follower               aIsLeader = mtk::Leader_Follower::LEADER ) override
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
         * evaluate the penalty parameter value
         */
        void eval_SP() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter derivative wrt to a leader dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dPPdLeaderDOF ( 1 x numDerDof )
         */
        void eval_dSPdLeaderDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter derivative wrt to a leader dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         * dPPdLeaderDV ( 1 x numDerDv )
         */
        void eval_dSPdLeaderDV( const Vector< gen::PDV_Type > &aDvTypes ) override
        {
            MORIS_ERROR( false, "SP_Ghost_Displacement::eval_dSPdLeaderDV: not implemented." );
        }

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_SP_GHOST_DISPLACEMENT_HPP_ */

