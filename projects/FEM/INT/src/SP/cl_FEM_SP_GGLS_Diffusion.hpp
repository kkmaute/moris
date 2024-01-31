/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_GGLS_Diffusion.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_GGLS_DIFFUSION_HPP_
#define SRC_FEM_CL_FEM_SP_GGLS_DIFFUSION_HPP_

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
         * Stabilization parameter for SUPG stabilization on diffusion-advection
         * tau_T =
         * from Tezduyar & Osawa (2000)
         */
        class SP_GGLS_Diffusion : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                // default tuple for element size to define cluster measure
                std::tuple<fem::Measure_Type,mtk::Primary_Void,mtk::Leader_Follower > mElementSizeTuple =
                        std::make_tuple(
                                fem::Measure_Type::CELL_LENGTH_MEASURE,
                                mtk::Primary_Void::PRIMARY,
                                mtk::Leader_Follower::LEADER );

                // default dof type
                MSI::Dof_Type mLeaderDofTemp = MSI::Dof_Type::TEMP;

                // internal threshold
                const real mEpsilon = 1e-18;

                // property type for the SP
                enum class Property_Type
                {
                    CONDUCTIVITY,
                    DENSITY,
                    HEAT_CAPACITY,
                    LATENT_HEAT,
                    PC_TEMP,
                    PHASE_STATE_FUNCTION,
                    PHASE_CHANGE_CONST,
                    MAX_ENUM
                };

            public:

                /*
                 * Rem: mParameters - no parameters needed
                 */

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_GGLS_Diffusion();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_GGLS_Diffusion(){};

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
                        moris::Cell< moris::Cell< ge::PDV_Type > > & aDvTypes,
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
                 * evaluate the stabilization parameter value
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
                void eval_dSPdLeaderDV( const moris::Cell< ge::PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_SUPG_Advection::eval_dSPdLeaderDV - not implemented." );
                }
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_GGLS_DIFFUSION_HPP_ */

