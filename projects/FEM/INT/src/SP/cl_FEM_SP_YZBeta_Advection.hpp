/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_YZBeta_Advection.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_YZBeta_ADVECTION_HPP_
#define SRC_FEM_CL_FEM_SP_YZBeta_ADVECTION_HPP_

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
         * Stabilization parameter for YZBeta discontinuity method
         * from from Bazilevs et al. (2007) DOI: 10.1002/fld.1484
         * see also: Tezduyar, Ramakrishnan and Sathe (2008) DOI: 10.1002/fld.174
         *
         * Note: contribution from strong form of residual is applied in IWG
         */
        class SP_YZBeta_Advection : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:
                // default dof type
                MSI::Dof_Type mLeaderDofScalarField = MSI::Dof_Type::TEMP;

                // property type for the SP
                enum class Property_Type
                {
                    BETA_CONSTANT,
                    REFERENCE_STATE,
                    MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                    DIFFUSION,
                    MAX_ENUM
                };

                // internal threshold
                const real mEpsilon = 1e-18;

                /*
                 * Rem: mParameters - no parameters needed
                 */

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_YZBeta_Advection();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_YZBeta_Advection(){}

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
                    MORIS_ERROR( false, "SP_YZBeta_Advection::eval_dSPdLeaderDV - not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_VELOCITY_DIRICHLET_NITSCHE_HPP_ */

