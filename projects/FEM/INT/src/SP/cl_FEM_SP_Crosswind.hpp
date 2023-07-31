/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Crosswind.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_CROSSWIND_HPP_
#define SRC_FEM_CL_FEM_SP_CROSSWIND_HPP_

#include <map>
//MRS/CNT/src
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

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /*
         * Stabilization parameter for crosswind stabilization
         * from Shadid et al. (2016) DOI:
         * see also: Codina (1992) DOI:
         *
         *
         * Note: contribution from strong form of residual is applied in IWG
         */
        class SP_Crosswind : public Stabilization_Parameter
        {

                //------------------------------------------------------------------------------
            private:

                /*
                 * Rem: mParameters( 0 ) - C = 0.70, for linear,
                 *                             0.35, for quadratic.
                 * see: Codina (1992) DOI:
                 *
                 *      mParameters( 1 ) - epsilon = 1.0e-10, default
                 *      check against zero
                 */

                // constant for crosswind stabilization
                real mC = 0.0;

                // internal threshold for zero
                real mEpsilon = 1.0e-10;

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                SP_Crosswind();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~SP_Crosswind(){}

                //------------------------------------------------------------------------------
                /**
                 * set dof types
                 * @param[ in ] aDofTypes a cell of cell of dof types
                 * @param[ in ] aDofStrings list of strings describing the dof types
                 * @param[ in ] aIsLeader enum for leader or follower
                 */
                void
                set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                        moris::Cell< std::string >&                  aDofStrings,
                        mtk::Leader_Follower                            aIsLeader = mtk::Leader_Follower::LEADER )
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
                 * set parameters
                 */
                void set_parameters( moris::Cell< Matrix< DDRMat > > aParameters );

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
                void eval_dSPdLeaderDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, "SP_Crosswind::eval_dSPdLeaderDV - not implemented." );
                }

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_CROSSWIND_HPP_ */

