/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_SP_SUPG_SPALART_ALLMARAS_TURBULENCE_HPP_

#include <map>
// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
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

        class SP_SUPG_Spalart_Allmaras_Turbulence : public Stabilization_Parameter
        {

            //------------------------------------------------------------------------------

          private:
            // default dof type
            MSI::Dof_Type mLeaderDofViscosity = MSI::Dof_Type::VISCOSITY;
            MSI::Dof_Type mLeaderDofVelocity  = MSI::Dof_Type::VX;

            // internal threshold
            const real mEpsilon = MORIS_REAL_EPS;

            // local constitutive enums
            enum class IWG_Constitutive_Type
            {
                SPALART_ALLMARAS_TURBULENCE,
                MAX_ENUM
            };

            // flag for evaluation
            bool                    mLengthScaleEval = true;
            moris::Matrix< DDBMat > mdLengthScaledLeaderDofEval;

            // storage
            real                            mLengthScale;
            moris::Cell< Matrix< DDRMat > > mdLengthScaledLeaderDof;

            /*
             * Rem: mParameters( 0 ) -
             * 0 - exponent for (tauA^r + tauK^r + tauK^r)^(-1/r)
             * 1 - use reaction term 0.0 for no, 1.0 for yes
             */

            // parameters
            bool mHasReaction = true;
            real mExponent    = 2.0;

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_SUPG_Spalart_Allmaras_Turbulence();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_SUPG_Spalart_Allmaras_Turbulence(){};

            //------------------------------------------------------------------------------
            /**
             * set parameters
             * @param[ in ] aParameters a list of parameters
             */
            void set_parameters( moris::Cell< Matrix< DDRMat > > aParameters );

            //------------------------------------------------------------------------------
            /**
             * set function pointers for evaluation
             */
            void set_function_pointers();

            //------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             * child implementation
             */
            void reset_eval_flags();

            //------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aDofStrings list of strings describing the dof types
             * @param[ in ] aIsLeader enum for leader or follower
             */
            void set_dof_type_list(
                    moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                    moris::Cell< std::string >&                  aDofStrings,
                    mtk::Leader_Follower                            aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             * child implementation
             */
            void build_global_dof_type_list();

            //------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes   a cell of group of dv types
             * @param[ in ] aDvStrings list of strings describing the dv types
             * @param[ in ] aIsLeader enum for leader or follower
             */
            void
            set_dv_type_list(
                    moris::Cell< moris::Cell< PDV_Type > >& aDvTypes,
                    moris::Cell< std::string >&             aDvStrings,
                    mtk::Leader_Follower                       aIsLeader = mtk::Leader_Follower::LEADER )
            {
                Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsLeader );
            }

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            void eval_SP();

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a leader dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dSPdLeaderDOF( const moris::Cell< MSI::Dof_Type >& aDofTypes );

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a leader dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            void
            eval_dSPdLeaderDV( const moris::Cell< PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdLeaderDV - not implemented." );
            }

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * return the length scale
             * @param[ out ] mLengthScale length scale for SUPG
             */
            real length_scale();

            /**
             * evaluate the length scale parameter
             */
            void eval_length_scale();

            //------------------------------------------------------------------------------
            /**
             * return the length scale derivative wrt to a leader dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            const Matrix< DDRMat >& dlengthscaledleaderu( const moris::Cell< MSI::Dof_Type >& aDofTypes );

            /**
             * evaluate the length scale derivative wrt to a leader dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            void eval_dlengthscaledleaderu( const moris::Cell< MSI::Dof_Type >& aDofTypes );
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_SPALART_ALLMARAS_TURBULENCE_HPP_ */

