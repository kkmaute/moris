/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Incompressible_Flow.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_
#define SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_FEM_Cluster.hpp"

namespace moris
{
    namespace fem
    {
    /*
     * Stabilization parameter for incompressible NS
     * from Schott et al. (2015) DOI: 10.1002/nme.4789
     * see also: Whiting (1999) PhD thesis
     *
     * Note: contribution from strong form of residual is applied in IWG
     */
        //------------------------------------------------------------------------------

        class SP_Incompressible_Flow : public Stabilization_Parameter
        {

          private:
            // default dof type
            MSI::Dof_Type mLeaderDofVelocity = MSI::Dof_Type::VX;
            MSI::Dof_Type mLeaderDofPressure = MSI::Dof_Type::P;

            // internal threshold
            const real mEpsilon = MORIS_REAL_EPS;

            // property type for the SP
            enum class Property_Type
            {
                DENSITY,             // fluid density
                VISCOSITY,           // fluid viscosity
                INV_PERMEABILITY,    // inverse of the permeability for flow through porous media
                MAX_ENUM
            };

            /*
             * Rem: mParameters( 0 ) - CI = 36,60,128 for linear, quadratic, cubic;
             * see: STABILIZED FINITE ELEMENT METHODS FOR FLUID DYNAMICS USING A HIERARCHICAL BASIS
             *      by Christian H. Whiting
             */

            // element inverse estimate parameter
            real mCI = 0.0;

            // is time solve parameter
            real mBetaTime    = 1.0;
            bool mSetBetaTime = false;

          public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            SP_Incompressible_Flow();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~SP_Incompressible_Flow(){};

            //------------------------------------------------------------------------------
            /**
             * set parameters
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
             * evaluate the penalty parameter value
             */
            void eval_SP();

            //------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt to a leader dof type
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
                MORIS_ERROR( false, "SP_Incompressible_Flow::eval_dSPdLeaderDV - not implemented." );
            }

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SP_INCOMPRESSIBLE_FLOW_HPP_ */

