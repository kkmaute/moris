/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Stab_Penalty_Contact.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"         //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_FEM_Cluster.hpp"                    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class SP_Stab_Penalty_Contact : public Stabilization_Parameter
    {

        //------------------------------------------------------------------------------

      private:
        real mLeaderVolume     = 0.5;    // volume on leader
        real mFollowerVolume   = 0.5;    // volume on follower
        real mInterfaceSurface = 1.0;    // surface on leader/follower interface

        enum class SP_Property_Type
        {
            MATERIAL,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        SP_Stab_Penalty_Contact();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~SP_Stab_Penalty_Contact() override {};

        //------------------------------------------------------------------------------
        /**
         * reset the cluster measures required for this SP
         */
        void reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mLeaderVolume     = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::INTERP, mtk::Leader_Follower::LEADER );
            mFollowerVolume   = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::INTERP, mtk::Leader_Follower::FOLLOWER );
            mInterfaceSurface = mCluster->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER );
            //                std::cout<<"mInterfaceSurface "<<mInterfaceSurface<<std::endl;
            //                std::cout<<"mInterfaceSurface in FEM "<<mCluster->compute_volume()<<std::endl;
        }

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
                Vector< gen::PDV_Type > &aDvTypes,
                Vector< std::string >   &aDvStrings,
                mtk::Leader_Follower     aIsLeader = mtk::Leader_Follower::LEADER ) override
        {
            Stabilization_Parameter::set_dv_type_list( aDvTypes, aIsLeader );
        }

        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter value
         */
        void eval_SP() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the stabilization parameter derivative wrt to a leader dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dSPdLeaderDOF ( 1 x numDerDof )
         */
        void eval_dSPdLeaderDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override
        {
            MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdLeaderDOF: not implemented." );
        }
        //------------------------------------------------------------------------------
        /**
         * evaluate the stabilization parameter derivative wrt to a follower dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dSPdFollowerDOF ( 1 x numDerDof )
         */
        void eval_dSPdFollowerDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override
        {
            MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdFollowerDOF: not implemented." );
        }
        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter derivative wrt to a leader dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         * dPPdLeaderDV ( 1 x numDerDv )
         */
        void eval_dSPdLeaderDV( const Vector< gen::PDV_Type > &aDvTypes ) override
        {
            MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdLeaderDV: not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * evaluate the penalty parameter derivative wrt to a follower dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         * dSPdFollowerDV ( 1 x numDerDv )
         */
        void eval_dSPdFollowerDV( const Vector< gen::PDV_Type > &aDvTypes ) override
        {
            MORIS_ERROR( false, "SP_Stab_Penalty_Contact::eval_dSPdFollowerDV: not implemented." );
        }
        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_SP_STAB_PENALTY_CONTACT_HPP_ */
