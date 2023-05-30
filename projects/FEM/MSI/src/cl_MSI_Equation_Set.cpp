/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Set.cpp
 *
 */

#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace MSI
    {
        //------------------------------------------------------------------------------

        moris::Cell< moris::Cell< MSI::Dof_Type > > & Equation_Set::get_dof_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderDofTypes;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerDofTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_type_list - can only be LEADER or FOLLOWER");
                    return mLeaderDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        Matrix< DDSMat > & Equation_Set::get_dof_type_map(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderDofTypeMap;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerDofTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dof_type_map - can only be LEADER or FOLLOWER");
                    return mLeaderDofTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dof_index_for_type(
                enum MSI::Dof_Type aDofType,
                mtk::Leader_Follower  aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // check if leader dof type in list
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mLeaderDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type - leader dof type does not exist in map." );

                    // get the index for leader dof type from map
                    sint tLeaderIndex = mLeaderDofTypeMap( static_cast< int >( aDofType ) );

                    //                    // check if leader dof type assigned in map
                    //                    MORIS_ASSERT( tLeaderIndex != -1,
                    //                                  "Equation_Set::get_dof_index_for_type - leader dof type not assigned in map." );

                    // return leader dof type index
                    return tLeaderIndex;
                }

                case mtk::Leader_Follower::FOLLOWER:
                {
                    // check if follower dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mFollowerDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type - follower dof type does not exist in map." );

                    // get the index for follower dof type from map
                    sint tFollowerIndex = mFollowerDofTypeMap( static_cast< int >( aDofType ) );

                    //                    // check if follower dof type assigned in map
                    //                    MORIS_ASSERT( tFollowerIndex != -1,
                    //                                  "Equation_Set::get_dof_index_for_type - follower dof type not assigned in map." );

                    if ( tFollowerIndex == -1 )
                    {
                        return tFollowerIndex;
                    }
                    else
                    {
                        // get maximum index from leader map
                        moris::sint tMaxLeaderIndex = mLeaderDofTypeMap.max();

                        // check if leader map was assigned
                        MORIS_ASSERT( tMaxLeaderIndex != -1,
                                "Equation_Set::get_dof_index_for_type - mLeaderDofTypeMap is empty." );

                        // return follower dof type index
                        return tFollowerIndex + tMaxLeaderIndex + 1;
                    }
                }

                default:
                {
                    MORIS_ERROR( false,
                            "Equation_Set::get_dof_index_for_type - can only be LEADER or FOLLOWER.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dof_index_for_type_1(
                enum MSI::Dof_Type     aDofType,
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // check if leader dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mLeaderDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type_1 - leader dof type does not exist in map." );

                    //                    // check if leader dof type assigned in map
                    //                    MORIS_ASSERT( mLeaderDofTypeMap( static_cast< int >( aDofType ) ) != -1,
                    //                                  "Equation_Set::get_dof_index_for_type_1 - leader dof type does not assigned in map." );

                    // return index for leader dof type
                    return mLeaderDofTypeMap( static_cast< int >( aDofType ) );
                }

                case mtk::Leader_Follower::FOLLOWER :
                {
                    // check if follower dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mFollowerDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type_1 - follower dof type does not exist in map." );

                    //                    // check if follower dof type assigned in map
                    //                    MORIS_ASSERT( mFollowerDofTypeMap( static_cast< int >( aDofType ) ) != -1,
                    //                                  "Equation_Set::get_dof_index_for_type_1 - follower dof type not assigned in map." );

                    // return index for follower dof type
                    return mFollowerDofTypeMap( static_cast< int >( aDofType ) );
                }

                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_index_for_type_1 - can only be LEADER or FOLLOWER." );
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > > & Equation_Set::get_dv_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderDvTypes;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerDvTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dv_type_list - can only be LEADER or FOLLOWER.");
                    return mLeaderDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Equation_Set::get_dv_type_map(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderDvTypeMap;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerDvTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dv_type_map - can only be LEADER or FOLLOWER");
                    return mLeaderDvTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< mtk::Field_Type > > & Equation_Set::get_field_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderFieldTypes;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerFieldTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_field_type_list - can only be LEADER or FOLLOWER.");
                    return mLeaderFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Equation_Set::get_field_type_map(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderFieldTypeMap;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerFieldTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_field_type_map - can only be LEADER or FOLLOWER");
                    return mLeaderFieldTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dv_index_for_type(
                enum PDV_Type     aDvType,
                mtk::Leader_Follower aIsLeader  )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mLeaderDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type - leader dv type does not exist in map." );

                    //                    // check if dv type assigned in map
                    //                    MORIS_ASSERT( mLeaderDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type - leader dv type not assigned in map." );

                    // return set index for dv type
                    return mLeaderDvTypeMap( static_cast< int >( aDvType ) );
                }

                case mtk::Leader_Follower::FOLLOWER:
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mFollowerDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type(), follower dv type does not exist in map." );

                    //                    // check if dv type assigned in map
                    //                    MORIS_ASSERT( mFollowerDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type - follower dv type not assigned in map." );

                    // get the set index for dv type
                    sint tFollowerIndex = mFollowerDvTypeMap( static_cast< int >( aDvType ) );

                    // if index is -1
                    if ( tFollowerIndex == -1 )
                    {
                        return tFollowerIndex;
                    }
                    else
                    {
                        // get the max set index for dv types
                        moris::sint tMaxLeaderIndex = mLeaderDvTypeMap.max();

                        // check if mLeaderDvTypeMap is set
                        MORIS_ASSERT( tMaxLeaderIndex != -1,
                                "Equation_Set::get_dv_index_for_type - mLeaderDvTypeMap is empty." );

                        // return set index for dv type
                        return tFollowerIndex + tMaxLeaderIndex + 1;
                    }
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_dv_index_for_type - can only be LEADER or FOLLOWER.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dv_index_for_type_1(
                enum PDV_Type            aDvType,
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mLeaderDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type_1 - leader dv type does not exist in map." );

                    //                    // check if dv type is set in map
                    //                    MORIS_ASSERT( mLeaderDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type_1 - leader dv type not assigned in map." );

                    // return set index for dv type
                    return mLeaderDvTypeMap( static_cast< int >( aDvType ) );
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mFollowerDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type_1 - follower dv type does not exist in map." );

                    //                    // check if dv type is set in map
                    //                    MORIS_ASSERT( mFollowerDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type_1 - follower dv type not assigned in map." );

                    // return set index for dv type
                    return mFollowerDvTypeMap( static_cast< int >( aDvType ) );
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_dv_index_for_type_1 - can only be LEADER or FOLLOWER.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_field_index_for_type_1(
                enum mtk::Field_Type aFieldType,
                mtk::Leader_Follower    aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aFieldType ) < mLeaderFieldTypeMap.numel(),
                            "Equation_Set::get_field_index_for_type_1 - leader field type does not exist in map." );

                    //// check if field type is set in map
                    //MORIS_ASSERT( mLeaderDvTypeMap( static_cast< int >( aFieldType ) ),
                    //          "Equation_Set::get_field_index_for_type_1 - leader field type not assigned in map." );

                    // return set index for dv type
                    return mLeaderFieldTypeMap( static_cast< int >( aFieldType ) );
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // check if field type is set in map
                    MORIS_ASSERT( static_cast< uint >( aFieldType ) < mFollowerFieldTypeMap.numel(),
                            "Equation_Set::get_field_index_for_type_1 - follower dv type does not exist in map." );

                    //// check if field type is set in map
                    //MORIS_ASSERT( mFollowerFieldTypeMap( static_cast< int >( aDvType ) ),
                    //        "Equation_Set::get_field_index_for_type_1 - follower field type not assigned in map." );

                    // return set index for dv type
                    return mFollowerFieldTypeMap( static_cast< int >( aFieldType ) );
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_field_index_for_type_1 - can only be LEADER or FOLLOWER.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Equation_Set::free_matrix_memory()
        {
            // if the Jacobian matrix was created
            if ( mJacobianExist )
            {
                // resize it to 0x0
                mJacobian.set_size( 0, 0 );

                // reset the exist flag
                mJacobianExist = false;
            }

            // if the residual cell of matrices was created
            if ( mResidualExist )
            {
                // resize each matrix to 0x0
                for( auto & tResidual : mResidual )
                {
                    tResidual.set_size( 0, 0 );
                }
                mResidual.clear();

                // reset the exist flag
                mResidualExist = false;
            }

            if ( mQIExist )
            {
//                // resize each matrix to 0x0
//                for( auto & tQI : mQI )
//                {
//                    tQI.resize( 0, 0 );
//                }
//                mQI.clear();

                // reset the exist flag
                mQIExist = false;
            }

            mIsStaggered = false;

            // free additional memory
            this->free_memory();
        }

        //------------------------------------------------------------------------------

        const moris::Cell< enum MSI::Dof_Type > & Equation_Set::get_requested_dof_types()
        {
            MORIS_ERROR( mModelSolverInterface != nullptr,
                    "Equation_Set::get_requested_dof_types - model solver interface not set yet." );

            return mModelSolverInterface->get_solver_interface()->get_requested_dof_types();
        }

        //------------------------------------------------------------------------------
        const moris::Cell< enum MSI::Dof_Type > & Equation_Set::get_secondary_dof_types()
        {
            MORIS_ERROR( mModelSolverInterface != nullptr,
                    "Equation_Set::get_requested_dof_types - model solver interface not set yet." );

            return mModelSolverInterface->get_solver_interface()->get_secondary_dof_types();
        }

        //------------------------------------------------------------------------------

        void Equation_Set::create_requested_IQI_type_map()
        {
            // get requested IQI names from the model
            const moris::Cell< std::string > & tIQINames =
                    mEquationModel->get_requested_IQI_names();

            // clear the requested IQI assembly map
            mRequestedIQINamesAssemblyMap.clear();

            // loop over the requested IQI names
            for( uint Ik = 0; Ik < tIQINames.size(); Ik++ )
            {
                // fill the IQI assembly map
                mRequestedIQINamesAssemblyMap[ tIQINames( Ik ) ] = Ik;
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell < enum PDV_Type > Equation_Set::get_requested_dv_types()
        {
            moris::Cell< enum PDV_Type > tDvTypes;
            mEquationModel->get_design_variable_interface()->get_ip_requested_dv_types( tDvTypes );
            return tDvTypes;
        }

        //------------------------------------------------------------------------------
        // FIXME this might be too slow.
        moris_index Equation_Set::get_QI_assembly_index( const std::string & aIQIName )
        {
            return mRequestedIQINamesAssemblyMap.find( aIQIName );
        }

        //-------------------------------------------------------------------------------------------------

    }/* end_namespace_msi */
}/* end_namespace_moris */

