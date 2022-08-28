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
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterDofTypes;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveDofTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_type_list - can only be MASTER or SLAVE");
                    return mMasterDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        Matrix< DDSMat > & Equation_Set::get_dof_type_map(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterDofTypeMap;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveDofTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return mMasterDofTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dof_index_for_type(
                enum MSI::Dof_Type aDofType,
                mtk::Master_Slave  aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // check if master dof type in list
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type - master dof type does not exist in map." );

                    // get the index for master dof type from map
                    sint tMasterIndex = mMasterDofTypeMap( static_cast< int >( aDofType ) );

                    //                    // check if master dof type assigned in map
                    //                    MORIS_ASSERT( tMasterIndex != -1,
                    //                                  "Equation_Set::get_dof_index_for_type - master dof type not assigned in map." );

                    // return master dof type index
                    return tMasterIndex;
                }

                case mtk::Master_Slave::SLAVE:
                {
                    // check if slave dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type - slave dof type does not exist in map." );

                    // get the index for slave dof type from map
                    sint tSlaveIndex = mSlaveDofTypeMap( static_cast< int >( aDofType ) );

                    //                    // check if slave dof type assigned in map
                    //                    MORIS_ASSERT( tSlaveIndex != -1,
                    //                                  "Equation_Set::get_dof_index_for_type - slave dof type not assigned in map." );

                    if ( tSlaveIndex == -1 )
                    {
                        return tSlaveIndex;
                    }
                    else
                    {
                        // get maximum index from master map
                        moris::sint tMaxMasterIndex = mMasterDofTypeMap.max();

                        // check if master map was assigned
                        MORIS_ASSERT( tMaxMasterIndex != -1,
                                "Equation_Set::get_dof_index_for_type - mMasterDofTypeMap is empty." );

                        // return slave dof type index
                        return tSlaveIndex + tMaxMasterIndex + 1;
                    }
                }

                default:
                {
                    MORIS_ERROR( false,
                            "Equation_Set::get_dof_index_for_type - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dof_index_for_type_1(
                enum MSI::Dof_Type     aDofType,
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // check if master dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type_1 - master dof type does not exist in map." );

                    //                    // check if master dof type assigned in map
                    //                    MORIS_ASSERT( mMasterDofTypeMap( static_cast< int >( aDofType ) ) != -1,
                    //                                  "Equation_Set::get_dof_index_for_type_1 - master dof type does not assigned in map." );

                    // return index for master dof type
                    return mMasterDofTypeMap( static_cast< int >( aDofType ) );
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // check if slave dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(),
                            "Equation_Set::get_dof_index_for_type_1 - slave dof type does not exist in map." );

                    //                    // check if slave dof type assigned in map
                    //                    MORIS_ASSERT( mSlaveDofTypeMap( static_cast< int >( aDofType ) ) != -1,
                    //                                  "Equation_Set::get_dof_index_for_type_1 - slave dof type not assigned in map." );

                    // return index for slave dof type
                    return mSlaveDofTypeMap( static_cast< int >( aDofType ) );
                }

                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_index_for_type_1 - can only be MASTER or SLAVE." );
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > > & Equation_Set::get_dv_type_list(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterDvTypes;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveDvTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dv_type_list - can only be MASTER or SLAVE.");
                    return mMasterDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Equation_Set::get_dv_type_map(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterDvTypeMap;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveDvTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dv_type_map - can only be MASTER or SLAVE");
                    return mMasterDvTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< mtk::Field_Type > > & Equation_Set::get_field_type_list(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterFieldTypes;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveFieldTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_field_type_list - can only be MASTER or SLAVE.");
                    return mMasterFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Equation_Set::get_field_type_map(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterFieldTypeMap;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlaveFieldTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_field_type_map - can only be MASTER or SLAVE");
                    return mMasterFieldTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dv_index_for_type(
                enum PDV_Type     aDvType,
                mtk::Master_Slave aIsMaster  )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type - master dv type does not exist in map." );

                    //                    // check if dv type assigned in map
                    //                    MORIS_ASSERT( mMasterDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type - master dv type not assigned in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }

                case mtk::Master_Slave::SLAVE:
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type(), slave dv type does not exist in map." );

                    //                    // check if dv type assigned in map
                    //                    MORIS_ASSERT( mSlaveDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type - slave dv type not assigned in map." );

                    // get the set index for dv type
                    sint tSlaveIndex = mSlaveDvTypeMap( static_cast< int >( aDvType ) );

                    // if index is -1
                    if ( tSlaveIndex == -1 )
                    {
                        return tSlaveIndex;
                    }
                    else
                    {
                        // get the max set index for dv types
                        moris::sint tMaxMasterIndex = mMasterDvTypeMap.max();

                        // check if mMasterDvTypeMap is set
                        MORIS_ASSERT( tMaxMasterIndex != -1,
                                "Equation_Set::get_dv_index_for_type - mMasterDvTypeMap is empty." );

                        // return set index for dv type
                        return tSlaveIndex + tMaxMasterIndex + 1;
                    }
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_dv_index_for_type - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_dv_index_for_type_1(
                enum PDV_Type            aDvType,
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type_1 - master dv type does not exist in map." );

                    //                    // check if dv type is set in map
                    //                    MORIS_ASSERT( mMasterDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type_1 - master dv type not assigned in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                            "Equation_Set::get_dv_index_for_type_1 - slave dv type does not exist in map." );

                    //                    // check if dv type is set in map
                    //                    MORIS_ASSERT( mSlaveDvTypeMap( static_cast< int >( aDvType ) ),
                    //                                  "Equation_Set::get_dv_index_for_type_1 - slave dv type not assigned in map." );

                    // return set index for dv type
                    return mSlaveDvTypeMap( static_cast< int >( aDvType ) );
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_dv_index_for_type_1 - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

        //------------------------------------------------------------------------------

        sint Equation_Set::get_field_index_for_type_1(
                enum mtk::Field_Type aFieldType,
                mtk::Master_Slave    aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aFieldType ) < mMasterFieldTypeMap.numel(),
                            "Equation_Set::get_field_index_for_type_1 - master field type does not exist in map." );

                    //// check if field type is set in map
                    //MORIS_ASSERT( mMasterDvTypeMap( static_cast< int >( aFieldType ) ),
                    //          "Equation_Set::get_field_index_for_type_1 - master field type not assigned in map." );

                    // return set index for dv type
                    return mMasterFieldTypeMap( static_cast< int >( aFieldType ) );
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // check if field type is set in map
                    MORIS_ASSERT( static_cast< uint >( aFieldType ) < mSlaveFieldTypeMap.numel(),
                            "Equation_Set::get_field_index_for_type_1 - slave dv type does not exist in map." );

                    //// check if field type is set in map
                    //MORIS_ASSERT( mSlaveFieldTypeMap( static_cast< int >( aDvType ) ),
                    //        "Equation_Set::get_field_index_for_type_1 - slave field type not assigned in map." );

                    // return set index for dv type
                    return mSlaveFieldTypeMap( static_cast< int >( aFieldType ) );
                }
                default:
                {
                    MORIS_ERROR(false,
                            "Equation_Set::get_field_index_for_type_1 - can only be MASTER or SLAVE.");
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

