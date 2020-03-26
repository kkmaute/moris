/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
 *
 *  Created on: Jan 27, 2020
 *      Author: ritzert
 */
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        IWG_Isotropic_Struc_Linear_Contact_Penalty::IWG_Isotropic_Struc_Linear_Contact_Penalty()
        {
            // set size for the constitutive model pointer cell
            // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "PenaltyContact" ]      = IWG_Stabilization_Type::NITSCHE_INTERFACE;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual( real tWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get indices for SP, CM, property
//            uint tPenaltyIndex      = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );

            moris::real tTempParam = 1e13;

            // evaluate gap
            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Master side

            if ( tGap(0,0) > 0 )
            {
                tGap(0,0) = 0;
            }
                // compute contact residual on slave side
                mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } )
                += (1) * tTempParam* trans( tFISlave -> N()) * mNormal * tGap(0,0) * tWStar;

                // compute contact residual on master side
                mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } )
                += (-1) * tTempParam* trans( tFIMaster -> N()) * mNormal * tGap(0,0) * tWStar;


//            moris::print(mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } ),"residual_s");
//            moris::print(mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } ),"residual_m");
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get indices for SP, CM, property
//            uint tPenaltyIndex      = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );

            // create quadrature of normal
            Matrix< DDRMat > tNormQ = mNormal * trans (mNormal);

            moris::real tTempParam = 1e13;

            // evaluate gap
            Matrix< DDRMat > tGap = trans(mNormal) * ( tFISlave->val() - tFIMaster->val()); // mNormal is normal on Master side

            if ( tGap(0,0) <= 0 && mResidualDofTypeRequested)
            {
                // compute the jacobian on slave side
//                std::cout<<"tDofIndexSlave   "<<tDofIndexSlave<<std::endl;
//                print(mSet->get_res_dof_assembly_map(),"mSet->get_res_dof_assembly_map()");
//                print(mSet->get_jac_dof_assembly_map(),"mSet->get_jac_dof_assembly_map()");
                // slave slave
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                        { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) } )
                        += tTempParam * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() * tWStar;

                // slave master
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                        { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) } )
                        += (-1) * tTempParam * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() * tWStar;

                // compute the jacobian on master side
                //master master
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                        { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) } )
                        += tTempParam* trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() * tWStar;

                //master slave
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                        { mSet->get_jac_dof_assembly_map()( tDofIndexMaster  )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } )
                        += (-1) * tTempParam * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() * tWStar;

                //    moris::print(mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                //            { mSet->get_jac_dof_assembly_map()( tDofIndexMaster  )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } ),"Jacobian4");

            }
        }
//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual - This function does nothing.");
        }


//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_drdpdv( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_drdpdv - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


