
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    IWG_Isotropic_Struc_Linear_Interface::IWG_Isotropic_Struc_Linear_Interface()
        {
            // FIXME set penalty parameter
            mGammaInterface = 1.0;

            // FIXME set weight master parameter
            mMasterWeight = 0.5;

            // FIXME set weight slave parameter
            mSlaveWeight = 0.5;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_residual( real tWStar )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

//            // check master and slave properties
//            this->check_properties( mtk::Master_Slave::MASTER );
//            this->check_properties( mtk::Master_Slave::SLAVE );
//
//            // check master and slave constitutive models
//            this->check_constitutive_models( mtk::Master_Slave::MASTER );
//            this->check_constitutive_models( mtk::Master_Slave::SLAVE );


            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            Field_Interpolator * tFIMaster = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            Field_Interpolator * tFISlave  = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // evaluate average traction
            Matrix< DDRMat > tTraction = mMasterWeight * mMasterCM( 0 )->traction( mNormal ) + mSlaveWeight * mSlaveCM( 0 )->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } )
                          += ( - trans( tFIMaster->N() ) * tTraction
                             + mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mGammaInterface * trans( tFIMaster->N() ) * tJump )                    *tWStar;

            // compute slave residual
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } )
                             +=  ( trans( tFISlave->N() ) * tTraction
                             + mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * tJump
                             - mGammaInterface * trans( tFISlave->N() ) * tJump )                    *tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian( real tWStar )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

//            // check master and slave properties
//            this->check_properties( mtk::Master_Slave::MASTER );
//            this->check_properties( mtk::Master_Slave::SLAVE );
//
//            // check master and slave constitutive models
//            this->check_constitutive_models( mtk::Master_Slave::MASTER );
//            this->check_constitutive_models( mtk::Master_Slave::SLAVE );

            Field_Interpolator * tFIMaster = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            Field_Interpolator * tFISlave  = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // evaluate displacement jump
            Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute the jacobian for direct dof dependencies
            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 2 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 3 ) } )
                    += (  mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * tFIMaster->N()
                        + mGammaInterface * trans( tFIMaster->N() ) * tFIMaster->N() ) * tWStar;

            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 0 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 2 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 3 ) } )
                    += ( - mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * tFISlave->N()
                         - mGammaInterface * trans( tFIMaster->N() ) * tFISlave->N() ) * tWStar;

            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 2 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 3 ) } )
                    += (  mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * tFIMaster->N()
                        - mGammaInterface * trans( tFISlave->N() ) * tFIMaster->N() ) * tWStar;

            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 2 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 3 ) } )
                    += ( - mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * tFISlave->N()
                         + mGammaInterface * trans( tFISlave->N() ) * tFISlave->N() ) * tWStar;

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) },
                                          { mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 2 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 3 ) } )
                            += ( - trans( tFIMaster->N() ) * mMasterWeight * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal )
                                 + mMasterWeight * mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) ) * tWStar;

                    mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) },
                                          { mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 2 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 3 ) } )
                              += trans( tFISlave->N() ) * mMasterWeight * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal ) * tWStar;
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );

                // if dependency on the dof type
                if ( mSlaveCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) },
                                          { mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 2 ), mSet->get_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 3 ) } )
                            += - trans( tFIMaster->N() ) * mSlaveWeight * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal ) * tWStar;

                    mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) },
                                          { mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 2 ), mSet->get_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 3 ) } )
                              += (   trans( tFISlave->N() ) * mSlaveWeight * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal )
                                   + mSlaveWeight * mSlaveCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) ) * tWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
