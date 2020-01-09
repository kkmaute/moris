
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
    //------------------------------------------------------------------------------
        IWG_Isotropic_Struc_Linear_Bulk::IWG_Isotropic_Struc_Linear_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ] = IWG_Property_Type::LOAD;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( real tWStar )
        {

#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for this residual (displacement)
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX);

            // get start and end index for residual assembly
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get property, CM, SP indices
            uint tLoadIndex        = static_cast< uint >( IWG_Property_Type::LOAD );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
            += trans( mMasterCM( tElastLinIsoIndex )->testStrain() ) * mMasterCM( tElastLinIsoIndex )->flux() * tWStar;

            // if body load
            if ( mMasterProp( tLoadIndex ) != nullptr )
            {
                // compute body load contribution
                mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                += - trans( tDisplacementFI->N() ) * mMasterProp( tLoadIndex )->val()( 0 ) * tWStar;
            }

            // if pressure dof
//            uint tElastLinIsoPressureIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO_PRESSURE );
//            if (mMasterCM(tElastLinIsoPressureIndex) != nullptr)
//            {
//                mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
//                -= trans(mMasterCM(tElastLinIsoIndex)->testStrain()) * mMasterCM(tElastLinIsoPressureIndex)->flux() * tWStar;
//            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for given dof type
            Field_Interpolator * tDisplacementFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ));

            // get property, CM, SP indices
            uint tLoadIndex        = static_cast< uint >( IWG_Property_Type::LOAD );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if body load
                if ( mMasterProp( tLoadIndex ) != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( tLoadIndex )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                        += - trans( tDisplacementFI->N() ) * mMasterProp( tLoadIndex )->dPropdDOF( tDofType ) * tWStar;
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                    += trans( mMasterCM( tElastLinIsoIndex )->testStrain() ) * mMasterCM( tElastLinIsoIndex )->dFluxdDOF( tDofType ) * tWStar;
                }

                //if (mMasterCM(tElastLinIsoPressureIndex) != nullptr && mMasterCM( tElastLinIsoPressureIndex )->check_dof_dependency( tDofType ))
//                if (mMasterCM(tElastLinIsoPressureIndex) != nullptr && tIndexDep != 0)
//                {
//                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
//                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
//                                          -= trans(mMasterCM(tElastLinIsoIndex)->testStrain()) * mMasterCM(tElastLinIsoPressureIndex)->strain() * tWStar;
//                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                             moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
