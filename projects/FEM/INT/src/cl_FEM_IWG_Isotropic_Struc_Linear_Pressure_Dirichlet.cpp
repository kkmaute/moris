#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::IWG_Isotropic_Struc_Linear_Pressure_Dirichlet()
        {

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;
            mPropertyMap[ "Select" ] = IWG_Property_Type::SELECT;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = IWG_Stabilization_Type::DIRICHLET_NITSCHE;

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get SP, CM and property indices
            uint tSelectIndex      = static_cast< uint >( IWG_Property_Type::SELECT );
            uint tDirichletIndex   = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
//            uint tNitscheIndex     = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // get CM
            moris::fem::CM_Struc_Linear_Isotropic* tLinearIso = static_cast <moris::fem::CM_Struc_Linear_Isotropic*> (mMasterCM(tElastLinIsoIndex).get());

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix if needed
            if ( mMasterProp( tSelectIndex ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();
                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( tSelectIndex )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

            // flatten normal
            Matrix< DDRMat > tFlattenedNormal(1, 1, 0.0);
            tLinearIso->flatten_normal(mNormal, tFlattenedNormal);

            // get start and end indices for residual assembly
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
            += trans(tLinearIso->dFluxdDOF( moris::Cell<MSI::Dof_Type> (1, MSI::Dof_Type::P) )) * trans(tFlattenedNormal) * tM * tJump * aWStar;

//            moris::print(trans(tLinearIso->dFluxdDOF( moris::Cell<MSI::Dof_Type> (1, MSI::Dof_Type::P) )), "dfluxddof");
//            moris::print(trans(tFlattenedNormal), "normal");
//            moris::print(tM, "tm");
//            moris::print(tJump, "jump");
//            std::cout << aWStar << std::endl;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for a given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for a given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get SP, CM and property indices
            uint tSelectIndex      = static_cast< uint >( IWG_Property_Type::SELECT );
            uint tDirichletIndex   = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
//            uint tNitscheIndex     = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // get CM
            moris::fem::CM_Struc_Linear_Isotropic* tLinearIso = static_cast <moris::fem::CM_Struc_Linear_Isotropic*> (mMasterCM(tElastLinIsoIndex).get());

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix
            if ( mMasterProp( tSelectIndex ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( tSelectIndex )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

            // flatten normal
            Matrix< DDRMat > tFlattenedNormal(1, 1, 0.0);
            tLinearIso->flatten_normal(mNormal, tFlattenedNormal);

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for this dof type
                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if direct dependency on the dof type
                if (tDofType(0) == MSI::Dof_Type::UX)
                {
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += trans(tLinearIso->dFluxdDOF( moris::Cell<MSI::Dof_Type> (1, MSI::Dof_Type::P) )) * trans(tFlattenedNormal) * tM * tFI->N() * aWStar;
                }

                // if property depends on dof type
                if ( mMasterProp( tDirichletIndex )->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                    += trans(tLinearIso->dFluxdDOF( moris::Cell<MSI::Dof_Type> (1, MSI::Dof_Type::P) )) * trans(tFlattenedNormal) * tM * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType ) * aWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Dirichlet::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
