//FEM/INT/src
#include "cl_FEM_IWG_Heat_Method_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Heat_Method_Bulk::IWG_Heat_Method_Bulk(){}

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                std::string tErrMsg =
                        "IWG_Heat_Method_Bulk::set_property - Unknown aPropertyString : " +
                        aPropertyString;

                MORIS_ERROR( false , tErrMsg.c_str() );
            }

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Heat_Method_Bulk::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            if ( mConstitutiveMap.find( aConstitutiveString ) == mConstitutiveMap.end() )
            {
                std::string tErrMsg =
                        "IWG_Heat_Method_Bulk::set_constitutive_model - Unknown aConstitutiveString: " +
                        aConstitutiveString;

                MORIS_ERROR( false, tErrMsg.c_str() );
            }

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Heat_Method_Bulk::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Heat_Method_Bulk::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get residual dof type field interpolator (here temperature)
            Field_Interpolator * tFITheta = mMasterFIManager->get_field_interpolators_for_type( mDofTheta );
            Field_Interpolator * tFIPhiD = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // compute normalized gradient of Theta
            Matrix< DDRMat > tGradTheta = tFITheta->gradx( 1 );
            real tNorm = norm( tGradTheta );
            Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            trans( tFIPhiD->dnNdxn( 1 ) ) * tFIPhiD->gradx( 1 ) +
                            trans( tFIPhiD->dnNdxn( 1 ) ) * tNormGradTheta );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Heat_Method_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get residual dof type field interpolator (here temperature)
            Field_Interpolator * tFITheta = mMasterFIManager->get_field_interpolators_for_type( mDofTheta );
            Field_Interpolator * tFIPhiD = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is density, add diagonal term (density - density DoF types)
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // add contibution
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIPhiD->dnNdxn( 1 ) ) * tFIPhiD->dnNdxn( 1 ) );
                }

                // if dof type is density, add diagonal term (density - density DoF types)
                if( tDofType( 0 ) == mDofTheta )
                {
                    // compute normalized gradient of Theta
                    Matrix< DDRMat > tBTheta = tFITheta->dnNdxn( 1 );
                    Matrix< DDRMat > tGradTheta = tFITheta->gradx( 1 );
                    real tNorm = norm( tGradTheta );
                    Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;
                    Matrix< DDRMat > tNormBTheta = tBTheta / tNorm;

                    // add contibution
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFIPhiD->dnNdxn( 1 ) ) * tNormBTheta -
                                    trans( tFIPhiD->dnNdxn( 1 ) ) *  tNormGradTheta * trans( tNormGradTheta ) * tNormBTheta );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Heat_Method_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Heat_Method_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Heat_Method_Bulk::compute_dRdp - Not implemented." );

        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
