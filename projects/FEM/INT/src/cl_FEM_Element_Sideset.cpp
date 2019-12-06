#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset( mtk::Cell const    * aCell,
                                          Set                * aSet,
                                          Cluster            * aCluster,
                                          moris::moris_index   aCellIndexInCluster) : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Sideset::compute_residual()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( tSideOrd ) );
            mSet->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                  tSideOrd ) );
            mSet->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // get number of field interpolator and properties
            uint tNumFI   = mSet->get_number_of_field_interpolators();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // get integration point location in the reference volume
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for IP geometry interpolator
                mSet->get_IP_geometry_interpolator()->set_space_time( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iFI = 0; iFI < tNumFI; iFI++ )
                {
                    mSet->mFieldInterpolatorManager->get_field_interpolators_for_type( mSet->mMasterDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER )
                                                   ->set_space_time( tGlobalIntegPoint );
                }

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    mSet->get_requested_IWGs()( iIWG )->compute_residual( tWStar );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );
                }
            }
//            // print residual for check
//            print( mSet->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian()
        {
            // get treated side ordinal
            uint tSideOrd = mCluster->mMasterListOfSideOrdinals( mCellIndexInCluster );

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_cell_physical_coords_on_side_ordinal( tSideOrd ) );
            mSet->get_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                  tSideOrd ) );
            mSet->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme default

            // get number of field interpolator and properties
            uint tNumFI   = mSet->get_number_of_field_interpolators();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_requested_IWGs();

            // loop over integration points
            uint tNumIntegPoints = mSet->get_number_of_integration_points();
            for( uint iGP = 0; iGP < tNumIntegPoints; iGP++ )
            {
                // get integration point location in the reference surface
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // get integration point location in the reference volume
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for IP geometry interpolator
                mSet->get_IP_geometry_interpolator()->set_space_time( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iFI = 0; iFI < tNumFI; iFI++ )
                {
                    mSet->mFieldInterpolatorManager->get_field_interpolators_for_type( mSet->mMasterDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER )
                                                   ->set_space_time( tGlobalIntegPoint );
//                    mSet->get_field_interpolators()( iFI )->set_space_time( tGlobalIntegPoint );
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // reset IWG
                    mSet->get_requested_IWGs()( iIWG )->reset_eval_flags();

                    // FIXME set BCs
                    mSet->get_requested_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_requested_IWGs()( iIWG )->set_normal( tNormal );

                    // compute jacobian at evaluation point
                    mSet->get_requested_IWGs()( iIWG )->compute_jacobian( tWStar );

//                    // check with finite difference
//                    real tPerturbation = 1E-4;
//                    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;
//                    mSet->get_IWGs()( iIWG )->compute_jacobian_FD( tJacobiansFD,
//                                                                   tPerturbation );
//                    print(tJacobiansFD(0)(0),"tJacobiansFD");
                }
            }
//            // print jacobian for check
//            print( mSet->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
