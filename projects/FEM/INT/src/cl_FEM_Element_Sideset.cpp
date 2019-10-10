#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
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
            uint tNumProp = mSet->get_number_of_properties();
            uint tNumCM   = mSet->get_number_of_constitutive_models();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_IWGs();

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
                    mSet->get_field_interpolators()( iFI )->set_space_time( tGlobalIntegPoint );
                }

                // reset properties
                for ( uint iProp = 0; iProp < tNumProp; iProp++ )
                {
                    mSet->get_properties()( iProp )->reset_eval_flags();
                }

                // reset constitutive models
                for ( uint iCM = 0; iCM < tNumCM; iCM++ )
                {
                    mSet->get_constitutive_models()( iCM )->reset_eval_flags();
                }

                // compute the integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // FIXME
                    mSet->get_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_IWGs()( iIWG )->set_normal( tNormal );

                    // compute residual at integration point
                    moris::Cell< Matrix< DDRMat > > tResidual;
                    mSet->get_IWGs()( iIWG )->compute_residual( tResidual );

                    // add contribution to jacobian from evaluation point
                    mSet->mResidual( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 0 ), mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 1 ) },
                                     { 0, 0 } ) += tWStar * tResidual( 0 );
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
            uint tNumProp = mSet->get_number_of_properties();
            uint tNumCM   = mSet->get_number_of_constitutive_models();

            // get number of IWGs
            uint tNumIWGs = mSet->get_number_of_IWGs();

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
                    mSet->get_field_interpolators()( iFI )->set_space_time( tGlobalIntegPoint );
                }

                // reset properties
                for ( uint iProp = 0; iProp < tNumProp; iProp++ )
                {
                    mSet->get_properties()( iProp )->reset_eval_flags();
                }

                // reset constitutive models
                for ( uint iCM = 0; iCM < tNumCM; iCM++ )
                {
                    mSet->get_constitutive_models()( iCM )->reset_eval_flags();
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // get the normal from mesh
                Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
                {
                    // FIXME set BCs
                    mSet->get_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // set the normal for the IWG
                    mSet->get_IWGs()( iIWG )->set_normal( tNormal );

                    // compute jacobian at evaluation point
                    moris::Cell< moris::Cell< Matrix< DDRMat > > > tJacobians;
                    mSet->get_IWGs()( iIWG )->compute_jacobian( tJacobians );
//                    print( tJacobians(0), "tJacobians" );
//
//                    // check with finite difference
//                    real tPerturbation = 1E-6;
//                    Cell< Matrix< DDRMat > > tJacobiansFD;
//                    mSet->get_IWGs()( iIWG )->compute_jacobian_FD( tJacobiansFD,
//                                                                   mSet->get_IWG_field_interpolators()( iIWG ),
//                                                                   tPerturbation );
//                    print(tJacobiansFD(0),"tJacobiansFD");

                    // loop over the IWG active dof types
                    uint tNumIWGDof = mSet->get_IWGs()( iIWG )->get_global_dof_type_list().size();
                    for ( uint iIWGFI = 0; iIWGFI < tNumIWGDof; iIWGFI++)
                    {
                        // add contribution to jacobian from evaluation point
                        mSet->mJacobian( { mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 0 ),      mSet->get_IWG_res_dof_assembly_map()( iIWG )( 0, 1 ) },
                                         { mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 0 ), mSet->get_IWG_jac_dof_assembly_map()( iIWG )( iIWGFI, 1 ) } )
                                       += tWStar * tJacobians( 0 )( iIWGFI );
                    }
                }
            }
//            // print jacobian for check
//            print( mCluster->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
