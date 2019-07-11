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
        { }

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
            // fixme param coeff from cluster
            mSet->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster,
                                                                                                                                  tSideOrd ) );
            mSet->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                // FIXME
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the number of active Dof_type for the ith IWG
                uint tNumOfIWGActiveDof = tTreatedIWG->get_active_dof_types().size();

                // get the field interpolators for the ith IWG in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_field_interpolator() );

                //get number of integration points
                uint tNumOfIntegPoints = mSet->get_num_integration_points();

                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                    // get integration point location in the reference volume
                    Matrix< DDRMat > tGlobalIntegPoint = mSet->get_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint );

                    // set integration point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                    }

                    // compute the integration point weight
                    real tWStar = mSet->get_integration_weights()( iGP )
                                * mSet->get_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                    // get the normal from mesh and set if for the IWG
                    Matrix< DDRMat > tNormal = mCluster->get_side_normal( mMasterCell, tSideOrd, tLocalIntegPoint );
                    tTreatedIWG->set_normal( tNormal );

                    // compute residual at integration point
                    Matrix< DDRMat > tResidual;
                    tTreatedIWG->compute_residual( tResidual, tIWGInterpolators );

                    // get the index of the residual dof type for the ith IWG in the list of element dof type
                    uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                    // get location of computed residual in global element residual
                    uint startDof = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                    uint stopDof  = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                    // add contribution to jacobian from evaluation point
                    mSet->mResidual( { startDof, stopDof }, { 0, 0 } )
                        = mSet->mResidual( { startDof, stopDof }, { 0, 0 } ) + tResidual * tWStar;
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

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mSet->get_IWGs()( iIWG );

                // FIXME
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_field_interpolator() );

                //get number of integration points
                uint tNumOfIntegPoints = mSet->get_num_integration_points();

                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get integration point location in the reference surface
                    Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );
                    //print(tLocalIntegPoint,"tLocalIntegPoint");

                    // get integration point location in the reference volume
                    Matrix< DDRMat > tGlobalIntegPoint = mSet->get_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint );
                    //print(tGlobalIntegPoint,"tGlobalIntegPoint");

                    // set integration point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                    }

                    // compute integration point weight
                    real tWStar = mSet->get_integration_weights()( iGP )
                                * mSet->get_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                    // get the normal from mesh and set if for the IWG
                    Matrix< DDRMat > tNormal = mMasterCell->compute_outward_side_normal( tSideOrd );
                    tTreatedIWG->set_normal( tNormal );

                    // compute jacobian at evaluation point
                    Cell< Matrix< DDRMat > > tJacobians;
                    tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                    // get location of computed jacobian in global element residual rows
                    uint startIDof = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                    uint stopIDof  = mSet->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                    // loop over the IWG active dof types
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                    {
                        // get the index of the active dof type
                        uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                        // get location of computed jacobian in global element residual columns
                        uint startJDof = mSet->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 0 );
                        uint stopJDof  = mSet->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 1 );

                        // add contribution to jacobian from evaluation point
                        mSet->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                            = mSet->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                            + tWStar * tJacobians( iIWGFI );
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
