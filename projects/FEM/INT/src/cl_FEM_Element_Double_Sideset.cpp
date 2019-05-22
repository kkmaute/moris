#include <iostream>

#include "cl_FEM_Element_Double_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Double_Sideset::Element_Double_Sideset( mtk::Cell    const * aMasterIGCell,
                                                        mtk::Cell    const * aSlaveIGCell,
                                                        Element_Block      * aElementBlock,
                                                        Cluster            * aCluster) : Element( aCell, aElementBlock, aCluster )
        { }

//------------------------------------------------------------------------------

        Element_Double_Sideset::~Element_Double_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_residual()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mCluster->mListOfSideOrdinals.numel();

            // set the geometry interpolator coefficients
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mCluster->mTime );
            // fixme param coeff from cluster
            mElementBlock->get_block_geometry_interpolator()->set_param_coeff();

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get the treatedSideOrdinal
                moris_index tTreatedSideOrdinal = mCluster->mListOfSideOrdinals( iSideset );

                // get the side phys and param coords
                mElementBlock->get_block_geometry_interpolator()->build_space_side_space_phys_coeff( tTreatedSideOrdinal );
                mElementBlock->get_block_geometry_interpolator()->build_space_side_space_param_coeff( tTreatedSideOrdinal );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the treated IWG
                    IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                    // FIXME
                    tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // get the number of active Dof_type for the ith IWG
                    uint tNumOfIWGActiveDof = tTreatedIWG->get_active_dof_types().size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = mElementBlock->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                    //get number of integration points
                    uint tNumOfIntegPoints = mElementBlock->get_num_integration_points();

                    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                    {
                        // get integration point location in the reference surface
                        Matrix< DDRMat > tSurfRefIntegPointI = mElementBlock->get_integration_points().get_column( iGP );

                        // get integration point location in the reference volume
                        Matrix< DDRMat > tVolRefIntegPointI
                            = mElementBlock->get_block_geometry_interpolator()->surf_val( tSurfRefIntegPointI );

                        // set integration point
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                        {
                            tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                        }

                        // compute integration point weight x detJ
                        real tSurfDetJ = mElementBlock->get_block_geometry_interpolator()->surf_det_J( tSurfRefIntegPointI );
                        Matrix< DDRMat > tNormal = mElementBlock->get_block_geometry_interpolator()->surf_normal( tSurfRefIntegPointI );

                        tTreatedIWG->set_normal( tNormal );

                        real tWStar = mElementBlock->get_integration_weights()( iGP ) * tSurfDetJ;

                        // compute jacobian at evaluation point
                        Matrix< DDRMat > tResidual;
                        tTreatedIWG->compute_residual( tResidual, tIWGInterpolators );

                        // get the index of the residual dof type for the ith IWG in the list of element dof type
                        uint tIWGResDofIndex
                            = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                        // get location of computed residual in global element residual
                        uint startDof = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                        uint stopDof  = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                        // add contribution to jacobian from evaluation point
                        mCluster->mResidual( { startDof, stopDof }, { 0, 0 } )
                            = mCluster->mResidual( { startDof, stopDof }, { 0, 0 } ) + tResidual * tWStar;
                    }
                }
            }
//            // print residual for check
//            print( mCluster->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_jacobian()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mCluster->mListOfSideOrdinals.numel();

            // set the geometry interpolator coefficients
            //FIXME: tHat are set by default but should come from solver
//            Matrix< DDRMat > tTHat = { {0.0}, {1.0} };
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mCluster->mTime );
            // fixme param coeff from cluster
            mElementBlock->get_block_geometry_interpolator()->set_param_coeff();

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get the treated side ordinal
                moris_index tTreatedSideOrdinal = mCluster->mListOfSideOrdinals( iSideset );

                // get the side phys and param coords
                mElementBlock->get_block_geometry_interpolator()->build_space_side_space_phys_coeff( tTreatedSideOrdinal );
                mElementBlock->get_block_geometry_interpolator()->build_space_side_space_param_coeff( tTreatedSideOrdinal );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the treated IWG
                    IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                    // FIXME
                    tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // get the index of the residual dof type for the ith IWG in the list of element dof type
                    uint tIWGResDofIndex
                        = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                    Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                    uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = mElementBlock->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                    //get number of integration points
                    uint tNumOfIntegPoints = mElementBlock->get_num_integration_points();

                    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                    {
                        // get integration point location in the reference surface
                        Matrix< DDRMat > tSurfRefIntegPointI = mElementBlock->get_integration_points().get_column( iGP );

                        // get integration point location in the reference volume
                        Matrix< DDRMat > tVolRefIntegPointI = mElementBlock->get_block_geometry_interpolator()->surf_val( tSurfRefIntegPointI );

                        // set integration point
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                        {
                            tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                        }

                        // compute integration point weight x detJ
//                        real tSurfDetJ;
//                        Matrix< DDRMat > tNormal;
//                        mElementBlock->get_block_geometry_interpolator()->surf_det_J( tSurfDetJ,
//                                                           tNormal,
//                                                           tSurfRefIntegPointI,
//                                                           tTreatedSideOrdinal );
                        real tSurfDetJ = mElementBlock->get_block_geometry_interpolator()->surf_det_J( tSurfRefIntegPointI );
                        real tWStar = mElementBlock->get_integration_weights()( iGP ) * tSurfDetJ;

                        // evaluate the normal
                        Matrix< DDRMat > tNormal = mElementBlock->get_block_geometry_interpolator()->surf_normal( tSurfRefIntegPointI );
                        tTreatedIWG->set_normal( tNormal );

                        // compute jacobian at evaluation point
                        Cell< Matrix< DDRMat > > tJacobians;
                        tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                        // get location of computed jacobian in global element residual rows
                        uint startIDof = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                        uint stopIDof  = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                        // loop over the IWG active dof types
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                        {
                            // get the index of the active dof type
                            uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                            // get location of computed jacobian in global element residual columns
                            uint startJDof = mElementBlock->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 0 );
                            uint stopJDof  = mElementBlock->get_interpolator_dof_assembly_map()( tIWGActiveDofIndex, 1 );

                            // add contribution to jacobian from evaluation point
                            mCluster->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                                = mCluster->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                                + tWStar * tJacobians( iIWGFI );
                        }
                    }
                }
            }
//            // print jacobian for check
//            print( mCluster->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Double_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Double_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
