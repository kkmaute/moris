#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"      //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset( mtk::Cell const * aCell,
                                          Set             * aElementBlock,
                                          Cluster         * aCluster) : Element( aCell, aElementBlock, aCluster )
        { }

//------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Sideset::compute_residual()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mCluster->mListOfSideOrdinals.numel();

            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_block_IG_geometry_interpolator()->set_space_coeff( mCell->get_vertex_coords());
            mSet->get_block_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );
            //print(mCell->get_vertex_coords(),"mCell->get_vertex_coords()");
            //print(mCluster->mTime,"mCluster->mTime");

            // set the geometry interpolator param space and time coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_block_IG_geometry_interpolator()->set_param_coeff();

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get the treatedSideOrdinal
                moris_index tTreatedSideOrdinal = mCluster->mListOfSideOrdinals( iSideset );
                //std::cout<<tTreatedSideOrdinal<<std::endl;

                // get the side phys and param coords
                // fixme remove as geo interp built on the side
                mSet->get_block_IG_geometry_interpolator()->build_space_side_space_phys_coeff( tTreatedSideOrdinal );
                mSet->get_block_IG_geometry_interpolator()->build_space_side_space_param_coeff( tTreatedSideOrdinal );

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
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_block_field_interpolator() );

                    //get number of integration points
                    uint tNumOfIntegPoints = mSet->get_num_integration_points();

                    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                    {
                        // get integration point location in the reference surface
                        Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                        // get integration point location in the reference volume
                        Matrix< DDRMat > tGlobalIntegPoint = mSet->get_block_IG_geometry_interpolator()->surf_val( tLocalIntegPoint );

                        // set integration point
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                        {
                            tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                        }

                        // compute the integration point weight
                        real tWStar = mSet->get_integration_weights()( iGP )
                                    * mSet->get_block_IG_geometry_interpolator()->surf_det_J( tLocalIntegPoint );

                        // compute the normal at integration point and set it for IWG
                        Matrix< DDRMat > tNormal = mSet->get_block_IG_geometry_interpolator()->surf_normal( tLocalIntegPoint );
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
            }
//            // print residual for check
//            print( mSet->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mCluster->mListOfSideOrdinals.numel();

            // set the geometry interpolator space and time physical coefficients for integration cell
            //moris::print(mCell->get_cell_physical_coords_on_side_ordinal(mCluster->mListOfSideOrdinals(0)),"coords on side");
            // FIXME: EITHER REMOVE LOOP OR MOVE THIS SET INTO LOOP
            //mSet->get_block_IG_geometry_interpolator()->set_space_coeff( mCell->get_cell_physical_coords_on_side_ordinal(mCluster->mListOfSideOrdinals(0)) );
            mSet->get_block_IG_geometry_interpolator()->set_space_coeff( mCell->get_vertex_coords() );
            mSet->get_block_IG_geometry_interpolator()->set_time_coeff( mCluster->mTime );

            // set the geometry interpolator space and time param coefficients for integration cell
            // fixme param coeff from cluster
            mSet->get_block_IG_geometry_interpolator()->set_param_coeff();

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
//                 get the treated side ordinal
                moris_index tTreatedSideOrdinal = mCluster->mListOfSideOrdinals( iSideset );
                //std::cout<<tTreatedSideOrdinal<<std::endl;

                // get the side phys and param coords
                // fixme remove as geo interp built on the side
                mSet->get_block_IG_geometry_interpolator()->build_space_side_space_phys_coeff( tTreatedSideOrdinal );
                mSet->get_block_IG_geometry_interpolator()->build_space_side_space_param_coeff( tTreatedSideOrdinal );

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
                        = mSet->get_IWG_field_interpolators( tTreatedIWG, mSet->get_block_field_interpolator() );

                    //get number of integration points
                    uint tNumOfIntegPoints = mSet->get_num_integration_points();

                    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                    {
                        // get integration point location in the reference surface
                        Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                        // get integration point location in the reference volume
                        Matrix< DDRMat > tGlobalIntegPoint = mSet->get_block_IG_geometry_interpolator()->surf_val( tLocalIntegPoint );

                        // set integration point
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                        {
                            tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                        }

                        // compute integration point weight
                        real tWStar = mSet->get_integration_weights()( iGP )
                                    * mSet->get_block_IG_geometry_interpolator()->surf_det_J( tLocalIntegPoint );

                        // evaluate the normal
                        // fixme should come from the mesh
                        Matrix< DDRMat > tNormal = mSet->get_block_IG_geometry_interpolator()->surf_normal( tLocalIntegPoint );
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
