#include <iostream>
#include "cl_FEM_Element_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Block.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset( mtk::Cell            const * aCell,
                                          moris::Cell< Node_Base* > & aNodes,
                                          Element_Block             * aElementBlock) : Cluster( aCell, aNodes, aElementBlock )
        {
        }

//------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset(){}

//------------------------------------------------------------------------------

        void Element_Sideset::compute_residual()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mListOfSideOrdinals.numel();

            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement();

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // set the geometry interpolator coefficients
            //FIXME: tHat are set by default but should come from solver
//            Matrix< DDRMat > tTHat = { {0.0}, {1.0} };
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mTime );

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get the treatedSideOrdinal
                moris_index tTreatedSideOrdinal = mListOfSideOrdinals( iSideset );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the treated IWG
                    IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                    // FIXME
                    tTreatedIWG->set_nodal_weak_bcs( this->get_weak_bcs() );

                    // get the number of active Dof_type for the ith IWG
                    uint tNumOfIWGActiveDof = tTreatedIWG->get_active_dof_types().size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = this->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                    // create an integration rule for the ith IWG
                    mtk::Geometry_Type tSideGeometryType = mElementBlock->get_block_geometry_interpolator()->get_side_geometry_type();
                    Integration_Rule tIntegrationRule( tSideGeometryType,
                                                       Integration_Type::GAUSS,
                                                       this->get_auto_integration_order( tSideGeometryType ),
                                                       Integration_Type::GAUSS,
                                                       Integration_Order::BAR_1 );

                    // create an integrator for the ith IWG
                    Integrator tIntegrator( tIntegrationRule );

                    //get number of integration points
                    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

                    // get integration points
                    Matrix< DDRMat > tSurfRefIntegPoints = tIntegrator.get_points();

                    // get integration weights
                    Matrix< DDRMat > tIntegWeights = tIntegrator.get_weights();

                    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                    {
                        // get integration point location in the reference surface
                        Matrix< DDRMat > tSurfRefIntegPointI = tSurfRefIntegPoints.get_column( iGP );

                        // get integration point location in the reference volume
                        Matrix< DDRMat > tVolRefIntegPointI
                            = mElementBlock->get_block_geometry_interpolator()->surf_val( tSurfRefIntegPointI, tTreatedSideOrdinal );

                        // set integration point
                        for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                        {
                            tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                        }

                        // compute integration point weight x detJ
                        real tSurfDetJ;
                        Matrix< DDRMat > tNormal;
                        mElementBlock->get_block_geometry_interpolator()->surf_det_J( tSurfDetJ,
                                                           tNormal,
                                                           tSurfRefIntegPointI,
                                                           tTreatedSideOrdinal );
                        tTreatedIWG->set_normal( tNormal );

                        real tWStar = tIntegWeights( iGP ) * tSurfDetJ;

                        // compute jacobian at evaluation point
                        Matrix< DDRMat > tResidual;
                        tTreatedIWG->compute_residual( tResidual, tIWGInterpolators );

                        // get the index of the residual dof type for the ith IWG in the list of element dof type
                        uint tIWGResDofIndex
                            = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                        // add contribution to jacobian from evaluation point
                        mResidualElement( tIWGResDofIndex )
                            = mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
                    }
                }
            }
            // residual assembly
            uint tCounterI = 0;
            uint startI, stopI;

            // loop over the field interpolators
            for ( uint iBuild = 0; iBuild < mElementBlock->get_num_interpolators(); iBuild++ )
            {
                // get the row position in the residual matrix
                startI = tCounterI;
                stopI  = tCounterI + mElementBlock->get_block_field_interpolator()( iBuild )->get_number_of_space_time_coefficients() - 1;

                // fill the global residual
                mResidual( { startI, stopI }, { 0 , 0 } ) = mResidualElement( iBuild ).matrix_data();

                // update the row counter
                tCounterI = stopI + 1;
            }
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian()
        {
            // get the number of side ordinals
            uint tNumOfSideSets = mListOfSideOrdinals.numel();

            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement();

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // set the geometry interpolator coefficients
            //FIXME: tHat are set by default but should come from solver
//            Matrix< DDRMat > tTHat = { {0.0}, {1.0} };
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mTime );

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get the treated side ordinal
                moris_index tTreatedSideOrdinal = mListOfSideOrdinals( iSideset );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the treated IWG
                    IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                    // FIXME
                    tTreatedIWG->set_nodal_weak_bcs( this->get_weak_bcs() );

                    // get the index of the residual dof type for the ith IWG in the list of element dof type
                    uint tIWGResDofIndex
                        = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                    Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                    uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = this->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                    // create an integration rule for the ith IWG
                    mtk::Geometry_Type tSideGeometryType = mElementBlock->get_block_geometry_interpolator()->get_side_geometry_type();
                    Integration_Rule tIntegrationRule( tSideGeometryType,
                                                       Integration_Type::GAUSS,
                                                       this->get_auto_integration_order( tSideGeometryType ),
                                                       Integration_Type::GAUSS,
                                                       Integration_Order::BAR_1 );

                    // create an integrator for the ith IWG
                    Integrator tIntegrator( tIntegrationRule );

                    //get number of integration points
                    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

                    // get integration points
                    Matrix< DDRMat > tSurfRefIntegPoints = tIntegrator.get_points();

                   // get integration weights
                   Matrix< DDRMat > tIntegWeights = tIntegrator.get_weights();

                   for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                   {
                       // get integration point location in the reference surface
                       Matrix< DDRMat > tSurfRefIntegPointI = tSurfRefIntegPoints.get_column( iGP );

                       // get integration point location in the reference volume
                       Matrix< DDRMat > tVolRefIntegPointI
                           = mElementBlock->get_block_geometry_interpolator()->surf_val( tSurfRefIntegPointI,
                                                              tTreatedSideOrdinal );

                       // set integration point
                       for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                       {
                           tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                       }

                       // compute integration point weight x detJ
                       real tSurfDetJ;
                       Matrix< DDRMat > tNormal;
                       mElementBlock->get_block_geometry_interpolator()->surf_det_J( tSurfDetJ,
                                                          tNormal,
                                                          tSurfRefIntegPointI,
                                                          tTreatedSideOrdinal );
                       tTreatedIWG->set_normal( tNormal );

                       real tWStar = tIntegWeights( iGP ) * tSurfDetJ;

                       // compute jacobian at evaluation point
                       Cell< Matrix< DDRMat > > tJacobians;
                       tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                       // add contribution to jacobian from evaluation point
                       for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                       {
                           uint tIWGActiveDofIndex
                               = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                           uint tJacIndex
                               = tIWGResDofIndex * mElementBlock->get_num_interpolators() + tIWGActiveDofIndex;

                           mJacobianElement( tJacIndex )
                               = mJacobianElement( tJacIndex ) + tWStar * tJacobians( iIWGFI );
                       }
                   }
                }

            }

            // jacobian assembly
            uint tCounterI = 0;
            uint tCounterJ = 0;
            uint startI, stopI, startJ, stopJ;

            for ( uint i = 0; i < mElementBlock->get_num_interpolators(); i++ )
            {
                startI = tCounterI;
                stopI  = tCounterI + mElementBlock->get_block_field_interpolator()( i )->get_number_of_space_time_coefficients() - 1;

                tCounterJ = 0;
                for ( uint j = 0; j < mElementBlock->get_num_interpolators(); j++ )
                {
                    startJ = tCounterJ;
                    stopJ  = tCounterJ + mElementBlock->get_block_field_interpolator()( j )->get_number_of_space_time_coefficients() - 1;

                    mJacobian({ startI, stopI },{ startJ, stopJ }) = mJacobianElement( i * mElementBlock->get_num_interpolators() + j ).matrix_data();

                    tCounterJ = stopJ + 1;
                }
                tCounterI = stopI + 1;
            }
            //print( mJacobian,"mJacobian" );
        }

//------------------------------------------------------------------------------

        void Element_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
