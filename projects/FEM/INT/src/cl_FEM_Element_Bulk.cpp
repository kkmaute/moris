#include <iostream>

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"   //FEM/INT/src
#include "cl_FEM_Element_Block.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Bulk::Element_Bulk( mtk::Cell    const * aCell,
                                    Element_Block      * aElementBlock,
                                    Cluster            * aCluster ) : Element( aCell, aElementBlock, aCluster )
        {

        }

//------------------------------------------------------------------------------

        Element_Bulk::~Element_Bulk(){}

//------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian()
        {
            // initialize mJacobianElement and mResidualElement
            mCluster->initialize_mJacobianElement_and_mResidualElement();

            // get pdofs values for the element
//            mCluster->get_my_pdof_values();                    // FIXME

            // set the field interpolators coefficients
            mCluster->set_field_interpolators_coefficients();

            // set the geometry interpolator coefficients
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mCluster->mTime );         // FIXME

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                // FIXME
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG
                // in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators = this->get_IWG_field_interpolators( tTreatedIWG,
                                                                                                   mElementBlock->get_block_field_interpolator() );

                // create an integration rule for the ith IWG
                //FIXME: set by default
                Integration_Rule tIntegrationRule( mCell->get_geometry_type(),
                                                   Integration_Type::GAUSS,
                                                   this->get_auto_integration_order( mCell->get_geometry_type() ),
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

                // create an integrator for the ith IWG
                Integrator tIntegrator( tIntegrationRule );

                // get number of integration points
                uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

                // get integration points
                Matrix< DDRMat > tIntegPoints = tIntegrator.get_points();

                // get integration weights
                Matrix< DDRMat > tIntegWeights = tIntegrator.get_weights();

                // loop over integration points
                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the iGP integration point
                    Matrix< DDRMat > tTreatedIntegPoint = tIntegPoints.get_column( iGP );

                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tTreatedIntegPoint );
                    }

                    // compute integration point weight x detJ
                    real tWStar = tIntegWeights( iGP ) * mElementBlock->get_block_geometry_interpolator()->det_J( tTreatedIntegPoint );

                    // compute jacobian at evaluation point
                    moris::Cell< Matrix< DDRMat > > tJacobians;

                    tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                    {
                        uint tIWGActiveDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                        uint tJacIndex = tIWGResDofIndex * mElementBlock->get_num_interpolators() + tIWGActiveDofIndex;

                        mCluster->mJacobianElement( tJacIndex ) = mCluster->mJacobianElement( tJacIndex ) + tWStar * tJacobians( iIWGFI );
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

                    mCluster->mJacobian({ startI, stopI },{ startJ, stopJ }) = mCluster->mJacobianElement( i * mElementBlock->get_num_interpolators() + j ).matrix_data();

                    tCounterJ = stopJ + 1;
                }
                tCounterI = stopI + 1;
            }
            // print jacobian for check
            print( mCluster->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_residual()
        {
            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement();

            // get pdofs values for the element
//            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_field_interpolators_coefficients();

            // set the geometry interpolator coefficients
            mElementBlock->get_block_geometry_interpolator()->set_coeff( mCell->get_vertex_coords(), mCluster->mTime );

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                // FIXME: enforced nodal weak bcs
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG
                // in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = this->get_IWG_field_interpolators( tTreatedIWG,
                            mElementBlock->get_block_field_interpolator() );

                // create an integration rule for the ith IWG
                //FIXME: set by default
                Integration_Rule tIntegrationRule( mCell->get_geometry_type(),
                                                   Integration_Type::GAUSS,
                                                   this->get_auto_integration_order( mCell->get_geometry_type() ),
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

                // create an integrator for the ith IWG
                Integrator tIntegrator( tIntegrationRule );

                //get number of integration points
                uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

                // get integration points
                Matrix< DDRMat > tIntegPoints = tIntegrator.get_points();

                // get integration weights
                Matrix< DDRMat > tIntegWeights = tIntegrator.get_weights();

                // loop over integration points
                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the kth integration point
                    Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tIntegPointI );
                    }

                    // compute integration point weight x detJ
                    real tWStar = mElementBlock->get_block_geometry_interpolator()->det_J( tIntegPointI ) * tIntegWeights( iGP );

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    mElementBlock->get_IWGs()( iIWG )->compute_residual( tResidual, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    mCluster->mResidualElement( tIWGResDofIndex ) = mCluster->mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
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
                mCluster->mResidual( { startI, stopI }, { 0 , 0 } ) = mCluster->mResidualElement( iBuild ).matrix_data();

                // update the row counter
                tCounterI = stopI + 1;
            }
//            // print residual for check
//            print( mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

//        real Element::compute_integration_error(
//                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
//        {
//            // create field interpolation rule
//            Interpolation_Rule tFieldInterpolationRule( this->get_geometry_type(),
//                                                        Interpolation_Type::LAGRANGE,
//                                                        this->get_interpolation_order() );
//            // <- add second type in order
//            //    to interpolate in space and time
//
//            // create geometry interpolation rule
//            Interpolation_Rule tGeometryInterpolationRule( this->get_geometry_type(),
//                                                           Interpolation_Type::LAGRANGE,
//                                                           this->get_interpolation_order() );
//
//            // create integration rule
//            Integration_Rule tIntegration_Rule( this->get_geometry_type(),
//                                                Integration_Type::GAUSS,
//                                                this->get_auto_integration_order() );
//
//            // set number of fields
//            uint tNumberOfFields = 1;
//
//            // create interpolator
//            Interpolator tInterpolator( this->get_node_coords(),
//                                        tNumberOfFields,
//                                        tFieldInterpolationRule,
//                                        tGeometryInterpolationRule,
//                                        tIntegration_Rule );
//
//            // get number of points
//            auto tNumberOfIntegrationPoints
//                = tInterpolator.get_number_of_integration_points();
//
//            real aError = 0.0;
//
//            mIWG->create_matrices( &tInterpolator );
//
//            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
//            {
//                 //evaluate shape function at given integration point
//                aError += mIWG->compute_integration_error( mPdofValues,
//                                                           aFunction,
//                                                           k )
//                        * tInterpolator.get_det_J( k )
//                        * tInterpolator.get_integration_weight( k );
//            }
//
//            //std::cout << "Element error " << aError << std::endl;
//            mIWG->delete_matrices();
//
//            return aError;
//        }

//------------------------------------------------------------------------------

//        real Element::compute_element_average_of_scalar_field()
//        {
//
//            // create field interpolation rule
//            Interpolation_Rule tFieldInterpolationRule( this->get_geometry_type(),
//                                                        Interpolation_Type::LAGRANGE,
//                                                        this->get_interpolation_order() );
//            // <- add second type in order
//            //    to interpolate in space and time
//
//            // create geometry interpolation rule
//            Interpolation_Rule tGeometryInterpolationRule( this->get_geometry_type(),
//                                                           Interpolation_Type::LAGRANGE,
//                                                           mtk::Interpolation_Order::LINEAR );
//
//            // create integration rule
//            Integration_Rule tIntegration_Rule( this->get_geometry_type(),
//                                                Integration_Type::GAUSS,
//                                                this->get_auto_integration_order() );
//
//            // set number of fields
//            uint tNumberOfFields = 1;
//
//            // create interpolator
//            Interpolator tInterpolator( this->get_node_coords(),
//                                        tNumberOfFields,
//                                        tFieldInterpolationRule,
//                                        tGeometryInterpolationRule,
//                                        tIntegration_Rule );
//
//            // get number of points
//            auto tNumberOfIntegrationPoints
//                = tInterpolator.get_number_of_integration_points();
//
//            mIWG->create_matrices( &tInterpolator );
//
//            real aValue  = 0.0;Cell< Field_Interpolator* > tFieldInterpolators
//            real tWeight = 0.0;
//
//            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
//            {
//                real tScale = tInterpolator.get_integration_weight( k )
//                            * tInterpolator.get_det_J( k );
//
//                aValue += mIWG->interpolate_scalar_at_point( mNodalWeakBCs, k )
//                        * tScale;
//
//                tWeight += tScale;
//
//            }
//
//            // close IWG object
//            mIWG->delete_matrices();
//
//            return aValue / tWeight;
//
//        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
