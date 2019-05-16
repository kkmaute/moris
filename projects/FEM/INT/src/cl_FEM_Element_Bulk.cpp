#include <iostream>

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"   //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Bulk::Element_Bulk( mtk::Cell    const * aCell,
                                    Set                * aElementBlock,
                                    Cluster            * aCluster ) : Element( aCell, aElementBlock, aCluster )
        {
        }

//------------------------------------------------------------------------------

        Element_Bulk::~Element_Bulk(){}

//------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mElementBlock->get_block_IG_geometry_interpolator()->set_space_coeff( mCell->get_vertex_coords());
            mElementBlock->get_block_IG_geometry_interpolator()->set_time_coeff ( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            //FIXME for global param coords from mesh
            mElementBlock->get_block_IG_geometry_interpolator()->set_param_coeff();

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                // FIXME set nodal weak BCs
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG
                // in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = mElementBlock->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                // get number of integration points
                uint tNumOfIntegPoints = mElementBlock->get_num_integration_points();

                // loop over integration points
                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the ith integration point in the IG param space
                    Matrix< DDRMat > tLocalIntegPoint = mElementBlock->get_integration_points().get_column( iGP );

                    // bring the ith integration point in the IP param space
                    Matrix< DDRMat > tGlobalIntegPoint = mElementBlock->get_block_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint );

                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                    }

                    // compute integration point weight
                    real tWStar = mElementBlock->get_integration_weights()( iGP )
                                * mElementBlock->get_block_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                    // compute jacobian at evaluation point
                    moris::Cell< Matrix< DDRMat > > tJacobians;
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
                        mElementBlock->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                            = mElementBlock->mJacobian( { startIDof, stopIDof }, { startJDof, stopJDof } )
                            + tWStar * tJacobians( iIWGFI );
                    }
                }
            }

//            // jacobian assembly
//            uint tCounterI = 0;
//            uint tCounterJ = 0;
//            uint startI, stopI, startJ, stopJ;
//
//            for ( uint i = 0; i < mElementBlock->get_num_interpolators(); i++ )
//            {
//                startI = tCounterI;
//                stopI  = tCounterI + mElementBlock->get_block_field_interpolator()( i )->get_number_of_space_time_coefficients() - 1;
//
//                tCounterJ = 0;
//                for ( uint j = 0; j < mElementBlock->get_num_interpolators(); j++ )
//                {
//                    startJ = tCounterJ;
//                    stopJ  = tCounterJ + mElementBlock->get_block_field_interpolator()( j )->get_number_of_space_time_coefficients() - 1;
//
//                    mCluster->mJacobian({ startI, stopI },{ startJ, stopJ }) = mCluster->mJacobian({ startI, stopI },{ startJ, stopJ }) +
//                                                                               mCluster->mJacobianElement( i * mElementBlock->get_num_interpolators() + j ).matrix_data();
//
//                    tCounterJ = stopJ + 1;
//                }
//                tCounterI = stopI + 1;
//            }
//            // print jacobian for check
//            print( mCluster->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_residual()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mElementBlock->get_block_IG_geometry_interpolator()->set_space_coeff( mCell->get_vertex_coords());
            mElementBlock->get_block_IG_geometry_interpolator()->set_time_coeff(  mCluster->mTime );

            // set the geometry interpolator global param space and time coefficients for integration cell
            //FIXME for global param coords from mesh
            mElementBlock->get_block_IG_geometry_interpolator()->set_param_coeff();

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mElementBlock->get_IWGs()( iIWG );

                // FIXME: enforced nodal weak bcs
                tTreatedIWG->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex = mInterpDofTypeMap( static_cast< int >( tTreatedIWG->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = tTreatedIWG->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = mElementBlock->get_IWG_field_interpolators( tTreatedIWG, mElementBlock->get_block_field_interpolator() );

                //get number of integration points
                uint tNumOfIntegPoints = mElementBlock->get_num_integration_points();

                // loop over integration points
                for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the ith integration point in the IG param space
                    Matrix< DDRMat > tLocalIntegPoint = mElementBlock->get_integration_points().get_column( iGP );

                    // bring the ith integration point in the IP param space
                    Matrix< DDRMat > tGlobalIntegPoint = mElementBlock->get_block_IG_geometry_interpolator()->map_integration_point( tLocalIntegPoint );

                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tGlobalIntegPoint );
                    }

                    // compute integration point weight
                    real tWStar = mElementBlock->get_integration_weights()( iGP )
                                * mElementBlock->get_block_IG_geometry_interpolator()->det_J( tLocalIntegPoint );

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    mElementBlock->get_IWGs()( iIWG )->compute_residual( tResidual, tIWGInterpolators );

                    // get location of computed residual in global element residual
                    uint startDof = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 0 );
                    uint stopDof  = mElementBlock->get_interpolator_dof_assembly_map()( tIWGResDofIndex, 1 );

                    // add contribution to residual from evaluation point
                    mElementBlock->mResidual( { startDof, stopDof }, { 0, 0 } )
                              = mElementBlock->mResidual( { startDof, stopDof }, { 0, 0 } ) + tResidual * tWStar;
                }
            }

//            // residual assembly
//            uint tCounterI = 0;
//            uint startI, stopI;
//
//            // loop over the field interpolators
//            for ( uint iBuild = 0; iBuild < mElementBlock->get_num_interpolators(); iBuild++ )
//            {
//                // get the row position in the residual matrix
//                startI = tCounterI;
//                stopI  = tCounterI + mElementBlock->get_block_field_interpolator()( iBuild )->get_number_of_space_time_coefficients() - 1;
//
//                // fill the global residual
//                mCluster->mResidual( { startI, stopI }, { 0 , 0 } ) = mCluster->mResidual( { startI, stopI }, { 0 , 0 } ) +
//                                                                      mCluster->mResidualElement( iBuild ).matrix_data();
//
//                // update the row counter
//                tCounterI = stopI + 1;
//            }
//            // print residual for check
//            print( mCluster->mResidual, " mResidual " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Bulk::compute_jacobian_and_residual - not implemented. ");
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
