#include <iostream>

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        Element_Bulk::Element_Bulk( mtk::Cell const  * aCell,
                                    Set              * aSet,
                                    Cluster          * aCluster,
                                    moris::moris_index aCellIndexInCluster ) : Element( aCell, aSet, aCluster, aCellIndexInCluster )
        {}

//------------------------------------------------------------------------------

        Element_Bulk::~Element_Bulk(){}

//------------------------------------------------------------------------------

        void Element_Bulk::compute_jacobian()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_IG_geometry_interpolator()->set_time_coeff ( mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} );//fixme

            // loop over integration points
            for( uint iGP = 0; iGP < mSet->get_num_integration_points(); iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // bring the ith integration point in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iInterp = 0; iInterp < mSet->get_num_interpolators(); iInterp++ )
                {
                    mSet->get_field_interpolator()( iInterp )->set_space_time( tGlobalIntegPoint );
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mSet->get_num_IWG(); iIWG++ )
                {
                    // FIXME set nodal weak BCs
                    mSet->get_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // compute jacobian at evaluation point
                    moris::Cell< Matrix< DDRMat > > tJacobians;
                    mSet->get_IWGs()( iIWG )->compute_jacobian( tJacobians, mSet->get_IWG_field_interpolators()( iIWG ) );
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
                    for ( uint iIWGFI = 0; iIWGFI < mSet->get_IWG_num_active_dof()( iIWG ); iIWGFI++)
                    {
                        // add contribution to jacobian from evaluation point
                        mSet->mJacobian( { mSet->get_IWG_dof_assembly_map()( iIWG )( iIWGFI, 0 ), mSet->get_IWG_dof_assembly_map()( iIWG )( iIWGFI, 1 ) },
                                         { mSet->get_IWG_dof_assembly_map()( iIWG )( iIWGFI, 2 ), mSet->get_IWG_dof_assembly_map()( iIWG )( iIWGFI, 3 ) } )
                                       += tWStar * tJacobians( iIWGFI );
                    }
                }
            }
//            // print jacobian for check
//            print( mCluster->mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_residual()
        {
            // set the geometry interpolator physical space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_coeff( mMasterCell->get_vertex_coords());
            mSet->get_IG_geometry_interpolator()->set_time_coeff(  mCluster->mTime );

            // set the geometry interpolator param space and time coefficients for integration cell
            mSet->get_IG_geometry_interpolator()->set_space_param_coeff( mCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( mCellIndexInCluster) );
            mSet->get_IG_geometry_interpolator()->set_time_param_coeff( {{-1.0}, {1.0}} ); //fixme

            // loop over integration points
            for( uint iGP = 0; iGP < mSet->get_num_integration_points(); iGP++ )
            {
                // get the ith integration point in the IG param space
                Matrix< DDRMat > tLocalIntegPoint = mSet->get_integration_points().get_column( iGP );

                // set the ith integration point in the IG param space for IG geometry interpolator
                mSet->get_IG_geometry_interpolator()->set_space_time( tLocalIntegPoint );

                // bring the ith integration point in the IP param space
                Matrix< DDRMat > tGlobalIntegPoint;
                mSet->get_IG_geometry_interpolator()->map_integration_point( tGlobalIntegPoint );

                // set evaluation point for field interpolator
                for ( uint iInterp = 0; iInterp < mSet->get_num_interpolators(); iInterp++ )
                {
                    mSet->get_field_interpolator()( iInterp )->set_space_time( tGlobalIntegPoint );
                }

                // compute integration point weight
                real tWStar = mSet->get_integration_weights()( iGP )
                            * mSet->get_IG_geometry_interpolator()->det_J();

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mSet->get_num_IWG(); iIWG++ )
                {
                    // FIXME: enforced nodal weak bcs
                    mSet->get_IWGs()( iIWG )->set_nodal_weak_bcs( mCluster->get_weak_bcs() );

                    // compute residual at evaluation point
                    Matrix< DDRMat > tResidual;
                    mSet->get_IWGs()( iIWG )->compute_residual( tResidual, mSet->get_IWG_field_interpolators()( iIWG ) );

                    // add contribution to residual from evaluation point
                    mSet->mResidual( { mSet->get_IWG_dof_assembly_map()( iIWG )( 0, 0 ), mSet->get_IWG_dof_assembly_map()( iIWG )( 0, 1 ) },
                                     { 0, 0 } ) += tWStar * tResidual;
                }
            }
//            // print residual for check
//            print( mSet->mResidual, " mResidual " );
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
