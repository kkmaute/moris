#include <iostream>
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"           //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Time_Sideset::Element_Time_Sideset( mtk::Cell                 * aCell,
                                                    moris::Cell< IWG* >       & aIWGs,
                                                    moris::Cell< Node_Base* > & aNodes )
                                                  : Element( aCell, aIWGs, aNodes )
        {
            //create a geometry interpolation rule
            //FIXME: forced interpolation type and order
            Interpolation_Rule tGeometryInterpolationRule( mCell->get_geometry_type(),
                                                           Interpolation_Type::LAGRANGE,
                                                           this->get_auto_interpolation_order(),
                                                           Interpolation_Type::LAGRANGE,
                                                           mtk::Interpolation_Order::LINEAR );

            // create a geometry intepolator
            mGeometryInterpolator = new Geometry_Interpolator( tGeometryInterpolationRule );

            // set the geometry interpolator coefficients xHat and tHat
            //FIXME: tHat are set by default but should come from solver
            Matrix< DDRMat > tTHat( 2, 1 ); tTHat( 0 ) = 0.0; tTHat( 1 ) = 1.0;
            mGeometryInterpolator->set_coeff( mCell->get_vertex_coords(), tTHat );

            // create field interpolators for the element
            mFieldInterpolators = this->create_element_field_interpolators( mGeometryInterpolator );
        }

//------------------------------------------------------------------------------

        Element_Time_Sideset::~Element_Time_Sideset()
        {
            // delete the geometry interpolator pointer
            if ( mGeometryInterpolator != NULL )
            {
                delete mGeometryInterpolator;
            }

            // delete the field interpolator pointers
            for ( uint i = 0; i < mNumOfInterp; i++ )
            {
                if ( mFieldInterpolators( i ) != NULL )
                {
                    delete mFieldInterpolators( i );
                }
            }
        }

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_residual()
        {
            // FIXME: get the ordinal of the faces
            moris::Cell< uint > tTimeSidesetOrdinals = { 0 };
            uint tNumOfSideSets = tTimeSidesetOrdinals.size();

            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_element_field_interpolators_coefficients( mFieldInterpolators );

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get time sideset coordinates in the parent parametric space
                Matrix< DDRMat > tTimeSidesetParamCoords
                    = mGeometryInterpolator->get_time_sideset_param_coords( tTimeSidesetOrdinals( iSideset ) );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the index of the residual dof type for the ith IWG in the list of element dof type
                    uint tIWGResDofIndex
                        = mInterpDofTypeMap( static_cast< int >( mIWGs( iIWG )->get_residual_dof_type()( 0 ) ) );

                    Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = mIWGs( iIWG )->get_active_dof_types();
                    uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = this->get_IWG_field_interpolators( mIWGs( iIWG ), mFieldInterpolators );

                    // create an integration rule for the ith IWG
                    //FIXME: set by default
                    Integration_Rule tIntegrationRule( mCell->get_geometry_type(),
                                                       Integration_Type::GAUSS,
                                                       this->get_auto_integration_order(),
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
                       print(tSurfRefIntegPointI,"tSurfRefIntegPointI");

                       // get integration point location in the reference volume
                       Matrix< DDRMat > tVolRefIntegPointI
                           = mGeometryInterpolator->NXi( tSurfRefIntegPointI ) * tTimeSidesetParamCoords;
                       tVolRefIntegPointI = trans( tVolRefIntegPointI );
                       print(tVolRefIntegPointI,"tVolRefIntegPointI");

                       // set integration point
                       for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                       {
                           tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                       }

                       // compute determinant of Jacobian mapping
                       real tSurfDetJ = mGeometryInterpolator->surf_det_J( tSurfRefIntegPointI );

                       // compute integration point weight x detJ
                       real tWStar = tIntegWeights( iGP ) * tSurfDetJ;

                       // compute jacobian at evaluation point
                       Matrix< DDRMat > tResidual;
                       mIWGs( iIWG )->compute_residual( tResidual, tIWGInterpolators );

                       // add contribution to jacobian from evaluation point
                       mResidualElement( tIWGResDofIndex )
                           = mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
                   }
                }
            }
        }

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian()
        {
            //MORIS_ERROR( false, " Element_Time_Sideset::compute_jacobian - not implemented. ");

            // FIXME: get the ordinal of the faces
            moris::Cell< uint > tTimeSidesetOrdinals = { 0 };
            uint tNumOfSideSets = tTimeSidesetOrdinals.size();

            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_element_field_interpolators_coefficients( mFieldInterpolators );

            // loop over the sideset faces
            for ( uint iSideset = 0; iSideset < tNumOfSideSets; iSideset++ )
            {
                // get time sideset coordinates in the parent parametric space
                Matrix< DDRMat > tTimeSidesetParamCoords
                    = mGeometryInterpolator->get_time_sideset_param_coords( tTimeSidesetOrdinals( iSideset ) );

                // loop over the IWGs
                for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
                {
                    // get the index of the residual dof type for the ith IWG in the list of element dof type
                    //uint tIWGResDofIndex
                    //    = mInterpDofTypeMap( static_cast< int >( mIWGs( iIWG )->get_residual_dof_type()( 0 ) ) );

                    Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = mIWGs( iIWG )->get_active_dof_types();
                    uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                    // get the field interpolators for the ith IWG in the list of element dof type
                    Cell< Field_Interpolator* > tIWGInterpolators
                        = this->get_IWG_field_interpolators( mIWGs( iIWG ), mFieldInterpolators );

                    // create an integration rule for the ith IWG
                    //FIXME: set by default
                    Integration_Rule tIntegrationRule( mCell->get_geometry_type(),
                                                       Integration_Type::GAUSS,
                                                       this->get_auto_integration_order(),
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
                       print(tSurfRefIntegPointI,"tSurfRefIntegPointI");

                       // get integration point location in the reference volume
                       Matrix< DDRMat > tVolRefIntegPointI
                           = mGeometryInterpolator->NXi( tSurfRefIntegPointI ) * trans( tTimeSidesetParamCoords );
                       tVolRefIntegPointI = trans( tVolRefIntegPointI );
                       print(tVolRefIntegPointI,"tVolRefIntegPointI");
                       tVolRefIntegPointI( 3 ) = -1.0;

                       // set integration point
                       for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                       {
                           tIWGInterpolators( iIWGFI )->set_space_time( tVolRefIntegPointI );
                       }

                       // compute determinant of Jacobian mapping
                       real tSurfDetJ = mGeometryInterpolator->surf_det_J( tSurfRefIntegPointI );

                       // compute integration point weight x detJ
                       real tWStar = tIntegWeights( iGP ) * tSurfDetJ;
                       std::cout<<tWStar<<std::endl;
                       std::cout<<"----------------"<<std::endl;

//                       // compute jacobian at evaluation point
//                       Matrix< DDRMat > tResidual;
//                       mIWGs( iIWG )->compute_residual( tResidual, tIWGInterpolators );
//
//                       // add contribution to jacobian from evaluation point
//                       mResidualElement( tIWGResDofIndex )
//                           = mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
                   }
                }
            }
        }

//------------------------------------------------------------------------------

        void Element_Time_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element_Time_Sideset::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

        moris::Cell< fem::Field_Interpolator* >
        Element_Time_Sideset::create_element_field_interpolators
        ( fem::Geometry_Interpolator* aGeometryInterpolator )
        {
            // cell of field interpolators
            Cell< Field_Interpolator* > tFieldInterpolators( mNumOfInterp, nullptr );

            // loop on the dof type groups and create a field interpolator for each
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the ith dof type group
                Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                // create the field interpolation rule for the ith dof type group
                //FIXME: space interpolation based on the mtk::Cell
                //FIXME: time  interpolation set to constant
                Interpolation_Rule tFieldInterpolationRule( mCell->get_geometry_type(),
                                                            Interpolation_Type::LAGRANGE,
                                                            this->get_auto_interpolation_order(),
                                                            Interpolation_Type::LAGRANGE,
                                                            mtk::Interpolation_Order::LINEAR );

                // get number of field interpolated by the ith field interpolator
                uint tNumOfFields = tDofTypeGroup.size();

                // create an interpolator for the ith dof type group
                tFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                   tFieldInterpolationRule,
                                                                   aGeometryInterpolator );
            }
            return tFieldInterpolators;
        }

//------------------------------------------------------------------------------

        void
        Element_Time_Sideset::set_element_field_interpolators_coefficients
            ( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // loop on the dof types
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the ith dof type group
                Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                //FIXME:forced coefficients
                // get the pdof values for the ith dof type group
                Matrix< DDRMat > tCoeff(16, 1, 0.0);
                //this->get_my_pdof_values( tDofTypeGroup, tCoeff );
                print(tCoeff,"tCoeff");

                // set the field coefficients
                aFieldInterpolators( i )->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void
		Element_Time_Sideset::initialize_mJacobianElement_and_mResidualElement
        ( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            mJacobianElement.resize( mNumOfInterp * mNumOfInterp );
            mResidualElement.resize( mNumOfInterp );

            uint tTotalDof = 0;
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get number of pdofs for the ith dof type
                uint tNumOfDofi = aFieldInterpolators( i )->get_number_of_space_time_coefficients();

                // get total number of dof
                tTotalDof = tTotalDof + tNumOfDofi;

                // set mResidualElement size
                mResidualElement( i ).set_size( tNumOfDofi, 1, 0.0 );

                for( uint j = 0; j < mNumOfInterp; j++ )
                {
                    // get number of pdofs for the ith dof type
                    uint tNumOfDofj = aFieldInterpolators( j )->get_number_of_space_time_coefficients();

                    // set mResidualElement size
                    mJacobianElement( i * mNumOfInterp + j ).set_size( tNumOfDofi, tNumOfDofj, 0.0 );
                }
            }

            mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );
            mResidual.set_size( tTotalDof, 1, 0.0 );
        }


//------------------------------------------------------------------------------
        moris::Cell< Field_Interpolator* >
        Element_Time_Sideset::get_IWG_field_interpolators( IWG*                               & aIWG,
                                                      moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // ask the IWG for its active dof types
            Cell< Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();

            // number of active dof type for the IWG
            uint tNumOfIWGActiveDof = tIWGActiveDof.size();

            // select associated active interpolators
            Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
            for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGDofIndex = mInterpDofTypeMap( static_cast< int >( tIWGActiveDof( i )( 0 ) ) );

                // select the corresponding interpolator
                tIWGFieldInterpolators( i ) = aFieldInterpolators( tIWGDofIndex );
            }
            return tIWGFieldInterpolators;
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
