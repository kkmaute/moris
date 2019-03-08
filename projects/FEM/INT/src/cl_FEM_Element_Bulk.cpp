#include <iostream>

#include "cl_FEM_Element_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Bulk::Element_Bulk( mtk::Cell          * aCell,
                                    Cell< IWG* >       & aIWGs,
                                    Cell< Node_Base* > & aNodes )
                                  : Element( aCell, aIWGs, aNodes )
        {
            // begin: create an element active dof type list from IWGs----------------------

            // get the number of IWGs
            mNumOfIWGs = mIWGs.size();

            // set the size of the element active dof type list
            uint tCounter = 0;
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                tCounter = tCounter + mIWGs( i )->get_residual_dof_type().size();
            }
            mEqnObjDofTypeList.resize( tCounter );

            // loop over the IWGs
            tCounter = 0;
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the residual dof type of the ith IWG
                Cell< MSI::Dof_Type > tDofType = mIWGs( i )->get_residual_dof_type();

                for ( uint j = 0; j < tDofType.size(); j++ )
               {
                   // get the residual dof type of the ith IWG
                   mEqnObjDofTypeList( tCounter ) = tDofType( j );
                   tCounter++;
                }
            }

            // use std::unique and std::distance to create a unique list containing all used dof types
            auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(),
                                     ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
            auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
            mEqnObjDofTypeList.resize( pos );

            //------------------------------------------------------------------------------
            // set the size of the element active dof type list
            mInterpDofTypeList.resize( mNumOfIWGs );

            // loop over the IWGs
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the residual dof type of the ith IWG
                mInterpDofTypeList( i ) = mIWGs( i )->get_residual_dof_type();
            }
            // end: create an element active dof type list from IWGs------------------------

            // begin: create a map of the element active dof type list----------------------
//            // set number of unique pdof type of the element
//            mNumOfElemDofTypes = mEqnObjDofTypeList.size();
//
//            // get maximal dof type enum number
//            sint tMaxDofTypeEnumNumber = 0;
//
//            // loop over all pdof types to get the highest enum index
//            for ( uint i = 0; i < mNumOfElemDofTypes; i++ )
//            {
//                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mEqnObjDofTypeList( i ) ) );
//            }
//
//            for ( uint i = 0; i < tNumOfInterp; i++ )
//            {
//                tMaxDofTypeEnumNumber2 = std::max( tMaxDofTypeEnumNumber2, static_cast< int >( mInterpDofTypeList( i )( 0 ) ) );
//            }
//
//            // +1 because c++ is 0 based
//            tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;
//
//            // set size of mapping matrix
//            mElemDofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );
//
//            // loop over all dof types to create the mapping matrix
//            for ( uint i = 0; i < mNumOfElemDofTypes; i++ )
//            {
//                mElemDofTypeMap( static_cast< int >( mEqnObjDofTypeList( i ) ), 0 ) = i;
//            }

            // set number of unique pdof type of the element
            mNumOfInterp = mInterpDofTypeList.size();

            // get maximal dof type enum number
            sint tMaxDofTypeEnumNumber = 0;

            // loop over all pdof types to get the highest enum index
            for ( uint i = 0; i < mNumOfInterp; i++ )
            {
                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mInterpDofTypeList( i )( 0 ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;

            // set size of mapping matrix
            mInterpDofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

            // loop over all dof types to create the mapping matrix
            for ( uint i = 0; i < mNumOfInterp; i++ )
            {
                mInterpDofTypeMap( static_cast< int >( mInterpDofTypeList( i )( 0 ) ), 0 ) = i;
            }
            //print( mInterpDofTypeMap, "mInterpDofTypeMap" );

            // end: create a map of the element active dof type list------------------------

            // begin: create a field interpolator for each element active dof type----------
            //create a geometry interpolation rule
            //FIXME: set values
            Interpolation_Rule tGeometryInterpolationRule( mCell->get_geometry_type(),
                                                           Interpolation_Type::LAGRANGE,
                                                           this->get_auto_interpolation_order(),
                                                           Interpolation_Type::LAGRANGE,
                                                           mtk::Interpolation_Order::LINEAR );

            // create a geometry intepolator
            mGeometryInterpolator = new Geometry_Interpolator( tGeometryInterpolationRule );

            // set the geometry interpolator coefficients xHat and THat
            //FIXME: tHat are set by default but should come from solver
            Matrix< DDRMat > tTHat( 2, 1 ); tTHat( 0 ) = 0.0; tTHat( 1 ) = 1.0;
            mGeometryInterpolator->set_coeff( mCell->get_vertex_coords(), tTHat );

            // create field interpolators for the element
            mFieldInterpolators = this->create_element_field_interpolators( mGeometryInterpolator );
            // end: create a field interpolator for each element active dof type------------
        }

//------------------------------------------------------------------------------

        Element_Bulk::~Element_Bulk()
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

        void Element_Bulk::compute_jacobian()
        {
            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_element_field_interpolators_coefficients( mFieldInterpolators );

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the treated IWG
                IWG* tTreatedIWG = mIWGs( iIWG );

                // FIXME
                tTreatedIWG->set_nodal_weak_bcs( this->get_weak_bcs() );

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
                                                         mFieldInterpolators );

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
                    real tWStar = mGeometryInterpolator->det_J( tTreatedIntegPoint ) * tIntegWeights( iGP );

                    // compute jacobian at evaluation point
                    Cell< Matrix< DDRMat > > tJacobians( tNumOfIWGActiveDof );
                    tTreatedIWG->compute_jacobian( tJacobians, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++)
                    {
                        uint tIWGActiveDofIndex
                            = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( iIWGFI )( 0 ) ) );

                        uint tJacIndex
                            = tIWGResDofIndex * mNumOfInterp + tIWGActiveDofIndex;

                        mJacobianElement( tJacIndex )
                            = mJacobianElement( tJacIndex ) + tWStar * tJacobians( iIWGFI );
                    }
                }
            }

            // jacobian assembly
            uint tCounterI = 0;
            uint tCounterJ = 0;
            uint startI, stopI, startJ, stopJ;

            for ( uint i = 0; i < mNumOfInterp; i++ )
            {
                startI = tCounterI;
                stopI  = tCounterI + mFieldInterpolators( i )->get_number_of_space_time_coefficients() - 1;

                tCounterJ = 0;
                for ( uint j = 0; j < mNumOfInterp; j++ )
                {
                    startJ = tCounterJ;
                    stopJ  = tCounterJ + mFieldInterpolators( j )->get_number_of_space_time_coefficients() - 1;

                    mJacobian({ startI, stopI },{ startJ, stopJ }) = mJacobianElement( i * mNumOfInterp + j ).matrix_data();

                    tCounterJ = stopJ + 1;
                }
                tCounterI = stopI + 1;
            }

//            // print residual for check
//            for ( uint iPrint = 0; iPrint < mNumOfInterp*mNumOfInterp; iPrint++ )
//            {
//                print( mJacobianElement( iPrint ), " mJacobianElement " );
//            }
//            print( mJacobian, " mJacobian " );
        }

//------------------------------------------------------------------------------

        void Element_Bulk::compute_residual()
        {
            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_element_field_interpolators_coefficients( mFieldInterpolators );

            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // FIXME
                mIWGs( iIWG )->set_nodal_weak_bcs( this->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mInterpDofTypeMap( static_cast< int >( mIWGs( iIWG )->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = mIWGs( iIWG )->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG
                // in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = this->get_IWG_field_interpolators( mIWGs( iIWG ),
                                                         mFieldInterpolators );

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
                    real tWStar = mGeometryInterpolator->det_J( tIntegPointI ) * tIntegWeights( iGP );

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    mIWGs( iIWG )->compute_residual( tResidual, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    mResidualElement( tIWGResDofIndex )
                        = mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
                }
            }

            // residual assembly
            uint tCounterI = 0;
            uint startI, stopI;

            // loop over the field interpolators
            for ( uint iBuild = 0; iBuild < mNumOfInterp; iBuild++ )
            {
                // get the row position in the residual matrix
                startI = tCounterI;
                stopI  = tCounterI + mFieldInterpolators( iBuild )->get_number_of_space_time_coefficients() - 1;

                // fill the global residual
                mResidual( { startI, stopI }, { 0 , 0 } ) = mResidualElement( iBuild ).matrix_data();

                // update the row counter
                tCounterI = stopI + 1;
            }

//            // print residual for check
//            for ( uint iPrint = 0; iPrint < mNumOfInterp; iPrint++ )
//            {
//                print( mResidualElement( iPrint ), " mResidualElement " );
//            }
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

        Cell< Field_Interpolator* >
        Element_Bulk::create_element_field_interpolators
        ( Geometry_Interpolator* aGeometryInterpolator )
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
                                                            Interpolation_Type::CONSTANT,
                                                            mtk::Interpolation_Order::CONSTANT );

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
        Element_Bulk::set_element_field_interpolators_coefficients
        ( Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // loop on the dof types
            for( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the ith dof type group
                Cell< MSI::Dof_Type > tDofTypeGroup = mInterpDofTypeList( i );

                //FIXME:forced coefficients
                // get the pdof values for the ith dof type group
                Matrix< DDRMat > tCoeff;
                this->get_my_pdof_values( tDofTypeGroup, tCoeff );
//                Matrix< DDRMat > tCoeff ( aFieldInterpolators( i )->get_number_of_space_time_bases(),
//                                          aFieldInterpolators( i )->get_number_of_fields(),
//                                          0.0 );
//                tCoeff( 0, 0 ) = 1.0; tCoeff( 1, 0 ) = 2.0; tCoeff( 2, 0 ) = 3.0; tCoeff( 3, 0 ) = 4.0;

                // set the field coefficients
                aFieldInterpolators( i )->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void
        Element_Bulk::initialize_mJacobianElement_and_mResidualElement
        ( Cell< Field_Interpolator* > & aFieldInterpolators )
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
        Cell< Field_Interpolator* >
        Element_Bulk::get_IWG_field_interpolators( IWG*                        & aIWG,
                                                   Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // ask the IWG for its active dof types
            Cell< Cell< MSI::Dof_Type > > tIWGActiveDof = aIWG->get_active_dof_types();

            // number of active dof type for the IWG
            uint tNumOfIWGActiveDof = tIWGActiveDof.size();

            // select associated active interpolators
            Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
            for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
            {
//                // find the index of active dof type in the list of element dof type
//                uint tIWGDofIndex = mElemDofTypeMap( static_cast< int >( tIWGActiveDof( i ) ) );

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
