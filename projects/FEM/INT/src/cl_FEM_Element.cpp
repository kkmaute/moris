#include <iostream>
#include "cl_FEM_Element.hpp" //FEM/INT/src

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "fn_det.hpp"
#include "fn_sort.hpp"
#include "fn_eye.hpp"

#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp"   //MTK/src

#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"            //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element::Element( mtk::Cell          * aCell,
                          Cell< IWG* >       & aIWGs,
                          Cell< Node_Base* > & aNodes ) : mCell( aCell ),
                                                          mIWGs( aIWGs )
        {
            // begin: restrict mNodeObj to the element nodes
            //------------------------------------------------------------------------------
            // get vertices from cell
            Cell< mtk::Vertex* > tVertices = aCell->get_vertex_pointers();

            // get number of nodes from cell
            uint tNumOfNodes = tVertices.size();

            // assign node object
            mNodeObj.resize( tNumOfNodes, nullptr );

            // fill node objects
            for( uint i = 0; i < tNumOfNodes; i++)
            {
                mNodeObj( i ) = aNodes( tVertices( i )->get_index() );
            }
            // end: restrict mNodeObj to the element nodes
            //------------------------------------------------------------------------------

            // set size of Weak BCs
            mNodalWeakBCs.set_size( tNumOfNodes, 1 );

            // FIXME: Mathias, please comment
            mTimeSteps.set_size( 1, 1, 0 );

            // begin: create an element active dof type list from IWGs
            //------------------------------------------------------------------------------
            // get the number of IWGs
            mNumOfIWGs = mIWGs.size();

            // set the size of the element active dof type list
            mEqnObjDofTypeList.resize( mNumOfIWGs );

            // loop over the IWGs
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the residual dof type of the ith IWG
                mEqnObjDofTypeList( i ) = mIWGs( i )->get_residual_dof_type();
            }

            // use std::unique and std::distance to create a unique list containing all used dof types
            auto last = std::unique( ( mEqnObjDofTypeList.data() ).data(), ( mEqnObjDofTypeList.data() ).data() + mEqnObjDofTypeList.size() );
            auto pos  = std::distance( ( mEqnObjDofTypeList.data() ).data(), last );
            mEqnObjDofTypeList.resize( pos );
            // end: create an element active dof type list from IWGs
            //------------------------------------------------------------------------------

            // begin: create a map of the element active dof type list
            //------------------------------------------------------------------------------
            // set number of unique pdof type of the element
            mNumOfElemDofTypes = mEqnObjDofTypeList.size();

            // get maximal dof type enum number
            sint tMaxDofTypeEnumNumber = 0;

            // loop over all pdof types to get the highest enum index
            for ( uint i = 0; i < mNumOfElemDofTypes; i++ )
            {
                tMaxDofTypeEnumNumber = std::max( tMaxDofTypeEnumNumber, static_cast< int >( mEqnObjDofTypeList( i ) ) );
            }

            // +1 because c++ is 0 based
            tMaxDofTypeEnumNumber = tMaxDofTypeEnumNumber + 1;

            // set size of mapping matrix
            mElemDofTypeMap.set_size( tMaxDofTypeEnumNumber, 1, -1 );

            // loop over all dof types to create the mapping matrix
            for ( uint i = 0; i < mNumOfElemDofTypes; i++ )
            {
                mElemDofTypeMap( static_cast< int >( mEqnObjDofTypeList( i ) ), 0 ) = i;
            }
            // end: create a map of the element active dof type list
            //------------------------------------------------------------------------------

            // begin: create a field interpolator for each element active dof type
            //------------------------------------------------------------------------------
            // get pdofs values for the element
            mPdofValues.set_size( 16, 1, 0.0 );
            //this->get_my_pdof_values();

            //create a geometry interpolation rule
            //FIXME: set values
            Interpolation_Rule tGeometryInterpolationRule( mCell->get_geometry_type(),
                                                           Interpolation_Type::LAGRANGE,
                                                           mtk::Interpolation_Order::LINEAR ,
                                                           Interpolation_Type::LAGRANGE,
                                                           mtk::Interpolation_Order::LINEAR );

            // create a geometry intepolator
            Geometry_Interpolator* tGeometryInterpolator
                = new Geometry_Interpolator( tGeometryInterpolationRule );

            // set the geometry interpolator coefficients xHat and THat
            //FIXME: tHat are set by default but should come from solver
            Matrix< DDRMat > tTHat( 2, 1); tTHat( 0 ) = 0.0; tTHat( 1 ) = 1.0;
            tGeometryInterpolator->set_coeff( mCell->get_vertex_coords(), tTHat );

            // create field interpolators for the element
            mFieldInterpolators = this->create_element_field_interpolators( tGeometryInterpolator );

            // end: create a field interpolator for each element active dof type
            //------------------------------------------------------------------------------

            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            // set the jacobian martrix to identity
            //FIXME not true for space time element
            eye( tNumOfNodes, tNumOfNodes, mJacobian );

        }

//------------------------------------------------------------------------------

        Integration_Order Element::get_auto_integration_order()
        {
            switch( mCell->get_geometry_type() )
            {
                case( mtk::Geometry_Type::QUAD ) :
                {
                     return Integration_Order::QUAD_3x3;
                     break;
                }

                case( mtk::Geometry_Type::HEX ) :
                {
                    return Integration_Order::HEX_3x3x3;
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, " Element::get_auto_integration_order - not defined for this geometry type. ");
                    return Integration_Order::UNDEFINED;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void Element::compute_jacobian()
        {
            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mElemDofTypeMap( static_cast< int >( mIWGs( iIWG )->get_residual_dof_type() ) );

                Cell< MSI::Dof_Type > tIWGActiveDofType = mIWGs( iIWG )->get_active_dof_types();
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
                                                   Integration_Order::BAR_2 );

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
                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tIntegPoints.get_column( iGP ) );
                    }

                    // compute jacobian at evaluation point
                    Cell< Matrix< DDRMat > > tJacobians( tNumOfIWGActiveDof );
                    mIWGs( iIWG )->compute_jacobian( tJacobians,
                                                     tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    for ( uint l = 0; l < tNumOfIWGActiveDof; l++)
                    {
                        uint tIWGActiveDofIndex
                            = mElemDofTypeMap( static_cast< int >( tIWGActiveDofType( l ) ) );

                        uint tJacIndex
                            = tIWGResDofIndex * mNumOfElemDofTypes + tIWGActiveDofIndex;

                        mJacobianElement( tJacIndex )
                            = mJacobianElement( tJacIndex )
                            + tJacobians( l ) * tIWGInterpolators( 0 )->det_J() * tIntegWeights( iGP );
                    }
                }
            }
            for ( uint iPrint = 0; iPrint < mJacobianElement.size(); iPrint++ )
            {
                print( mJacobianElement( iPrint ), " mJacobianElement " );
            }
        }

//------------------------------------------------------------------------------
        void Element::compute_residual()
        {
            // loop over the IWGs
            for( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mElemDofTypeMap( static_cast< int >( mIWGs( i )->get_residual_dof_type() ) );

                Cell< MSI::Dof_Type > tIWGActiveDofType = mIWGs( i )->get_active_dof_types();
                uint tNumOfIWGActiveDof = tIWGActiveDofType.size();

                // get the field interpolators for the ith IWG
                // in the list of element dof type
                Cell< Field_Interpolator* > tIWGInterpolators
                    = this->get_IWG_field_interpolators( mIWGs( i ),
                                                         mFieldInterpolators );

                // create an integration rule for the ith IWG
                //FIXME: set by default
                Integration_Rule tIntegrationRule( mCell->get_geometry_type(),
                                                   Integration_Type::GAUSS,
                                                   this->get_auto_integration_order(),
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_2 );

                // create an integrator for the ith IWG
                Integrator tIntegrator( tIntegrationRule );

                //get number of integration points
                uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

                // get integration points
                Matrix< DDRMat > tIntegPoints = tIntegrator.get_points();

                // get integration weights
                Matrix< DDRMat > tIntegWeights = tIntegrator.get_weights();

                // loop over integration points
                for( uint k = 0; k < tNumOfIntegPoints; k++ )
                {
                    // set evaluation point
                    for ( uint l = 0; l < tNumOfIWGActiveDof; l++ )
                    {
                        tIWGInterpolators( l )->set_space_time( tIntegPoints.get_column( k ) );
                    }

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    mIWGs( i )->compute_residual( tResidual, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    mResidualElement( tIWGResDofIndex )
                        = mResidualElement( tIWGResDofIndex )
                        + tResidual * tIWGInterpolators( 0 )->det_J() * tIntegWeights( k );
                }
                print( mResidualElement( i ), " mResidualElement " );
            }
        }

//------------------------------------------------------------------------------

        void Element::compute_jacobian_and_residual()
        {
            // loop over the IWGs
            for( uint iIWG = 0; iIWG < mNumOfIWGs; iIWG++ )
            {
                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mElemDofTypeMap( static_cast< int >( mIWGs( iIWG )->get_residual_dof_type() ) );

                Cell< MSI::Dof_Type > tIWGActiveDofType = mIWGs( iIWG )->get_active_dof_types();
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
                                                   Integration_Order::BAR_2 );

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
                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tIntegPoints.get_column( iGP ) );
                    }

                    // compute jacobian at evaluation point
                    Cell< Matrix< DDRMat > > tJacobians( tNumOfIWGActiveDof );
                    Matrix< DDRMat > tResidual;
                    mIWGs( iIWG )->compute_jacobian_and_residual( tJacobians,
                                                                  tResidual,
                                                                  tIWGInterpolators );
                    // add contribution to residual from evaluation point
                    mResidualElement( tIWGResDofIndex )
                        = mResidualElement( tIWGResDofIndex )
                        + tResidual * tIWGInterpolators( 0 )->det_J() * tIntegWeights( iGP );

                    // add contribution to jacobian from evaluation point
                    for ( uint l = 0; l < tNumOfIWGActiveDof; l++)
                    {
                        uint tIWGActiveDofIndex
                            = mElemDofTypeMap( static_cast< int >( tIWGActiveDofType( l ) ) );

                        uint tJacIndex
                            = tIWGResDofIndex * mNumOfElemDofTypes + tIWGActiveDofIndex;

                        mJacobianElement( tJacIndex )
                            = mJacobianElement( tJacIndex )
                            + tJacobians( l ) * tIWGInterpolators( 0 )->det_J() * tIntegWeights( iGP );
                    }
                }
                print( mResidualElement( iIWG ), " mResidualElement " );
            }
            for ( uint iPrint = 0; iPrint < mJacobianElement.size(); iPrint++ )
            {
                print( mJacobianElement( iPrint ), " mJacobianElement " );
            }
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

        /*Mat< moris_index >
        Element::get_adof_indices()
        {
            return sort( mCell->get_adof_indices() );
        }*/

//------------------------------------------------------------------------------

        Cell< Field_Interpolator* >
        Element::create_element_field_interpolators( Geometry_Interpolator* aGeometryInterpolator )
        {
            // cell of field interpolators
            Cell< Field_Interpolator* > tFieldInterpolators( mNumOfElemDofTypes, nullptr );

            // loop on the pdof types and create a field interpolator for each
            for( uint i = 0; i < mNumOfElemDofTypes; i++ )
            {
                // create the field interpolation rule for pdof type i
                //FIXME: space interpolation is based on the mtk::Cell
                //       time interpolation is set
                Interpolation_Rule tFieldInterpolationRule( mCell->get_geometry_type(),
                                                            Interpolation_Type::LAGRANGE,
                                                            mtk::Interpolation_Order::LINEAR ,
                                                            Interpolation_Type::LAGRANGE,
                                                            mtk::Interpolation_Order::LINEAR );

                //set the number of fields for pdof type i
                //FIXME: set to one
                Cell< MSI::Dof_Type > tDofType = { mEqnObjDofTypeList( i ) };
                uint tNumOfFields = tDofType.size();

                // create an interpolator for pdof type i
                tFieldInterpolators( i ) = new Field_Interpolator( tNumOfFields,
                                                                   tFieldInterpolationRule,
                                                                   aGeometryInterpolator );
                // get the pdof values for pdof type i
                // FIXME: set by default,
                // but should come from this->get_my_pdof_values();
                // which pdof type, which order
                //Matrix< DDRMat > tCoeff;
                //this->get_my_pdof_values( tDofType, tCoeff );

                uint tNumOfBases = tFieldInterpolators( i )->get_number_of_space_time_bases();
                Matrix< DDRMat > tCoeff( tNumOfBases, tNumOfFields, 0.0 );
                //tCoeff( 0 ) = 1.0; tCoeff( 1 ) = 1.0; tCoeff( 2 ) = 1.0; tCoeff( 3 ) = 1.0;
                //tCoeff( 4 ) = 2.0; tCoeff( 5 ) = 3.0; tCoeff( 6 ) = 2.0; tCoeff( 7 ) = 2.0;

                // set the field coefficients
                tFieldInterpolators( i )->set_coeff( tCoeff );
            }
            return tFieldInterpolators;
        }

//------------------------------------------------------------------------------
        void
        Element::initialize_mJacobianElement_and_mResidualElement( Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            mJacobianElement.resize( mNumOfElemDofTypes*mNumOfElemDofTypes );
            mResidualElement.resize( mNumOfElemDofTypes );

            for( uint i = 0; i < mNumOfElemDofTypes; i++ )
            {
                // get number of pdofs for the ith dof type
                uint tNumOfDofi = aFieldInterpolators( i )->get_number_of_space_time_bases();

                // set mResidualElement size
                mResidualElement( i ).set_size( tNumOfDofi, 1, 0.0 );

                for( uint j = 0; j < mNumOfElemDofTypes; j++ )
                {
                    // get number of pdofs for the ith dof type
                    uint tNumOfDofj = aFieldInterpolators( j )->get_number_of_space_time_bases();

                    // set mResidualElement size
                    mJacobianElement( i * mNumOfElemDofTypes + j ).set_size( tNumOfDofi, tNumOfDofj, 0.0 );
                }
            }
        }


//------------------------------------------------------------------------------
        Cell< Field_Interpolator* >
        Element::get_IWG_field_interpolators( IWG*                        & aIWG,
                                              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // ask the IWG for its active dof types
            Cell< MSI::Dof_Type> tIWGActiveDof = aIWG->get_active_dof_types();

            // number of active dof type for the IWG
            uint tNumOfIWGActiveDof = tIWGActiveDof.size();

            // select associated active interpolators
            Cell< Field_Interpolator* > tIWGFieldInterpolators( tNumOfIWGActiveDof, nullptr );
            for( uint i = 0; i < tNumOfIWGActiveDof; i++ )
            {
                // find the index of active dof type in the list of element dof type
                uint tIWGDofIndex = mElemDofTypeMap( static_cast< int >( tIWGActiveDof( i ) ) );

                // select the corresponding interpolator
                tIWGFieldInterpolators( i ) = aFieldInterpolators( tIWGDofIndex );
            }
            return tIWGFieldInterpolators;
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
