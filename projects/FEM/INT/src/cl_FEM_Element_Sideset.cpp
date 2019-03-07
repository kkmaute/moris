#include <iostream>
#include "cl_FEM_Element_Sideset.hpp"           //FEM/INT/src

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "fn_det.hpp"
#include "fn_eye.hpp"

#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp"   //MTK/src

#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"            //FEM/INT/src

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Sideset::Element_Sideset( mtk::Cell          * aCell,
                                          Cell< IWG* >       & aIWGs,
                                          Cell< Node_Base* > & aNodes ) : Element( aCell, aIWGs, aNodes )
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
            // end: create a map of the element active dof type list------------------------

            // begin: create a field interpolator for each element active dof type----------
            //create a geometry interpolation rule for the bulk
            //FIXME: set values
            Interpolation_Rule tGeometryInterpolationRuleBulk( mCell->get_geometry_type(),
                                                               Interpolation_Type::LAGRANGE,
                                                               this->get_auto_interpolation_order(),
                                                               Interpolation_Type::LAGRANGE,
                                                               mtk::Interpolation_Order::LINEAR );

            // create a geometry intepolator
            mGeometryInterpolator = new Geometry_Interpolator( tGeometryInterpolationRuleBulk );

            // set the geometry interpolator coefficients xHat and THat
            //FIXME: tHat are set by default but should come from solver
            Matrix< DDRMat > tTHat( 2, 1 ); tTHat( 0 ) = 0.0; tTHat( 1 ) = 1.0;
            mGeometryInterpolator->set_coeff( mCell->get_vertex_coords(), tTHat );

            // create field interpolators for the element
            mFieldInterpolators = this->create_element_field_interpolators( mGeometryInterpolator );
            // end: create a field interpolator for each element active dof type------------

        }

//------------------------------------------------------------------------------

        Element_Sideset::~Element_Sideset()
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

        void Element_Sideset::compute_jacobian()
        {
            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            //FIXME forced values
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
                    // set evaluation point
                    for ( uint iIWGFI = 0; iIWGFI < tNumOfIWGActiveDof; iIWGFI++ )
                    {
                        tIWGInterpolators( iIWGFI )->set_space_time( tIntegPoints.get_column( iGP ) );
                    }

                    // compute Integration point weight x detJ
                    real tWStar = mGeometryInterpolator->det_J( tIntegPoints.get_column( iGP ) )
                                * tIntegWeights( iGP );

                    // compute jacobian at evaluation point
                    Cell< Matrix< DDRMat > > tJacobians( tNumOfIWGActiveDof );
                    mIWGs( iIWG )->compute_jacobian( tJacobians,
                                                     tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    for ( uint l = 0; l < tNumOfIWGActiveDof; l++)
                    {
                        uint tIWGActiveDofIndex
                            = mInterpDofTypeMap( static_cast< int >( tIWGActiveDofType( l )( 0 ) ) );

                        uint tJacIndex
                            = tIWGResDofIndex * mNumOfInterp + tIWGActiveDofIndex;

                        mJacobianElement( tJacIndex )
                            = mJacobianElement( tJacIndex ) + tJacobians( l ) * tWStar;
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
                //std::cout<<"startI:"; std::cout<<startI<<std::endl;
                //std::cout<<"stopI: "; std::cout<<stopI <<std::endl;

                tCounterJ = 0;
                for ( uint j = 0; j < mNumOfInterp; j++ )
                {
                    startJ = tCounterJ;
                    stopJ  = tCounterJ + mFieldInterpolators( j )->get_number_of_space_time_coefficients() - 1;
                    //std::cout<<"startJ:"; std::cout<<startJ<<std::endl;
                    //std::cout<<"stopJ: "; std::cout<<stopJ <<std::endl;

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

        void Element_Sideset::compute_residual()
        {
            // initialize mJacobianElement and mResidualElement
            this->initialize_mJacobianElement_and_mResidualElement( mFieldInterpolators );

            //FIXME: forced values
            // get pdofs values for the element
            this->get_my_pdof_values();

            // set field interpolators coefficients
            this->set_element_field_interpolators_coefficients( mFieldInterpolators );

            // loop over the IWGs
            for( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // FIXME
                mIWGs( i )->set_nodal_weak_bcs( this->get_weak_bcs() );

                // get the index of the residual dof type for the ith IWG
                // in the list of element dof type
                uint tIWGResDofIndex
                    = mInterpDofTypeMap( static_cast< int >( mIWGs( i )->get_residual_dof_type()( 0 ) ) );

                Cell< Cell< MSI::Dof_Type > > tIWGActiveDofType = mIWGs( i )->get_active_dof_types();
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
                for( uint k = 0; k < tNumOfIntegPoints; k++ )
                {
                    // set evaluation point
                    for ( uint l = 0; l < tNumOfIWGActiveDof; l++ )
                    {
                        tIWGInterpolators( l )->set_space_time( tIntegPoints.get_column( k ) );
                    }

                    // compute integration point weight x detJ
                    real tWStar = mGeometryInterpolator->det_J( tIntegPoints.get_column( k ) )
                                * tIntegWeights( k );

                    // compute jacobian at evaluation point
                    Matrix< DDRMat > tResidual;
                    mIWGs( i )->compute_residual( tResidual, tIWGInterpolators );

                    // add contribution to jacobian from evaluation point
                    mResidualElement( tIWGResDofIndex )
                        = mResidualElement( tIWGResDofIndex ) + tResidual * tWStar;
                }
            }

            // residual assembly
            uint tCounterI = 0;
            uint startI, stopI;

            // loop over the field interpolators
            for ( uint i = 0; i < mNumOfInterp; i++ )
            {
                // get the row position in the residual matrix
                startI = tCounterI;
                stopI  = tCounterI + mFieldInterpolators( i )->get_number_of_space_time_coefficients() - 1;

                // fill the global residual
                mResidual( { startI, stopI }, { 0 , 0 } ) = mResidualElement( i ).matrix_data();

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

        void Element_Sideset::compute_jacobian_and_residual()
        {
            MORIS_ERROR( false, " Element::compute_jacobian_and_residual - not implemented. ");
        }

//------------------------------------------------------------------------------

        Cell< Field_Interpolator* >
        Element_Sideset::create_element_field_interpolators
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
        Element_Sideset::set_element_field_interpolators_coefficients
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

                // set the field coefficients
                aFieldInterpolators( i )->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void
        Element_Sideset::initialize_mJacobianElement_and_mResidualElement
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
        Element_Sideset::get_IWG_field_interpolators( IWG*                        & aIWG,
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
