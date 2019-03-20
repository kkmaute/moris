/*
 * cl_FEM_Element.hpp
 *
 *  Created on: Mar 07, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_HPP_

#include "assert.h"
#include <cmath>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Cell.hpp"                  //MTK/src

#include "cl_MSI_Equation_Object.hpp"       //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Integrator.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Element : public MSI::Equation_Object
    {

    protected:

        //! pointer to cell on mesh
        const mtk::Cell * mCell;

        // cell of pointers to IWG objects
        moris::Cell< IWG* > mIWGs;

        //! node indices of this element
        //  @node: MTK interface returns copy of vertices. T
        //         storing the indices in private matrix is faster,
        //         but might need more memory
        moris::Matrix< IndexMat > mNodeIndices;

        // map of the element active dof types
        moris::Cell< moris::Cell< DDRMat > > mElemDofTypeList;
        moris::Matrix< DDSMat >              mElemDofTypeMap;
        uint                                 mNumOfElemDofTypes;
        uint                                 mNumOfIWGs;

        Geometry_Interpolator*               mGeometryInterpolator;

        moris::Cell< Field_Interpolator* >   mFieldInterpolators;
        moris::Cell< Cell< MSI::Dof_Type > > mInterpDofTypeList;
        moris::Matrix< DDSMat >              mInterpDofTypeMap;
        uint                                 mNumOfInterp;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         */
        Element( mtk::Cell                 * aCell,
                 moris::Cell< IWG* >       & aIWGs,
                 moris::Cell< Node_Base* > & aNodes )
        {

            // fill the bulk mtk::Cell pointer
            mCell = aCell;

            // fill the cell of IWGs pointers
            mIWGs = aIWGs;

            // select the element nodes from aNodes and fill mNodeObj
            // get vertices from cell
            moris::Cell< mtk::Vertex* > tVertices = aCell->get_vertex_pointers();

            // get number of nodes from cell
            uint tNumOfNodes = tVertices.size();

            // assign node object
            mNodeObj.resize( tNumOfNodes, nullptr );

            // fill node objects
            for( uint i = 0; i < tNumOfNodes; i++)
            {
                mNodeObj( i ) = aNodes( tVertices( i )->get_index() );
            }

            // set size of Weak BCs
            mNodalWeakBCs.set_size( tNumOfNodes, 1 );

            // FIXME: Mathias, please comment
            mTimeSteps.set_size( 1, 1, 1 );

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
            // end: create a map of the element active dof type list------------------------


        };
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        virtual ~Element(){};

//------------------------------------------------------------------------------

        virtual void compute_jacobian() = 0;

//------------------------------------------------------------------------------

        virtual void compute_residual() = 0;

//------------------------------------------------------------------------------

        virtual void compute_jacobian_and_residual() = 0;

//------------------------------------------------------------------------------

//        real compute_integration_error( real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

//        real compute_element_average_of_scalar_field();

//------------------------------------------------------------------------------

        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                           moris::Cell< MSI::Dof_Type > aDofType )
        {
            // get pdofs values for the element
            this->get_my_pdof_values();

            // get a specific dof type profs values
            Matrix< DDRMat > tPdofValues;
            this->get_my_pdof_values( aDofType, tPdofValues );

            // select the required nodal value
            Matrix< IndexMat > tElemVerticesIndices = mCell->get_vertex_inds();
            uint tElemNumOfVertices = mCell->get_number_of_vertices();

            moris_index tVertexIndex;
            for( uint i = 0; i < tElemNumOfVertices; i++ )
            {
                if ( tElemVerticesIndices( i ) == aVertexIndex )
                {
                    tVertexIndex =  i ;
                    break;
                }
            }
            return tPdofValues( tVertexIndex );

        }

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
        /**
         * compute element volume
         */
        real compute_element_volume( Geometry_Interpolator* aGeometryInterpolator )
        {
            //FIXME: enforced Intergation_Type and Integration_Order
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

            // init volume
            real tVolume = 0;

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // compute integration point weight x detJ
                real tWStar = aGeometryInterpolator->det_J( tIntegPoints.get_column( iGP ) )
                            * tIntegWeights( iGP );

                // add contribution to jacobian from evaluation point
                //FIXME: include a thickness if 2D
                tVolume = tVolume + tWStar;
            }

            // FIXME: compute the element size + switch 1D, 2D, 3D
            //real he = std::pow( 6*tVolume/M_PI, 1.0/3.0 );
            //real he = std::pow( 4*tVolume/M_PI, 1.0/2.0 );
            //std::cout<<he<<std::endl;

            return tVolume;
        }

//------------------------------------------------------------------------------
        /**
          * auto detect interpolation scheme
          */
        fem::Integration_Order get_auto_integration_order()
        {
            switch( mCell->get_geometry_type() )
            {
                case( mtk::Geometry_Type::LINE ) :
                    return fem::Integration_Order::BAR_3;
                    break;

                case( mtk::Geometry_Type::QUAD ) :
                     return fem::Integration_Order::QUAD_3x3;
                     break;

                case( mtk::Geometry_Type::HEX ) :
                    return fem::Integration_Order::HEX_3x3x3;
                    break;

                default :
                    MORIS_ERROR( false, " Element::get_auto_integration_order - not defined for this geometry type. ");
                    return Integration_Order::UNDEFINED;
                    break;
            }
        }

//------------------------------------------------------------------------------
        /**
         * auto detect full integration scheme
         */
        //FIXME: works for Lagrange only
        mtk::Interpolation_Order get_auto_interpolation_order()
        {
            switch( mCell->get_geometry_type() )
            {
                case( mtk::Geometry_Type::LINE ) :
                    switch( mCell->get_number_of_vertices() )
                    {
                       case( 2 ) :
                           return mtk::Interpolation_Order::LINEAR;
                           break;

                       case( 3 ) :
                           return mtk::Interpolation_Order::QUADRATIC;
                           break;

                       default :
                           MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for LINE and number of vertices. ");
                           return mtk::Interpolation_Order::UNDEFINED;
                           break;
                    }

                case( mtk::Geometry_Type::QUAD ) :
                    switch( mCell->get_number_of_vertices() )
                    {
                        case( 4 ) :
                            return mtk::Interpolation_Order::LINEAR;
                            break;

                        case( 8 ) :
                            return mtk::Interpolation_Order::SERENDIPITY;
                            break;

                        case( 9 ) :
                            return mtk::Interpolation_Order::QUADRATIC;
                            break;

                        case( 16 ) :
                            return mtk::Interpolation_Order::CUBIC;
                            break;

                        default :
                            MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for QUAD and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                            break;
                    }

                case( mtk::Geometry_Type::HEX ) :
                    switch( mCell->get_number_of_vertices() )
                    {
                        case( 8 ) :
                            return mtk::Interpolation_Order::LINEAR;
                            break;

                        case( 20 ) :
                            return mtk::Interpolation_Order::SERENDIPITY;
                            break;

                        case( 27 ) :
                            return mtk::Interpolation_Order::QUADRATIC;
                            break;

                        case( 64 ) :
                            return mtk::Interpolation_Order::CUBIC;
                            break;

                        default :
                            MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for HEX and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                            break;
                    }

                default :
                    MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for this geometry type. ");
                    return mtk::Interpolation_Order::UNDEFINED;
                    break;
            }
        }

//------------------------------------------------------------------------------
        /**
         * create the field interpolators for the element
         */
        virtual Cell< Field_Interpolator* >
        create_element_field_interpolators( Geometry_Interpolator* aGeometryInterpolator ) = 0;

//------------------------------------------------------------------------------
        /**
         * set the field interpolators coefficients
         */
        virtual void
        set_element_field_interpolators_coefficients( Cell< Field_Interpolator* > & aFieldInterpolators ) = 0;

//------------------------------------------------------------------------------
        /**
         * get the field interpolators for an IWG
         */
        virtual Cell< Field_Interpolator* >
        get_IWG_field_interpolators( IWG*                        & aIWG,
                                     Cell< Field_Interpolator* > & aFieldInterpolators ) = 0;

//------------------------------------------------------------------------------
        /**
         * set the initial sizes and values for mJacobianElement and mResidualElement
         */
        virtual void
        initialize_mJacobianElement_and_mResidualElement( Cell< Field_Interpolator* > & aFieldInterpolators ) = 0;

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */




#endif /* SRC_FEM_CL_FEM_ELEMENT_HPP_ */
