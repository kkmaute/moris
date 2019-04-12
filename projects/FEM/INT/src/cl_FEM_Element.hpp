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
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"                  //MTK/src
#include "cl_MSI_Equation_Object.hpp"       //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Integrator.hpp"            //FEM/INT/src

#include "cl_FEM_Element_Block.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Element_Block;
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
        uint                                 mNumOfInterp;

        moris::Matrix< DDSMat >              mInterpDofTypeMap;

        Element_Block * mElementBlock;

        bool mIsMaster = false;
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

            mIsMaster = true;

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

            //create an element active dof type list from IWGs----------------------
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

            // create a list of the groups of dof types provided by the IWGs----------------
            // FIXME works as long as the dof type are always grouped in the same way
            moris::Cell< MSI::Dof_Type > tInterpDofTypeListBuild( mNumOfIWGs );

            // loop over the IWGs
            for ( uint i = 0; i < mNumOfIWGs; i++ )
            {
                // get the first dof type of each group
                tInterpDofTypeListBuild( i ) = mIWGs( i )->get_residual_dof_type()( 0 );
            }

            // get a unique list of the first dof type of each group
            Cell<moris::moris_index> tUniqueDofTypeGroupsIndices = unique_index( tInterpDofTypeListBuild );

            // get the number of unique dof type groups
            uint tNumOfUniqueDofTypeGroupsIndices = tUniqueDofTypeGroupsIndices.size();

            // set the size of the list of unique dof type groups
            mInterpDofTypeList.resize( tNumOfUniqueDofTypeGroupsIndices );

            // loop over the list of unique dof type groups
            for ( uint i = 0; i < tNumOfUniqueDofTypeGroupsIndices; i++ )
            {
                // get the unique residual dof type groups
                mInterpDofTypeList( i ) = mIWGs( tUniqueDofTypeGroupsIndices( i ) )->get_residual_dof_type();
            }

            // create a map of the element active dof type list------------------------
            // set number of unique pdof type of the element
            mNumOfInterp = tNumOfUniqueDofTypeGroupsIndices;

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
        };

        Element( mtk::Cell                 * aCell,
                 moris::Cell< IWG* >       & aIWGs,
                 moris::Cell< Node_Base* > & aNodes,
                 Element_Block      * aElementBlock) : mElementBlock(aElementBlock)
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

            // get the number of IWGs
            mNumOfIWGs = mIWGs.size();

            mGeometryInterpolator = mElementBlock->get_block_geometry_interpolator();
            mFieldInterpolators   = mElementBlock->get_block_field_interpolator();
            mEqnObjDofTypeList    = mElementBlock->get_unique_dof_type_list();
            mNumOfInterp          = mElementBlock->get_num_interpolators();
            mInterpDofTypeList    = mElementBlock->get_interpolator_dof_type_list();
            mInterpDofTypeMap     = mElementBlock->get_interpolator_dof_type_map();
        };
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
         ~Element()
        {
            if(mIsMaster)
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
        };

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
        fem::Integration_Order get_auto_integration_order( const mtk::Geometry_Type aGeometryType )
        {
            switch( aGeometryType )
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

                case( mtk::Geometry_Type::TRI ) :
                    return fem::Integration_Order::TRI_6;
                    break;

                case( mtk::Geometry_Type::TET ) :
                    return fem::Integration_Order::TET_5;
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

                case( mtk::Geometry_Type::TRI ) :
                    switch( mCell->get_number_of_vertices() )
                    {
                        case( 3 ) :
                            return mtk::Interpolation_Order::LINEAR;
                            break;

                        case( 6 ) :
                            return mtk::Interpolation_Order::QUADRATIC;
                            break;

                        case( 10 ) :
                            return mtk::Interpolation_Order::CUBIC;
                            break;

                        default :
                            MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TRI and number of vertices. ");
                            return mtk::Interpolation_Order::UNDEFINED;
                            break;
                    }

                case( mtk::Geometry_Type::TET ) :
                    switch( mCell->get_number_of_vertices() )
                    {
                        case( 4 ) :
                            return mtk::Interpolation_Order::LINEAR;
                            break;

                        case( 10 ) :
                            return mtk::Interpolation_Order::QUADRATIC;
                            break;

                        case( 20 ) :
                            return mtk::Interpolation_Order::CUBIC;
                            break;

                        default :
                            MORIS_ERROR( false, " Element::get_auto_interpolation_order - not defined for TRI and number of vertices. ");
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
        Cell< Field_Interpolator* > create_field_interpolators( Geometry_Interpolator* aGeometryInterpolator )
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
        /**
         * set the field interpolators coefficients
         */
        void set_field_interpolators_coefficients
             ( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
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
        /**
         * set the initial sizes and values for mJacobianElement and mResidualElement
         */
         void initialize_mJacobianElement_and_mResidualElement
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

//             std::cout<<tTotalDof<<std::endl;
             mJacobian.set_size( tTotalDof, tTotalDof, 0.0 );
             mResidual.set_size( tTotalDof, 1, 0.0 );
         }


 //------------------------------------------------------------------------------
         /**
          * get the field interpolators for an IWG
          */
         moris::Cell< Field_Interpolator* > get_IWG_field_interpolators
             ( IWG*                               & aIWG,
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
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_HPP_ */
