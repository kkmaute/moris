/*
 * cl_FEM_Cluster.hpp
 *
 *  Created on: Apr 20, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_CLUSTER_HPP_
#define SRC_FEM_CL_FEM_CLUSTER_HPP_

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

#include "cl_FEM_Element_Factory.hpp"            //FEM/INT/src
#include "cl_FEM_Set.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Set;
//------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Cluster : public MSI::Equation_Object
    {

    protected:

        // pointer to the mesh interpolation cell
        const mtk::Cell * mInterpolationCell;

        // list of pointers to the mesh integration cells
        moris::Cell< mtk::Cell const * > mIntegrationCells;

        // list of pointers to element
        moris::Cell< fem::Element * > mElements;

        //! node indices of this element
        //  @node: MTK interface returns copy of vertices. T
        //         storing the indices in private matrix is faster,
        //         but might need more memory
        moris::Matrix< IndexMat > mNodeIndices;

        uint mNumOfIWGs;

        Set * mSet;

        Element_Type mElementType;

        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Cluster( const Element_Type                aElementType,
                 const mtk::Cell                 * aCell,
                       moris::Cell< Node_Base* > & aNodes,
                       Set                       * aElementBlock) : MSI::Equation_Object( aElementBlock ),
                                                                    mSet( aElementBlock )
        {
            // fill the bulk mtk::Cell pointer //FIXME
            mInterpolationCell = aCell;

            mIntegrationCells.resize( 1, mInterpolationCell );       //Fixme

            mElementType = aElementType;

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

            // get the number of IWGs
            mNumOfIWGs = mSet->get_num_IWG(); //FIXME

            fem::Element_Factory tElementFactory;

            mElements.resize( mIntegrationCells.size(), nullptr );

            for( moris::uint Ik = 0; Ik < mElements.size(); Ik++)
            {
                // create the element
                mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                  mIntegrationCells( Ik ),
                                                                  mSet,
                                                                  this );
            }
        };

        Cluster( const Element_Type                aElementType,
                 const mtk::Cell_Cluster         * aCellCluster,
                       moris::Cell< Node_Base* > & aNodes,
                       Set                       * aSet ) : MSI::Equation_Object( aSet ),
                                                            mSet( aSet )
        {
            // fill the interpolation cell
            mInterpolationCell = & aCellCluster->get_interpolation_cell();

            // fill the integration cells
            mIntegrationCells = aCellCluster->get_primary_cells_in_cluster();

            // fill the element type
            mElementType = aElementType;

            // select the element nodes from aNodes and fill mNodeObj
            // get vertices from cell
            moris::Cell< mtk::Vertex* > tVertices = mInterpolationCell->get_vertex_pointers();

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

            // get the number of IWGs
            mNumOfIWGs = mSet->get_num_IWG(); //FIXME

            // element factory
            fem::Element_Factory tElementFactory;

            mElements.resize( mIntegrationCells.size(), nullptr );

            for( moris::uint Ik = 0; Ik < mIntegrationCells.size(); Ik++)
            {
                // create an element
                mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                  mIntegrationCells( Ik ),
                                                                  mSet,
                                                                  this );
            }
        };

        Cluster( const Element_Type                aElementType,
                 const mtk::Side_Cluster         * aSideCluster,
                       moris::Cell< Node_Base* > & aNodes,
                       Set                       * aSet ) : MSI::Equation_Object( aSet ),
                                                            mSet( aSet )
        {
            // fill the interpolation cell
            mInterpolationCell = & aSideCluster->get_interpolation_cell();

            // fill the integration cells
            mIntegrationCells = aSideCluster->get_cells_in_side_cluster();

            // set the side ordinals for the IG cells in the cluster
            mListOfSideOrdinals = aSideCluster->get_cell_side_ordinals();

            // fill the element type
            mElementType = aElementType;

            // select the element nodes from aNodes and fill mNodeObj
            // get vertices from cell
            moris::Cell< mtk::Vertex* > tVertices = mInterpolationCell->get_vertex_pointers();

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

            // get the number of IWGs
            mNumOfIWGs = mSet->get_num_IWG(); //FIXME

            // element factory
            fem::Element_Factory tElementFactory;

            mElements.resize( mIntegrationCells.size(), nullptr );

            for( moris::uint Ik = 0; Ik < mIntegrationCells.size(); Ik++)
            {
                // create an element
                mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                  mIntegrationCells( Ik ),
                                                                  mSet,
                                                                  this );
            }
        };
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster()
        {
            mElements.clear();
        };

//------------------------------------------------------------------------------

        void compute_jacobian();

//------------------------------------------------------------------------------

        void compute_residual();

//------------------------------------------------------------------------------

        void compute_jacobian_and_residual() {};

//------------------------------------------------------------------------------

        Matrix< DDRMat > & get_weak_bcs()
        {
            return mNodalWeakBCs;
        }

//------------------------------------------------------------------------------

        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                           moris::Cell< MSI::Dof_Type > aDofType )
        {
            // get pdofs values for the element
            this->get_my_pdof_values();

            // get a specific dof type pdofs values
            Matrix< DDRMat > tPdofValues;
            this->get_my_pdof_values( aDofType, tPdofValues );

            // select the required nodal value
            Matrix< IndexMat > tElemVerticesIndices = mInterpolationCell->get_vertex_inds();
            uint tElemNumOfVertices = mInterpolationCell->get_number_of_vertices();

            moris_index tVertexIndex = MORIS_INDEX_MAX;
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
            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_num_integration_points();

            // init volume
            real tVolume = 0;

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // compute integration point weight x detJ
                real tWStar = aGeometryInterpolator->det_J( mSet->get_integration_points().get_column( iGP ) )
                            * mSet->get_integration_weights()( iGP );

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
         * set the field interpolators coefficients
         */
        void set_field_interpolators_coefficients( )
         {
             // loop on the dof types
             for( uint i = 0; i < mSet->get_num_interpolators(); i++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_interpolator_dof_type_list()( i );

                 //FIXME:forced coefficients
                 // get the pdof values for the ith dof type group
                 Matrix< DDRMat > tCoeff;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff );

                 // set the field coefficients
                 mSet->get_block_field_interpolator()( i )->set_coeff( tCoeff );
             }
         }

 //------------------------------------------------------------------------------
        /**
         * @Brief set the initial sizes and values for mJacobian
         */
         void initialize_mJacobian();

//------------------------------------------------------------------------------

         /**
          * @Brief set the initial sizes and values for mResidual
          */
         void initialize_mResidual();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
