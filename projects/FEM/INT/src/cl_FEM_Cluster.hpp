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

        //! pointer to cell on mesh
        const mtk::Cell * mInterpolationCell;

        moris::Cell< mtk::Cell const * > mIntegrationCell;

        moris::Cell< fem::Element * > mInterpElements;

        //! node indices of this element
        //  @node: MTK interface returns copy of vertices. T
        //         storing the indices in private matrix is faster,
        //         but might need more memory
        moris::Matrix< IndexMat > mNodeIndices;

        uint mNumOfIWGs;

        Set * mElementBlock;

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
                       Set             * aElementBlock) : MSI::Equation_Object( aElementBlock ),
                                                                    mElementBlock( aElementBlock )
        {
            // fill the bulk mtk::Cell pointer //FIXME
            mInterpolationCell = aCell;

            mIntegrationCell.resize( 1, mInterpolationCell );       //Fixme

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
            mNumOfIWGs = mElementBlock->get_num_IWG(); //FIXME

            fem::Element_Factory tElementFactory;

            mInterpElements.resize( mIntegrationCell.size(), nullptr );

            for( moris::uint Ik = 0; Ik < mInterpElements.size(); Ik++)
            {
                // create the element
                mInterpElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                        mIntegrationCell( Ik ),
                                                                        mElementBlock,
                                                                        this );
            }
        };
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster()
        {
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

            // get a specific dof type profs values
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
            uint tNumOfIntegPoints = mElementBlock->get_num_integration_points();

            // init volume
            real tVolume = 0;

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // compute integration point weight x detJ
                real tWStar = aGeometryInterpolator->det_J( mElementBlock->get_integration_points().get_column( iGP ) )
                            * mElementBlock->get_integration_weights()( iGP );

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
             for( uint i = 0; i < mElementBlock->get_num_interpolators(); i++ )
             {
                 // get the ith dof type group
                 Cell< MSI::Dof_Type > tDofTypeGroup = mElementBlock->get_interpolator_dof_type_list()( i );

                 //FIXME:forced coefficients
                 // get the pdof values for the ith dof type group
                 Matrix< DDRMat > tCoeff;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff );

                 // set the field coefficients
                 mElementBlock->get_block_field_interpolator()( i )->set_coeff( tCoeff );
             }
         }

 //------------------------------------------------------------------------------
        /**
         * @Brief set the initial sizes and values for mJacobianElement
         */
         void initialize_mJacobianElement();

//------------------------------------------------------------------------------

         /**
          * @Brief set the initial sizes and values formResidualElement
          */
         void initialize_mResidualElement();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
