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
#include "cl_FEM_Element_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_Set.hpp"                   //FEM/INT/src

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
         // pointer to the mesh cluster
         const mtk::Cluster* mMeshCluster = nullptr;

//        // pointer to the mesh cell cluster
//        const mtk::Cell_Cluster* mCellCluster = nullptr;
//
//        // pointer to the mesh side cluster
//        const mtk::Side_Cluster* mSideCluster = nullptr;

        // pointer to the mesh interpolation cell
        const mtk::Cell * mInterpolationCell;

        // list of pointers to the mesh integration cells
        moris::Cell< mtk::Cell const * > mIntegrationCells;

        // sideset information
        Matrix< IndexMat > mListOfSideOrdinals;
        Matrix< IndexMat > mListOfTimeOrdinals;

        // list of pointers to element
        moris::Cell< fem::Element * > mElements;

        //! node indices of this element
        //  @node: MTK interface returns copy of vertices. T
        //         storing the indices in private matrix is faster,
        //         but might need more memory
        moris::Matrix< IndexMat > mNodeIndices;

        // number of IWG on the cluster
        uint mNumOfIWGs;

        // pointer to the fem set
        Set * mSet;

        // element type
        Element_Type mElementType;

//        // mesh double side cluster
//        mtk::Double_Side_Cluster mDoubleSideCluster;

        // pointer to the left and right mesh interpolation cell
        const mtk::Cell * mLeftInterpolationCell;
        const mtk::Cell * mRightInterpolationCell;

        // list of pointers to the left and right mesh integration cells
        moris::Cell< mtk::Cell const * > mLeftIntegrationCells;
        moris::Cell< mtk::Cell const * > mRightIntegrationCells;

        // left and right side ordinal information
        Matrix< IndexMat > mLeftListOfSideOrdinals;
        Matrix< IndexMat > mRightListOfSideOrdinals;

        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        Cluster( const Element_Type                aElementType,
                 const mtk::Cluster              * aMeshCluster,
                       moris::Cell< Node_Base* > & aNodes,
                       Set                       * aSet ) : MSI::Equation_Object( aSet ),
                                                            mSet( aSet ),
                                                            mElementType( aElementType )
        {
            // fill the cell cluster pointer
            mMeshCluster = aMeshCluster;

            if( mElementType == fem::Element_Type::BULK )
            {
                // fill the interpolation cell
                mInterpolationCell = & aMeshCluster->get_interpolation_cell();

                // fill the integration cells
                mIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

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
                                                                      this,
                                                                      Ik);
                }
            }
            else if( mElementType == fem::Element_Type::SIDESET )
            {
                // fill the interpolation cell
                mInterpolationCell = & aMeshCluster->get_interpolation_cell();

                // fill the integration cells
                mIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

                // set the side ordinals for the IG cells in the cluster
                mListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals();

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
                                                                      this,
                                                                      Ik );

                }
            }
            else if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // fill the left and right interpolation cell
                mLeftInterpolationCell  = & aMeshCluster->get_interpolation_cell( 0 );
                mRightInterpolationCell = & aMeshCluster->get_interpolation_cell( 1 );

                // fill the left and right integration cells
                mLeftIntegrationCells  = aMeshCluster->get_primary_cells_in_cluster( 0 );
                mRightIntegrationCells = aMeshCluster->get_primary_cells_in_cluster( 1 );

                // set the side ordinals for the left and right IG cells
                mLeftListOfSideOrdinals  = aMeshCluster->get_cell_side_ordinals( 0 );
                mRightListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals( 1 );

                // select the element nodes from aIPNodes and fill mNodeObj
                // get vertices from cell
                moris::Cell< mtk::Vertex* > tLeftVertices  = mLeftInterpolationCell->get_vertex_pointers();
                moris::Cell< mtk::Vertex* > tRightVertices = mRightInterpolationCell->get_vertex_pointers();

                // get number of nodes from cell
                uint tNumOfNodes = tLeftVertices.size() + tRightVertices.size();

                // assign node object
                mNodeObj.resize( tNumOfNodes, nullptr );

                // fill node objects
                uint tNodeCounter = 0;
                for( uint i = 0; i < tLeftVertices.size(); i++)
                {
                    mNodeObj( tNodeCounter ) = aNodes( tLeftVertices( i )->get_index() );
                    tNodeCounter++;
                }
                for( uint i = 0; i < tRightVertices.size(); i++)
                {
                    mNodeObj( tNodeCounter ) = aNodes( tRightVertices( i )->get_index() );
                    tNodeCounter++;
                }

                // set size of Weak BCs
                mNodalWeakBCs.set_size( tNumOfNodes, 1 );

                // get the number of IWGs
                mNumOfIWGs = mSet->get_num_IWG(); //FIXME

                // element factory
                fem::Element_Factory tElementFactory;

                mElements.resize( mLeftIntegrationCells.size(), nullptr );

                for( moris::uint Ik = 0; Ik < mLeftIntegrationCells.size(); Ik++)
                {
                    // create an element
                    mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                      mLeftIntegrationCells( Ik ),
                                                                      mRightIntegrationCells( Ik ),
                                                                      mSet,
                                                                      this,
                                                                      Ik );
                }
            }
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster();

//------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                      moris::moris_index aSideOrdinal )
        {
            // check we are working with a side cluster
            MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - not a side cluster.");


            // if trivial cluster IP cell = IG cell
            if( mMeshCluster->is_trivial() )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_IP_geometry_interpolator()->extract_space_side_space_param_coeff( aSideOrdinal,
                                                                                                   mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                // get the side param coords from the side cluster
                return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster );
            }
        }

//------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat> get_primary_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aPrimaryCellIndexInCluster )
        {

            // check we are working with a side cluster
            MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a cell cluster.");

            // if trivial cluster IP cell = IG cell
            if( mMeshCluster->is_trivial() )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_IP_geometry_interpolator()->extract_space_param_coeff( mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                // get the side param coords from the side cluster
                return mMeshCluster->get_primary_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellIndexInCluster );
            }
        }

//------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > get_left_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                           moris::moris_index aSideOrdinal )
        {
            // check we are working with a side cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                          "Cluster::get_left_cell_local_coords_on_side_wrt_interp_cell - not a double side cluster.");

            // if trivial cluster IP cell = IG cell
            if( mMeshCluster->is_trivial( 0 ) )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_left_IP_geometry_interpolator()->extract_space_side_space_param_coeff( aSideOrdinal,
                                                                                                        mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                // get the side param coords from the side cluster
                return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, 0 );
            }
        }

        moris::Matrix< moris::DDRMat > get_right_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                            moris::moris_index aSideOrdinal )
         {
             // check we are working with a side cluster
             MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                           "Cluster::get_right_cell_local_coords_on_side_wrt_interp_cell - not a double side cluster.");

             // if trivial cluster IP cell = IG cell
             if( mMeshCluster->is_trivial( 1 ) )
             {
                 // get the side param coords from the IG geometry interpolator
                 return mSet->get_right_IP_geometry_interpolator()->extract_space_side_space_param_coeff( aSideOrdinal,
                                                                                                          mSet->get_IG_space_interpolation_order() );
             }
             // if non trivial cluster
             else
             {
                 // get the side param coords from the side cluster
                 return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, 1 );
             }
         }

//------------------------------------------------------------------------------

        Matrix< DDRMat > get_side_normal( const mtk::Cell        * aCell,
                                          moris::moris_index aSideOrdinal,
                                          Matrix< DDRMat > & aParamPoint )
        {
            // init normal
            Matrix < DDRMat > tNormal;

            // if interpolation cell is linear
            if( mSet->get_IG_space_interpolation_order() == mtk::Interpolation_Order::LINEAR )
            {
                // egt normal from the mesh
                tNormal = aCell->compute_outward_side_normal( aSideOrdinal );
            }
            // if integration cell is higher order
            else
            {
                if( mElementType == fem::Element_Type::SIDESET )
                {
                    // get normal from the integration cell geometry interpolator
                    tNormal = mSet->get_IG_geometry_interpolator()->get_normal( aParamPoint );
                }
                else if( mElementType == fem::Element_Type::DOUBLE_SIDESET )
                {
                    // get normal from the integration cell geometry interpolator
                    tNormal = mSet->get_left_IG_geometry_interpolator()->get_normal( aParamPoint );
                }
            }

            return tNormal;

		}

//------------------------------------------------------------------------------

        moris::moris_index get_left_vertex_pair( moris::mtk::Vertex const * aLeftVertex )
        {
            // check we are working with a side cluster
            MORIS_ASSERT( mElementType == fem::Element_Type::DOUBLE_SIDESET,
                          "Cluster::get_left_vertex_pair - not a double side cluster.");

            // get the paired vertex on the right
            moris::mtk::Vertex const * tRightVertex = mMeshCluster->get_left_vertex_pair( aLeftVertex );

            // return the index of the paired vertex on the right
            return mMeshCluster->get_vertex_cluster_index( tRightVertex, 1 );

        }

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

                 // get the pdof values for the ith dof type group
                 Matrix< DDRMat > tCoeff;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff );

                 // set the field coefficients
                 mSet->get_field_interpolator()( i )->set_coeff( tCoeff );
             }
         }

        /**
         * set the field interpolators coefficients
         */
        void set_field_interpolators_coefficients_double( )
         {
             // loop on the dof types
             for( uint i = 0; i < mSet->get_num_interpolators(); i++ )
             {
                 // get the ith dof type group
                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_interpolator_dof_type_list()( i );

                 // get the pdof values for the ith dof type group
                 Matrix< DDRMat > tCoeff;
                 this->get_my_pdof_values( tDofTypeGroup, tCoeff );

                 // FIXME ok only if both left and right have same interpolation
                 // set the field coefficients for the left interpolation cell
                 mSet->get_left_field_interpolator()( i )->set_coeff( tCoeff({0,tCoeff.numel()/2 -1},{0,0}) );

                 // set the field coefficients for the right interpolation cell
                 mSet->get_right_field_interpolator()( i )->set_coeff( tCoeff({tCoeff.numel()/2, tCoeff.numel()-1},{0,0}) );
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
