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

        // time sideset information
        Matrix< IndexMat > mListOfTimeOrdinals;

        // pointer to the master and slave mesh interpolation cell
        const mtk::Cell * mMasterInterpolationCell;
        const mtk::Cell * mSlaveInterpolationCell;

        // list of pointers to the master and slave mesh integration cells
        moris::Cell< mtk::Cell const * > mMasterIntegrationCells;
        moris::Cell< mtk::Cell const * > mSlaveIntegrationCells;

        // master and slave side ordinal information
        Matrix< IndexMat > mMasterListOfSideOrdinals;
        Matrix< IndexMat > mSlaveListOfSideOrdinals;

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

            // fill the master interpolation cell
            mMasterInterpolationCell = & aMeshCluster->get_interpolation_cell();

            // fill the master integration cells
            mMasterIntegrationCells = aMeshCluster->get_primary_cells_in_cluster();

            // get the number of IWGs //FIXME
            mNumOfIWGs = mSet->get_num_IWG();

            // switch on the element type
            switch ( mElementType )
            {
                case ( fem::Element_Type::BULK ):
                {
                    // select the element nodes from aNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tVertices = mMasterInterpolationCell->get_vertex_pointers();

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

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik);
                    }
                    break;
                }
                case ( fem::Element_Type::SIDESET ):
                {
                    // set the side ordinals for the IG cells in the cluster
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals();

                    // select the element nodes from aNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tVertices = mMasterInterpolationCell->get_vertex_pointers();

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

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik );

                    }
                    break;
                }
                case( fem::Element_Type::DOUBLE_SIDESET ):
                {
                    // fill the slave interpolation cell
                    mSlaveInterpolationCell  = & aMeshCluster->get_interpolation_cell( mtk::Master_Slave::SLAVE );

                    // fill the slave integration cells
                    mSlaveIntegrationCells  = aMeshCluster->get_primary_cells_in_cluster( mtk::Master_Slave::SLAVE );

                    // set the side ordinals for the master and slave IG cells
                    mMasterListOfSideOrdinals = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::MASTER );
                    mSlaveListOfSideOrdinals  = aMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::SLAVE );

                    // select the element nodes from aIPNodes and fill mNodeObj
                    // get vertices from cell
                    moris::Cell< mtk::Vertex* > tMasterVertices = mMasterInterpolationCell->get_vertex_pointers();
                    moris::Cell< mtk::Vertex* > tSlaveVertices  = mSlaveInterpolationCell->get_vertex_pointers();

                    // get number of nodes from cell
                    uint tNumOfNodes = tMasterVertices.size() + tSlaveVertices.size();

                    // assign node object
                    mNodeObj.resize( tNumOfNodes, nullptr );

                    // fill node objects
                    uint tNodeCounter = 0;
                    for( uint i = 0; i < tMasterVertices.size(); i++)
                    {
                        mNodeObj( tNodeCounter ) = aNodes( tMasterVertices( i )->get_index() );
                        tNodeCounter++;
                    }
                    for( uint i = 0; i < tSlaveVertices.size(); i++)
                    {
                        mNodeObj( tNodeCounter ) = aNodes( tSlaveVertices( i )->get_index() );
                        tNodeCounter++;
                    }

                    // set size of Weak BCs
                    mNodalWeakBCs.set_size( tNumOfNodes, 1 );

                    // element factory
                    fem::Element_Factory tElementFactory;

                    mElements.resize( mMasterIntegrationCells.size(), nullptr );

                    for( moris::uint Ik = 0; Ik < mMasterIntegrationCells.size(); Ik++)
                    {
                        // create an element
                        mElements( Ik ) = tElementFactory.create_element( mElementType,
                                                                          mMasterIntegrationCells( Ik ),
                                                                          mSlaveIntegrationCells( Ik ),
                                                                          mSet,
                                                                          this,
                                                                          Ik );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "Cluster::Cluster - No element type specified" );
                    break;
            }
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster();

//------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                      moris::moris_index aSideOrdinal,
                                                                                      mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            // check we are working with a side cluster
            MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");
            MORIS_ASSERT( (( mElementType == fem::Element_Type::DOUBLE_SIDESET ) || ( mElementType == fem::Element_Type::SIDESET )),
                          "Cluster::get_cell_local_coords_on_side_wrt_interp_cell - not a side or double side cluster.");

            // is trivial master or slave?
            bool tIsTrivial = mSet->mIsTrivialMaster;
            if ( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                tIsTrivial = mSet->mIsTrivialSlave;
            }

            // if trivial cluster IP cell = IG cell
            if( tIsTrivial )
            {
                // get the side param coords from the IG geometry interpolator
                return mSet->get_IP_geometry_interpolator( aIsMaster )->extract_space_side_space_param_coeff( aSideOrdinal,
                                                                                                              mSet->get_IG_space_interpolation_order() );
            }
            // if non trivial cluster
            else
            {
                // get the side param coords from the side cluster
                return mMeshCluster->get_cell_local_coords_on_side_wrt_interp_cell( aCellIndexInCluster, aIsMaster );
            }
        }

//------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat> get_primary_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aPrimaryCellIndexInCluster )
        {
            // check we are working with a bulk cluster
            MORIS_ASSERT( mMeshCluster != NULL, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - empty cluster.");
            MORIS_ASSERT( mElementType == fem::Element_Type::BULK, "Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell - not a bulk cluster.");

            // if trivial cluster IP cell = IG cell
            if( mSet->mIsTrivialMaster )
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

        Matrix< DDRMat > get_side_normal( const mtk::Cell        * aCell,
                                          moris::moris_index aSideOrdinal )
        {
            // init normal
            Matrix < DDRMat > tNormal;

            // if interpolation cell is linear
            if( mSet->get_IG_space_interpolation_order() == mtk::Interpolation_Order::LINEAR )
            {
                // get normal from the mesh
                tNormal = aCell->compute_outward_side_normal( aSideOrdinal );
            }
            // if integration cell is higher order
            else
            {
                // get normal from the integration cell geometry interpolator
                mSet->get_IG_geometry_interpolator( mtk::Master_Slave::MASTER )->get_normal( tNormal );
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
            return mMeshCluster->get_vertex_cluster_index( tRightVertex, mtk::Master_Slave::SLAVE );

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
            Matrix< IndexMat > tElemVerticesIndices = mMasterInterpolationCell->get_vertex_inds();
            uint tElemNumOfVertices = mMasterInterpolationCell->get_number_of_vertices();

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
                // set integration point for geometry interpolator
                aGeometryInterpolator->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute integration point weight x detJ
                real tWStar = aGeometryInterpolator->det_J() * mSet->get_integration_weights()( iGP );

                // add contribution to jacobian from evaluation point
                //FIXME: include a thickness if 2D
                tVolume += tWStar;
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

                 switch ( mElementType )
                 {
                     case ( fem::Element_Type::BULK ) :
                     case ( fem::Element_Type::SIDESET ):
                     {
                         // set the field coefficients
                         mSet->get_field_interpolator()( i )->set_coeff( tCoeff );
                         break;
                     }
                     case ( fem::Element_Type::DOUBLE_SIDESET ) :
                     {
                         // set the field coefficients for the master and the slave interpolation cell
                         mSet->get_field_interpolator( mtk::Master_Slave::MASTER )( i )->set_coeff( tCoeff({0,tCoeff.numel()/2 -1},{0,0}) );
                         mSet->get_field_interpolator( mtk::Master_Slave::SLAVE )( i )->set_coeff( tCoeff({tCoeff.numel()/2, tCoeff.numel()-1},{0,0}) );
                         break;
                     }
                     default :
                     {
                         MORIS_ERROR( false, "Cluster::set_field_interpolators_coefficients - unknown element type ");
                         break;
                     }
                 }
             }
         }

//        /**
//         * set the field interpolators coefficients
//         */
//        void set_field_interpolators_coefficients_double( )
//         {
//             // loop on the dof types
//             for( uint i = 0; i < mSet->get_num_interpolators(); i++ )
//             {
//                 // get the ith dof type group
//                 moris::Cell< MSI::Dof_Type > tDofTypeGroup = mSet->get_interpolator_dof_type_list()( i );
//
//                 // get the pdof values for the ith dof type group
//                 Matrix< DDRMat > tCoeff;
//                 this->get_my_pdof_values( tDofTypeGroup, tCoeff );
//
//                 // FIXME ok only if both left and right have same interpolation
//                 // set the field coefficients for the left interpolation cell
//                 mSet->get_field_interpolator( mtk::Master_Slave::MASTER )( i )->set_coeff( tCoeff({0,tCoeff.numel()/2 -1},{0,0}) );
//
//                 // set the field coefficients for the right interpolation cell
//                 mSet->get_field_interpolator( mtk::Master_Slave::SLAVE )( i )->set_coeff( tCoeff({tCoeff.numel()/2, tCoeff.numel()-1},{0,0}) );
//             }
//         }

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
