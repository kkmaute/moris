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

        // pointer to the visualization mesh cluster
        const mtk::Cluster * mVisMeshCluster = nullptr;

        // list of pointers to the master and slave mesh visualization cells
        moris::Cell< mtk::Cell const * > mMasterVisCells;
        moris::Cell< mtk::Cell const * > mSlaveVisCells;

        // master and slave side ordinal information
        Matrix< IndexMat > mMasterVisListOfSideOrdinals;
        Matrix< IndexMat > mSlaveVisListOfSideOrdinals;

        // list of pointers to element
        moris::Cell< fem::Element * > mVisElements;

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
        /**
         * trivial constructor
         */
        Cluster(){};

//        /**
//         * constructor
//         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
//         * @param[ in ] aMeshCluster cluster pointer from mtk mesh
//         * @param[ in ] aNodes       cell of node pointers
//         * @param[ in ] aSet         a fem set
//         */
//        Cluster( const Element_Type                aElementType,
//                 const mtk::Cluster              * aMeshCluster,
//                       moris::Cell< Node_Base* > & aNodes,
//                       Set                       * aSet );

        /**
         * constructor
         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster cluster pointer from mtk mesh
         * @param[ in ] aNodes       cell of node pointers for IP vertices
         * @param[ in ] aIGNodes     cell of node pointers for IG vertices
         * @param[ in ] aSet         a fem set
         */
        Cluster( const Element_Type                aElementType,
                 const mtk::Cluster              * aMeshCluster,
                       moris::Cell< Node_Base* > & aNodes,
                       moris::Cell< Node_Base* > & aIGNodes,
                       Set                       * aSet );

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster();

//------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates on the side wrt to the IP cell
         * @param[ in ] aCellIndexInCluster index of the IG cell within the cluster
         * @param[ in ] aSideOrdinal        ordinal for the side
         * @param[ in ] aIsMaster           enum for master or slave
         */
        moris::Matrix< moris::DDRMat > get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                      moris::moris_index aSideOrdinal,
                                                                                      mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER,
                                                                                      bool aIsVis = false );

//------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates wrt IP cell
         * @param[ in ] aPrimaryCellIndexInCluster index of the IG cell within the cluster
         */
        moris::Matrix< moris::DDRMat > get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellIndexInCluster,
                                                                                              bool aIsVis = false );

//------------------------------------------------------------------------------
        /**
         * get side normal
         * @param[ in ] aCell        mesh cell pointer
         * @param[ in ] aSideOrdinal ordinal of the side where normal is evaluated
         */
        Matrix< DDRMat > get_side_normal( const mtk::Cell    * aCell,
                                          moris::moris_index   aSideOrdinal );

//------------------------------------------------------------------------------
        /**
         * get the index of the vertex associated with a given master vertex
         * @param[ in ] aLeftVertex mesh vertex pointer
         */
        moris::mtk::Vertex const * get_left_vertex_pair( moris::mtk::Vertex const * aLeftVertex );

//------------------------------------------------------------------------------
        /**
         * get the ordinal of the right vertex on the facet
         * @param[ in ] aCellIndexInCluster an index for the cell in the cluster
         * @param[ in ] aVertex             a vertex pointer
         */
        moris::moris_index get_right_vertex_ordinal_on_facet( moris_index                aCellIndexInCluster,
                                                              moris::mtk::Vertex const * aVertex );

//------------------------------------------------------------------------------
        /**
         * compute the jacobian on cluster
         */
        void compute_jacobian();

//------------------------------------------------------------------------------
        /**
         * compute the residual on cluster
         */
        void compute_residual();

//------------------------------------------------------------------------------
        /**
         * compute the jacobian and the residual on cluster
         */
        void compute_jacobian_and_residual();

//------------------------------------------------------------------------------
        /**
         * compute the quantity of interest on cluster
         * @param[ in ] aOutputType an enum for output type
         * @param[ in ] aFieldType  an enum for computation/return type
         *                          GLOBAL, NODAL, ELEMENTAL
         */
        //void compute_quantity_of_interest( fem::QI_Compute_Type aQIComputeType );
        void compute_quantity_of_interest( enum vis::Output_Type aOutputType,
                                           enum vis::Field_Type  aFieldType );

//------------------------------------------------------------------------------
        /**
         * FIXME get rid of mNodalWeakBCs
         * gets the nodal values for the weak BCs
         */
        Matrix< DDRMat > & get_weak_bcs()
        {
            return mNodalWeakBCs;
        }

//------------------------------------------------------------------------------
        /**
         * get the nodal pdof values on an element
         * @param[ in ] aVertexIndex a vertex index
         * @param[ in ] aDofType     a dof type
         */
        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                           moris::Cell< MSI::Dof_Type > aDofType )
        {
            // get pdofs values for the element
            this->compute_my_pdof_values();

            // get a specific dof type pdofs values
            Matrix< DDRMat > tPdofValues;

            moris::Cell< Matrix< DDRMat > > tPdofValues_Original;

            this->get_my_pdof_values( aDofType, tPdofValues_Original );

            // reshape tCoeffs into the order the cluster expects them
            this->reshape_pdof_values( tPdofValues_Original, tPdofValues );

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
        /**
         * set visualization cluster
         * @param[ in ] aVisMeshCluster a pointer to a visualization mesh cluster
         */
        void set_visualization_cluster( const mtk::Cluster * aVisMeshCluster )
        {
            // set a visualization cluster
            mVisMeshCluster = aVisMeshCluster;

            // get the visualization cells
            mMasterVisCells = mVisMeshCluster->get_primary_cells_in_cluster();

            // create an element factory
            fem::Element_Factory tElementFactory;

            // set size for the visualization element list
            uint tNumMasterVisCells = mMasterVisCells.size();
            mVisElements.resize( tNumMasterVisCells, nullptr );

             // switch on the element type
             switch ( mElementType )
             {
                 case ( fem::Element_Type::BULK ) :
                 {
                     // loop over the visualization cells
                     for( moris::uint Ik = 0; Ik < tNumMasterVisCells; Ik++)
                     {
                         // create an element
                         mVisElements( Ik )
                         = tElementFactory.create_element( mElementType,
                                                           mMasterVisCells( Ik ),
                                                           mSet,
                                                           this,
                                                           Ik );
                     }
                     break;
                 }
                 case ( fem::Element_Type::SIDESET ) :
                 {
                     // set the side ordinals for the IG cells in the cluster
                     mMasterVisListOfSideOrdinals = mVisMeshCluster->get_cell_side_ordinals();

                     // loop over the visualization cells
                     for( moris::uint Ik = 0; Ik < tNumMasterVisCells; Ik++)
                     {
                         // create an element
                         mVisElements( Ik )
                         = tElementFactory.create_element( mElementType,
                                                           mMasterVisCells( Ik ),
                                                           mSet,
                                                           this,
                                                           Ik );
                     }
                     break;
                 }
                 case ( fem::Element_Type::DOUBLE_SIDESET ) :
                 {
                     // fill the slave visualization cells
                     mSlaveVisCells  = mVisMeshCluster->get_primary_cells_in_cluster( mtk::Master_Slave::SLAVE );

                     // set the side ordinals for the master and slave vis cells
                     mMasterVisListOfSideOrdinals = mVisMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::MASTER );
                     mSlaveVisListOfSideOrdinals  = mVisMeshCluster->get_cell_side_ordinals( mtk::Master_Slave::SLAVE );

                     // loop over the visualization cells
                     for( moris::uint Ik = 0; Ik < tNumMasterVisCells; Ik++)
                     {
                         // create an element
                         mVisElements( Ik )
                         = tElementFactory.create_element( mElementType,
                                                           mMasterVisCells( Ik ),
                                                           mSlaveVisCells( Ik ),
                                                           mSet,
                                                           this,
                                                           Ik );
                     }
                     break;
                 }
                 default :
                 {
                     MORIS_ERROR( false, "Cluster::set_visualization_cluster - undefined element type" );
                     break;
                 }
             }
        }

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
        /**
         * compute the cluster volume
         */
        real compute_volume();

//------------------------------------------------------------------------------
        /**
         * set the field interpolators coefficients
         */
        void set_field_interpolators_coefficients();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
