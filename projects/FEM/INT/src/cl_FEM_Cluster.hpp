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

    class Cluster
    {

    protected:

        // pointer to the mesh cluster
        const mtk::Cluster* mMeshCluster = nullptr;

        MSI::Equation_Object* mInterpolationElement = nullptr;

        // time sideset information
        Matrix< IndexMat > mListOfTimeOrdinals;

        // list of pointers to the master and slave mesh integration cells
        moris::Cell<const mtk::Cell *> mMasterIntegrationCells;
        moris::Cell<const mtk::Cell *> mSlaveIntegrationCells;

        // master and slave side ordinal information
        Matrix< IndexMat > mMasterListOfSideOrdinals;
        Matrix< IndexMat > mSlaveListOfSideOrdinals;

        // list of pointers to element
        moris::Cell< fem::Element * > mElements;

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

        /**
         * constructor
         * @param[ in ] aElementType    enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster    cluster pointer from mtk mesh
         * @param[ in ] aSet            a fem set
         * @param[ in ] aEquationObject pointer to the corresponding interpolation element
         */
        Cluster( const Element_Type                aElementType,
                 const mtk::Cluster              * aMeshCluster,
                       Set                       * aSet,
                       MSI::Equation_Object      * aEquationObject );

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster();

//------------------------------------------------------------------------------
        /**
         * get mesh cluster
         * @param[ out ] mMeshCluster a mesh cluster
         */
        const mtk::Cluster * get_mesh_cluster()
        {
            return mMeshCluster;
        }

//------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates on the side wrt to the IP cell
         * @param[ in ] aCellIndexInCluster index of the IG cell within the cluster
         * @param[ in ] aSideOrdinal        ordinal for the side
         * @param[ in ] aIsMaster           enum for master or slave
         */
        moris::Matrix< moris::DDRMat > get_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aCellIndexInCluster,
                                                                                      moris::moris_index aSideOrdinal,
                                                                                      mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates wrt IP cell
         * @param[ in ] aPrimaryCellIndexInCluster index of the IG cell within the cluster
         */
        moris::Matrix< moris::DDRMat > get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellIndexInCluster );

//------------------------------------------------------------------------------
        /**
         * get side normal
         * @param[ in ] aCell        mesh cell pointer
         * @param[ in ] aSideOrdinal ordinal of the side where normal is evaluated
         */
        Matrix< DDRMat > get_side_normal( const mtk::Cell          * aCell,
                                                moris::moris_index   aSideOrdinal );

//------------------------------------------------------------------------------
        /**
         * get the index of the vertex associated with a given master vertex
         * @param[ in ] aLeftVertex mesh vertex pointer
         */
        const moris::mtk::Vertex * get_left_vertex_pair( const moris::mtk::Vertex * aLeftVertex );

//------------------------------------------------------------------------------
        /**
         * get the ordinal of the right vertex on the facet
         * @param[ in ] aCellIndexInCluster an index for the cell in the cluster
         * @param[ in ] aVertex             a vertex pointer
         */
        moris::moris_index get_right_vertex_ordinal_on_facet(       moris_index          aCellIndexInCluster,
                                                              const moris::mtk::Vertex * aVertex);

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
        void compute_quantity_of_interest( const uint aMeshIndex,
                                           enum vis::Output_Type aOutputType,
                                           enum vis::Field_Type  aFieldType );

//------------------------------------------------------------------------------
        /**
         * get the nodal pdof values on an element
         * @param[ in ] aVertexIndex a vertex index
         * @param[ in ] aDofType     a dof type
         */
//        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
//                                           moris::Cell< MSI::Dof_Type > aDofType )
//        {
//            // get pdofs values for the element
//            this->compute_my_pdof_values();
//
//            // get a specific dof type pdofs values
//            Matrix< DDRMat > tPdofValues;
//
//            moris::Cell< Matrix< DDRMat > > tPdofValues_Original;
//
//            this->get_my_pdof_values( aDofType, tPdofValues_Original );
//
//            // reshape tCoeffs into the order the cluster expects them
//            this->reshape_pdof_values( tPdofValues_Original, tPdofValues );
//
//            // select the required nodal value
//            Matrix< IndexMat > tElemVerticesIndices = mMasterInterpolationCell->get_vertex_inds();
//            uint tElemNumOfVertices = mMasterInterpolationCell->get_number_of_vertices();
//
//            moris_index tVertexIndex = MORIS_INDEX_MAX;
//            for( uint i = 0; i < tElemNumOfVertices; i++ )
//            {
//                if ( tElemVerticesIndices( i ) == aVertexIndex )
//                {
//                    tVertexIndex =  i ;
//                    break;
//                }
//            }
//            return tPdofValues( tVertexIndex );
//        }
        
//------------------------------------------------------------------------------
        /**
         * compute the cluster volume
         */
        real compute_volume();

//------------------------------------------------------------------------------
        /*
         * Compute the measure (volume 3d or area 2d) of the cells in the void or primary phase
         */
        moris::real
        compute_cluster_cell_measure(const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                                     const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;
//------------------------------------------------------------------------------

        /*
         * Compute the side measure (surface area 3d or length 2d) of the cells in the void or primary phase on the side set.
         * Only valid on side cluster type mtk clusters
         */
        moris::real
        compute_cluster_cell_side_measure(const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                                          const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

//------------------------------------------------------------------------------
    protected:


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
