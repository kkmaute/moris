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
        /**
         * constructor
         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster cluster pointer from mtk mesh
         * @param[ in ] aNodes       cell of node pointers
         * @param[ in ] aSet         a fem set
         */
        Cluster( const Element_Type                aElementType,
                 const mtk::Cluster              * aMeshCluster,
                       moris::Cell< Node_Base* > & aNodes,
                       Set                       * aSet );
        Cluster(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Cluster();

//------------------------------------------------------------------------------
        /**
         * gets the IG cell local coordinates on the side wrt to the IP cell
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
        Matrix< DDRMat > get_side_normal( const mtk::Cell    * aCell,
                                          moris::moris_index   aSideOrdinal );

//------------------------------------------------------------------------------
        /**
         * gets the index of the vertex associated with a given master vertex
         * @param[ in ] aLeftVertex mesh vertex pointer
         */
        moris::mtk::Vertex const * get_left_vertex_pair( moris::mtk::Vertex const * aLeftVertex );

        moris::moris_index get_right_vertex_ordinal_on_facet( moris_index aCellIndexInCluster,
                                                             moris::mtk::Vertex const * aVertex );

//------------------------------------------------------------------------------
        /**
         * computes the jacobian on cluster
         */
        void compute_jacobian();

//------------------------------------------------------------------------------
        /**
         * computes the residual on cluster
         */
        void compute_residual();

//------------------------------------------------------------------------------
        /**
         * computes the jacobian and the residual on cluster
         */
        void compute_jacobian_and_residual();

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
        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                           moris::Cell< MSI::Dof_Type > aDofType )
        {
            // get pdofs values for the element
            this->get_my_pdof_values();

            // get a specific dof type pdofs values
            Matrix< DDRMat > tPdofValues;

            Cell< Matrix< DDRMat > > tPdofValues_Original;

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
    protected:
//------------------------------------------------------------------------------
        /**
         * computes the cluster volume
         */
        real compute_volume();

//------------------------------------------------------------------------------
        /**
         * sets the field interpolators coefficients
         */
        void set_field_interpolators_coefficients();

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
    private:
//------------------------------------------------------------------------------

         void reshape_pdof_values( const Cell< Matrix< DDRMat > > & aPdofValues,
                                         Matrix< DDRMat >         & aReshapedPdofValues);

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
