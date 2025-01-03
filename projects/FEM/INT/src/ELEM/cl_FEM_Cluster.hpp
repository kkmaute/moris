/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Cluster.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CLUSTER_HPP_
#define SRC_FEM_CL_FEM_CLUSTER_HPP_

#include "assert.h"
#include <cmath>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"                     //MTK/src
#include "cl_MSI_Equation_Object.hpp"          //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                    //FEM/INT/src
#include "cl_FEM_Node.hpp"                     //FEM/INT/src
#include "cl_FEM_IWG.hpp"                      //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"       //FEM/INT/src
#include "cl_MTK_Integrator.hpp"               //MTK/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_Set.hpp"                      //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp"    //FEM/INT/src
#include "cl_FEM_Cluster_Measure.hpp"          //FEM/INT/src

namespace moris::fem
{
    class Set;
    //------------------------------------------------------------------------------

    class Cluster
    {

      protected:
        // pointer to the mesh cluster
        const mtk::Cluster *mMeshCluster = nullptr;

        MSI::Equation_Object *mInterpolationElement = nullptr;

        // time sideset information
        Matrix< IndexMat > mListOfTimeOrdinals;

        // list of pointers to the leader and follower mesh integration cells
        Vector< const mtk::Cell * > mLeaderIntegrationCells;
        Vector< const mtk::Cell * > mFollowerIntegrationCells;

        // leader and follower side ordinal information
        Matrix< IndexMat > mLeaderListOfSideOrdinals;
        Matrix< IndexMat > mFollowerListOfSideOrdinals;

        Matrix< DDRMat > mQuadraturePoints;

        Matrix< DDRMat > mQuadratureWeights;

        // list of pointers to element
        Vector< fem::Element * > mElements;

        // flag for all IG element whether or not to compute residual and QI (and derivatives)
        Vector< bool > mComputeResidualAndIQI;

        // pointer to the fem set
        Set *mSet = nullptr;

        // element type
        Element_Type mElementType = Element_Type::END_ELEMENT_TYPE;

        // acceptable volume error
        // IG cells with a relative volume below this threshold are ignored
        const real mVolumeError = 1e-6;

        // cluster measures
        std::map< Cluster_Measure::ClusterMeasureSpecification, std::shared_ptr< Cluster_Measure > > mClusterMeasures;

        // flag to easily differentiate between ACTUAL FEM and VIS clusters for debugging
        bool mIsVisCluster = false;

        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element_Time_Boundary;
        friend class Element;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aElementType    enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster    cluster pointer from mtk mesh
         * @param[ in ] aSet            a fem set
         * @param[ in ] aEquationObject pointer to the corresponding interpolation element
         * @param[ in ] aClusterMeasureTuples cell of cluster measure specifications
         */
        Cluster(
                const Element_Type    aElementType,
                const mtk::Cluster   *aMeshCluster,
                Set                  *aSet,
                MSI::Equation_Object *aEquationObject,
                bool                  aIsVisCluster = false );

        //------------------------------------------------------------------------------

        /**
         * trivial constructor for testing purpose
         */
        Cluster();

        //------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Cluster();

        //------------------------------------------------------------------------------

        bool
        is_VIS_cluster() const
        {
            return mIsVisCluster;
        }

        //------------------------------------------------------------------------------

        Element_Type
        get_element_type() const
        {
            return mElementType;
        }

        Vector< fem::Element * >
        get_elements() const
        {
            return mElements;
        }

        //------------------------------------------------------------------------------

        /**
         * get mesh cluster
         * @param[ out ] mMeshCluster a mesh cluster
         */
        const mtk::Cluster *
        get_mesh_cluster()
        {
            return mMeshCluster;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > get_quadrature_points() const ;

        //------------------------------------------------------------------------------
          
        Matrix< DDRMat > get_quadrature_weights() const ;

        //------------------------------------------------------------------------------

        void set_quadrature_points() ;

        //------------------------------------------------------------------------------

        void set_quadrature_weights() ;

        //------------------------------------------------------------------------------
        /**
         * get side ordinal information
         * @param[ out ] mMeshCluster a mesh cluster
         */
        Matrix< IndexMat > &get_side_ordinal_info(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get the vertices local coordinates on the IP cell
         */
        moris::Matrix< moris::DDRMat > get_vertices_local_coordinates_wrt_interp_cell(
                mtk::Leader_Follower aLeaderFollower = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get the vertices indices in cluster for sensitivity
         * @param[ in ] aVerticesIndices matrix of vertex indices on cluster to be filled
         */
        void get_vertex_indices_in_cluster_for_sensitivity(
                moris::Matrix< moris::IndexMat > &aVerticesIndices );

        //------------------------------------------------------------------------------
        /**
         * get the vertices indices in cluster for visualization
         * @param[ in ] aVerticesIndices matrix of vertex indices on cluster to be filled
         */
        void get_vertex_indices_in_cluster_for_visualization(
                moris::Matrix< moris::IndexMat > &aVerticesIndices,
                mtk::Leader_Follower              aLeaderFollower = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates on the side wrt to the IP cell
         * @param[ in ] aCellIndexInCluster index of the IG cell within the cluster
         * @param[ in ] aSideOrdinal        ordinal for the side
         * @param[ in ] aIsLeader           enum for leader or follower
         */
        moris::Matrix< moris::DDRMat > get_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index   aCellIndexInCluster,
                moris::moris_index   aSideOrdinal,
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get the IG cell local coordinates wrt IP cell
         * @param[ in ] aPrimaryCellIndexInCluster index of the IG cell within the cluster
         */
        moris::Matrix< moris::DDRMat > get_primary_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index aPrimaryCellIndexInCluster );

        //------------------------------------------------------------------------------
        /**
         * get side normal
         * @param[ in ] aCell        mesh cell pointer
         * @param[ in ] aSideOrdinal ordinal of the side where normal is evaluated
         */
        Matrix< DDRMat > get_side_normal(
                const mtk::Cell   *aCell,
                moris::moris_index aSideOrdinal );

        //------------------------------------------------------------------------------
        /**
         * get the index of the vertex associated with a given leader vertex
         * @param[ in ] aLeftVertex mesh vertex pointer
         */
        const moris::mtk::Vertex *get_left_vertex_pair(
                const moris::mtk::Vertex *aLeftVertex );

        //------------------------------------------------------------------------------
        /**
         * get the ordinal of the right vertex on the facet
         * @param[ in ] aCellIndexInCluster an index for the cell in the cluster
         * @param[ in ] aVertex             a vertex pointer
         */
        moris::moris_index get_right_vertex_ordinal_on_facet(
                moris_index               aCellIndexInCluster,
                const moris::mtk::Vertex *aVertex );

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
         * @param[ in ] aFemMeshIndex mesh index for used IG mesh
         * @param[ in ] aFieldType enum for computation/return type
         *                         GLOBAL, NODAL, ELEMENTAL_INT, ELEMENTAL_AVG
         */
        void compute_quantity_of_interest(
                const uint           aFemMeshIndex,
                enum vis::Field_Type aFieldType );

        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest on cluster
         * @param[ in ] aValues field values
         * @param[ in ] aFieldType enum for computation/return type
         *                         ELEMENTAL
         * @param[ in ] aIQIIndex IQI index
         */
        void compute_quantity_of_interest(
                Matrix< DDRMat >      &aValues,
                mtk::Field_Entity_Type aFieldType,
                uint                   aIQIIndex,
                real                  &aSpaceTimeVolume );

        //------------------------------------------------------------------------------
        /**
         * compute dRdp by analytical formulation
         */
        void compute_dRdp();

        //------------------------------------------------------------------------------
        /**
         * compute the quantities of interest on cluster
         */
        void compute_QI();

        //------------------------------------------------------------------------------
        /**
         * compute dQIdp by analytical formulation
         */
        void compute_dQIdp_explicit();

        //------------------------------------------------------------------------------
        /**
         * compute dRdp and dQIdp by analytical formulation
         */
        void compute_dRdp_and_dQIdp();

        //------------------------------------------------------------------------------
        /**
         * compute dQIdu
         */
        void compute_dQIdu();

        //------------------------------------------------------------------------------
        /**
         * get cluster measure
         */
        std::shared_ptr< Cluster_Measure > get_cluster_measure(
                fem::Measure_Type    aMeasureType,
                mtk::Primary_Void    aIsPrimary,
                mtk::Leader_Follower aIsLeader );

        //------------------------------------------------------------------------------
        /**
         * get cluster measures
         */
        std::map< Cluster_Measure::ClusterMeasureSpecification, std::shared_ptr< Cluster_Measure > >
        get_cluster_measures()
        {
            return mClusterMeasures;
        }

        //------------------------------------------------------------------------------
        /**
         * reset cluster measure
         */
        void reset_cluster_measure();

        //------------------------------------------------------------------------------
        /**
         * reset cluster measure derivative
         */
        void reset_cluster_measure_derivatives();

        //------------------------------------------------------------------------------
        /*
         * Compute the measure (volume 3d or area 2d) of the cells in the void or primary phase
         */
        moris::real compute_cluster_cell_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /*
         * Compute the measure derivatives of the cells in the void or primary phase
         */
        moris::Matrix< DDRMat > compute_cluster_cell_measure_derivative(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * Compute the side measure (surface area 3d or length 2d) of the cells in the void or primary phase on the side set.
         * Only valid on side cluster type mtk clusters
         */
        moris::real compute_cluster_cell_side_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /*
         * Compute the side measure derivatives in the void or primary phase on the side set.
         * Only valid on side cluster type mtk clusters
         */
        moris::Matrix< DDRMat > compute_cluster_cell_side_measure_derivative(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * Compute the element size (length) of the cells in the void or primary phase
         */
        moris::real compute_cluster_cell_length_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /*
         * Compute the element size (length) derivatives of the cells in the void or primary phase
         */
        moris::Matrix< DDRMat > compute_cluster_cell_length_measure_derivative(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * Compute the ip element size (length) leader or follower
         * @param[ in ] aIsLeader enum for leader or follower
         */
        moris::real compute_ip_cell_length_measure(
                const mtk::Leader_Follower aIsLeader ) const;

        //------------------------------------------------------------------------------
        /**
         * compute the cluster volume by calling MTK cell function
         */
        real compute_volume();

        //------------------------------------------------------------------------------
        /**
         * compute the IG cells volumes by calling MTK cell function
         */
        Matrix< DDRMat > compute_element_volumes();

        //------------------------------------------------------------------------------
        /**
         * compute the cluster volume by numerical integrating IG elements in FEM
         */
        real compute_volume_in_fem();

        /**
         * @brief The cluster might have some mtk cells in their elements added multiple times to it (e.g. in a nonconformal set, a cell might have two different follower cells).
         * This method will return only the unique cells that have been used in the elements of the cluster.
         * @note  The (nonconformal) cluster contains cells [ 1, 2, 3 ] and the follower cluster contains cells [ 5, 6, 7 ].
         * If the leader-follower pairing is as follows (the elements are stored at the 'local index' position in the vector of elements):
         * @verbatim
         *  Local | Leader | Follower
         *  Index | Index  |  Index
         *    0       1   ->   5
         *    1       1   ->   6
         *    2       2   ->   6
         *    3       3   ->   7
         *  @endverbatim
         * Cells 1 and 6 will be used twice in the cluster (in the elements 0, 1 and 2) such that they would be counted twice e.g. in the volume computation.
         * For the leader elements, the method would return the following map:
         * @verbatim
         * Cell 1 -> [ 0, 1 ]
         * Cell 2 -> [ 2 ]
         * Cell 3 -> [ 3 ]
         * @endverbatim
         * And for the follower elements, the method would return the following map:
         * @verbatim
         * Cell 5 -> [ 0 ]
         * Cell 6 -> [ 1, 2 ]
         * Cell 7 -> [ 3 ]
         * @endverbatim
         * @param aLeaderFollowerType whether to return the unique leader or follower elements
         * @return A map from the cells in the cluster (unique) to the indices of the elements (possibly multiple) in the cluster that use the cell.
         */
        std::map< mtk::Cell const *, Vector< moris_index > > get_cell_to_element_map( mtk::Leader_Follower aLeaderFollowerType = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /**
         * compute volume for each element in cluster relative to total volume of all elements;
         * if volume is zero or negative the return matrix is filled with negative one
         *
         * @return vector with relative volumes
         */
        Matrix< DDRMat > compute_relative_volume();

        //------------------------------------------------------------------------------
        /**
         * compute elements to be considered/ignored for residual and IQI computation
         */
        void determine_elements_for_residual_and_iqi_computation();

        //------------------------------------------------------------------------------
        /*
         * Compute the threshold such that the total volume of all elements with a relative element volume
         * smaller than this threshold is less than a give volume error
         *
         * @param[ in ] tRelativeElementVolume vector of element volumes
         * @param[ in ] tVolumeError           acceptable volume erro
         *
         * @return threshold value
         */
        real compute_volume_drop_threshold(
                const Matrix< DDRMat > &tRelativeCellVolume,
                const real             &tVolumeError );

        //------------------------------------------------------------------------------
    };

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CLUSTER_HPP_ */
