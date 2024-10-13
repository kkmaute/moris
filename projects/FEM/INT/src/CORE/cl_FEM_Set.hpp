/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Set.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_SET_HPP_
#define SRC_FEM_CL_FEM_SET_HPP_

#include <set>

#include "assert.h"
#include "cl_Communication_Tools.hpp"

// MTK
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include <cl_FEM_Cluster_Measure.hpp>
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_FEM_IQI.hpp"
// FEM/MSI/src
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"
// FEM/VIS/src
#include "cl_VIS_Output_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        class Cell;
        class Set;
    }    // namespace mtk
    namespace MSI
    {
        class Model_Solver_Interface;
    }
    namespace fem
    {
        class FEM_Model;
        class IWG;
        class IQI;
        class Field;
        class Field_Interpolator;
        class Geometry_Interpolator;
        class Field_Interpolator_Manager;

        //------------------------------------------------------------------------------
        /**
         * FEM set
         */
        class Set : public MSI::Equation_Set
        {
          private:
            // FEM model
            fem::FEM_Model* mFemModel = nullptr;

            // pointer to the corresponding mesh set
            moris::mtk::Set* mMeshSet = nullptr;

          private:
            // interpolation mesh geometry type
            mtk::Geometry_Type mIPGeometryType = mtk::Geometry_Type::UNDEFINED;

            // integration mesh geometry type
            mtk::Geometry_Type mIGGeometryType = mtk::Geometry_Type::UNDEFINED;

            Vector< Node_Base* > mIPNodes;

            // space interpolation order for IP cells
            mtk::Interpolation_Order mIPSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

            // space interpolation order for IG cells
            mtk::Interpolation_Order mIGSpaceInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

            // time interpolation order for IP cells
            mtk::Interpolation_Order mIPTimeInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

            // space interpolation order for IG cells
            mtk::Interpolation_Order mIGTimeInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

            // field interpolator manager pointers
            Field_Interpolator_Manager* mLeaderFIManager          = nullptr;
            Field_Interpolator_Manager* mFollowerFIManager        = nullptr;
            Field_Interpolator_Manager* mLeaderPreviousFIManager  = nullptr;
            Field_Interpolator_Manager* mLeaderEigenFIManager     = nullptr;
            Field_Interpolator_Manager* mLeaderAdjointFIManager   = nullptr;
            Field_Interpolator_Manager* mFollowerAdjointFIManager = nullptr;

            // number of eigen vectors
            uint mNumEigenVectors = 0;

            // cell of pointers to IWG objects
            Vector< std::shared_ptr< IWG > > mIWGs;
            Vector< std::shared_ptr< IWG > > mRequestedIWGs;

            // cell of pointer to IQI objects
            Vector< std::shared_ptr< IQI > > mIQIs;
            Vector< std::shared_ptr< IQI > > mRequestedIQIs;

            moris::map< std::string, moris_index > mIQINameToIndexMap;

            // cell of pointer to IQI objects for vis
            Vector< std::shared_ptr< IQI > > mRequestedGlobalIQIs;
            Vector< moris_index >            mRequestedGlobalIQIsGlobalIndices;
            Vector< std::shared_ptr< IQI > > mRequestedElementalIQIs;
            Vector< moris_index >            mRequestedElementalIQIsGlobalIndices;
            Vector< std::shared_ptr< IQI > > mRequestedNodalIQIs;
            Vector< moris_index >            mRequestedNodalIQIsGlobalIndices;
            Vector< std::shared_ptr< IQI > > mRequestedFieldIQIs;
            Vector< moris_index >            mRequestedFieldIQIsGlobalIndices;

            // enum for element type
            fem::Element_Type mElementType = fem::Element_Type::UNDEFINED;

            // integration points
            Matrix< DDRMat > mIntegPoints;

            // integration weights
            Matrix< DDRMat > mIntegWeights;

            // map for the dof type
            Matrix< DDSMat > mUniqueDofTypeMap;
            Matrix< DDSMat > mUniqueDvTypeMap;
            Matrix< DDSMat > mUniqueFieldTypeMap;

            // map visualization cell id to position in vector
            Vector< Matrix< DDSMat > > mCellAssemblyMap;    // input: VIS mesh index, VIS cell index || output: position of cell within list of cells on VIS set
            Vector< uint >             mNumIgCellsOnSet;    // input: VIS mesh index || output: number of IG cells in that VIS set

            Vector< Matrix< DDSMat > > mFacetAssemblyMap;    // input: VIS mesh index, IG cell index, side ordinal  || output: position of facet within the output (dbl) side set
            Vector< uint >             mNumFacetsOnSet;      // input: VIS mesh index || output: number of facets in that VIS (dbl) side set

            std::set< Cluster_Measure::ClusterMeasureSpecification > mClusterMeasuresSpecs;

            bool mEvalClusterMeasuresSpecs = true;

            // bool for time continuity
            bool mTimeContinuity = false;

            // bool for time boundary integral
            bool mTimeBoundary = false;

            // bool for analytical/FD SA
            bool mIsAnalyticalFA = true;

            // enum for FD scheme used for FD SA
            fem::FDScheme_Type mFDSchemeForFA = fem::FDScheme_Type::UNDEFINED;

            // real for FD perturbation size
            real mFDPerturbationFA = 0.0;

            // bool for analytical/FD SA
            bool mIsAnalyticalSA = false;

            // enum for FD scheme used for FD SA
            fem::FDScheme_Type mFDSchemeForSA = fem::FDScheme_Type::UNDEFINED;

            // real for FD perturbation size
            real mFDPerturbation = 0.0;

            // enum for perturbation strategy used for FD (FA and SA)
            fem::Perturbation_Type mPerturbationStrategy = fem::Perturbation_Type::RELATIVE;

            Vector< std::shared_ptr< gen::Design_Extraction_Operator > > mIGExtractionOperators;
            map< sint, sint >                                            mVertexMeshIndexToClusterIndexMap;

            Vector< sint >             mAdvGeoAssemblyVector;
            Vector< Matrix< DDRMat > > mAdvGeoWeights;           // Vector of each node with Adv weights
            Matrix< DDRMat >           mCurrrentAdvGeoWeight;    // Adv weight of current node

            friend class MSI::Equation_Object;
            friend class Cluster;
            friend class Element_Bulk;
            friend class Element_Sideset;
            friend class Element_Time_Sideset;
            friend class Element_Double_Sideset;
            friend class Element_Nonconformal_Sideset;
            friend class Element_Time_Continuity;
            friend class Element_Time_Boundary;
            friend class Element;
            friend class Field_Interpolator_Manager;
            friend class Interpolation_Element;

            void gather_requested_IQIs(
                    Vector< std::string > const &     aNames,
                    Vector< std::shared_ptr< IQI > >& aListOfRequestedIQIs,
                    Vector< moris_index >&            aListOfIQIGlobalIndices );

            //------------------------------------------------------------------------------

          public:
            void create_fem_clusters();

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aFemModel a FEM model pointer
             * @param[ in ] aMeshSet  a set from the mesh
             * @param[ in ] aSetInfo  user defined info for set
             * @param[ in ] aIPNodes  cell of node pointers
             */
            Set(
                    fem::FEM_Model*             aFemModel,
                    moris::mtk::Set*            aMeshSet,
                    const fem::Set_User_Info&   aSetInfo,
                    const Vector< Node_Base* >& aIPNodes );

            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Set()
            {
                mIsEmptySet = true;    // FIXME this flag is a hack. find better solution
            }

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~Set() override;

            //------------------------------------------------------------------------------
            /**
             * delete the pointers built on the set
             */
            void delete_pointers();

            //------------------------------------------------------------------------------
            /**
             * initialize the set
             * @param[ in ] aIsStaggered            bool true if computation of a residual
             * @param[ in ] aTimeContinuityOnlyFlag bool true if only time continuity IWGs should be considered
             */
            void initialize_set(
                    const bool                 aIsStaggered            = false,
                    const Time_Continuity_Flag aTimeContinuityOnlyFlag = Time_Continuity_Flag::DEFAULT ) override;

            void update() override;

            //------------------------------------------------------------------------------

            fem::FEM_Model*
            get_fem_model()
            {
                return mFemModel;
            };

            //------------------------------------------------------------------------------
            /**
             * finalize the fem sets
             * create and set the field interpolators
             * @param[ in ] aModelSolverInterface model solver interface pointer
             */
            void finalize( MSI::Model_Solver_Interface* aModelSolverInterface ) override;

            //------------------------------------------------------------------------------
            /**
             * free the memory on the set
             */
            void free_memory() override;

            //------------------------------------------------------------------------------
            /**
             * create integrator
             */
            void create_integrator( MSI::Model_Solver_Interface* aModelSolverInterface );

            //------------------------------------------------------------------------------

            void set_custom_integration_rule( MSI::Model_Solver_Interface* aModelSolverInterface );

            //------------------------------------------------------------------------------

            [[nodiscard]] moris::mtk::Set* get_mesh_set() const
            {
                return mMeshSet;
            }

            //------------------------------------------------------------------------------

            void set_mesh_set( moris::mtk::Set* const aMeshSet )
            {
                mMeshSet = aMeshSet;
            }
            //------------------------------------------------------------------------------
            /**
             * set visualization mesh set
             * @param[ in ] aVisMeshSet a mesh set pointer for visualization
             */
            void set_visualization_set(
                    const uint       aMeshIndex,
                    moris::mtk::Set* aVisMeshSet,
                    const bool       aOnlyPrimaryCells ) override;

            //------------------------------------------------------------------------------

            void
            construct_cell_assembly_map_for_VIS_set(
                    const uint       aMeshIndex,
                    moris::mtk::Set* aVisMeshSet,
                    const bool       aOnlyPrimaryCells );

            //------------------------------------------------------------------------------

            void
            construct_facet_assembly_map_for_VIS_set(
                    const uint       aMeshIndex,
                    moris::mtk::Set* aVisMeshSet );

            //------------------------------------------------------------------------------
            /**
             * get the element on the set
             * @param[ out ] aElementType element type for the set
             */
            fem::Element_Type
            get_element_type() const override
            {
                return mElementType;
            }

            //------------------------------------------------------------------------------
            /**
             * set forward analysis type flag on the set (analytical or finite difference)
             * @param[ in ] aIsAnalyticalFA bool true if analytical forward analysis
             *                                   false if finite difference
             */
            void
            set_is_analytical_forward_analysis( bool aIsAnalyticalFA )
            {
                mIsAnalyticalFA = aIsAnalyticalFA;
            }

            //------------------------------------------------------------------------------
            /**
             * get flag for forward analysis on the set (analytical or finite difference)
             * @param[ out ] mIsAnalyticalFA bool true if analytical forward analysis
             *                                    false if finite difference
             */
            bool
            get_is_analytical_forward_analysis() const
            {
                return mIsAnalyticalFA;
            }

            //------------------------------------------------------------------------------
            /**
             * set FD scheme enum for forward analysis on the set
             * @param[ in ] aFDSchemeForFA enum for FD scheme used for forward analysis
             */
            void
            set_finite_difference_scheme_for_forward_analysis(
                    enum fem::FDScheme_Type aFDSchemeForFA )
            {
                mFDSchemeForFA = aFDSchemeForFA;
            }

            //------------------------------------------------------------------------------
            /**
             * get enum for FD scheme for forward analysis on the set
             * @param[ out ] mFDSchemeForFA enum for FD scheme used for forward analysis
             */
            enum fem::FDScheme_Type
            get_finite_difference_scheme_for_forward_analysis() const
            {
                return mFDSchemeForFA;
            }

            //------------------------------------------------------------------------------
            /**
             * set perturbation size for finite difference for forward analysis
             * @param[ in ] aFDPerturbationFA perturbation size
             */
            void
            set_finite_difference_perturbation_size_forward( real aFDPerturbationFA )
            {
                mFDPerturbationFA = aFDPerturbationFA;
            }

            //------------------------------------------------------------------------------
            /**
             * get perturbation size for finite difference for forward analysis
             * @param[ out ] mFDPerturbation perturbation size
             */
            real
            get_finite_difference_perturbation_size_forward()
            {
                return mFDPerturbationFA;
            }

            //------------------------------------------------------------------------------
            /**
             * set sensitivity analysis type flag on the set (analytical or finite difference)
             * @param[ in ] aIsAnalyticalSA bool true if analytical sensitivity analysis
             *                                   false if finite difference
             */
            void
            set_is_analytical_sensitivity_analysis( bool aIsAnalyticalSA )
            {
                mIsAnalyticalSA = aIsAnalyticalSA;
            }

            //------------------------------------------------------------------------------
            /**
             * get flag for sensitivity analysis on the set (analytical or finite difference)
             * @param[ out ] mIsAnalyticalSA bool true if analytical sensitivity analysis
             *                                    false if finite difference
             */
            bool
            is_analytical_sensitivity_analysis() const
            {
                return mIsAnalyticalSA;
            }

            //------------------------------------------------------------------------------
            /**
             * set FD scheme enum for sensitivity analysis on the set
             * @param[ in ] aFDSchemeForSA enum for FD scheme used for
             */
            void
            set_finite_difference_scheme_for_sensitivity_analysis(
                    enum fem::FDScheme_Type aFDSchemeForSA )
            {
                mFDSchemeForSA = aFDSchemeForSA;
            }

            //------------------------------------------------------------------------------
            /**
             * get enum for FD scheme for sensitivity analysis on the set
             * @param[ out ] mFDSchemeForSA enum for FD scheme used for
             */
            enum fem::FDScheme_Type
            get_finite_difference_scheme_for_sensitivity_analysis() const
            {
                return mFDSchemeForSA;
            }

            //------------------------------------------------------------------------------
            /**
             * set perturbation size for finite difference
             * @param[ in ] aFDPerturbation perturbation size
             */
            void
            set_finite_difference_perturbation_size( real aFDPerturbation )
            {
                mFDPerturbation = aFDPerturbation;
            }

            //------------------------------------------------------------------------------
            /**
             * get perturbation size for finite difference
             * @param[ out ] mFDPerturbation perturbation size
             */
            real
            get_finite_difference_perturbation_size()
            {
                return mFDPerturbation;
            }

            //------------------------------------------------------------------------------
            /**
             * set perturbation strategy enum for sensitivity analysis on the set
             * @param[ in ] aPerturbationStrategy enum for perturbation strategy used for
             * FA and SA if computed by FD
             */
            void
            set_perturbation_strategy(
                    enum fem::Perturbation_Type aPerturbationStrategy )
            {
                mPerturbationStrategy = aPerturbationStrategy;
            }

            //------------------------------------------------------------------------------
            /**
             * get perturbation strategy enum for sensitivity analysis on the set
             * @param[ out ] aPerturbationStrategy enum for perturbation strategy used for
             * FA and SA if computed by FD
             */
            enum fem::Perturbation_Type
            get_perturbation_strategy() const
            {
                return mPerturbationStrategy;
            }

            //------------------------------------------------------------------------------
            /**
             * get the clusters on the set
             * @param[ out ] aClusters cell of mesh cluster pointers
             */
            // Vector< mtk::Cluster const* > const&
            // get_clusters_on_set() const
            // {
            //     return mMeshSet->get_clusters_on_set();
            // }

            //------------------------------------------------------------------------------
            /**
             * get bool for time continuity
             * param[ out ] mTimeContinuity bool true if time continuity
             */
            bool
            get_time_continuity() const
            {
                return mTimeContinuity;
            }

            //------------------------------------------------------------------------------
            /**
             * create a unique dof type list for the solver
             * Vector< MSI::Dof_Type >, no group of dof type
             * one for both leader and follower
             */
            void create_unique_dof_and_dv_type_lists();

            //------------------------------------------------------------------------------
            /**
             * create a unique group of dof type list for the set
             * Vector< Vector< MSI::Dof_Type > > list of groups of dof type
             * one for the leader, one for the follower
             */
            void create_dof_and_dv_type_lists();

            //------------------------------------------------------------------------------
            /**
             * create a map of the dof type for the set
             * one for the leader, one for the follower
             */
            void create_dof_and_dv_type_maps();

            //------------------------------------------------------------------------------
            /**
             * create field interpolator managers for the set
             * @param[ in ] aModelSolverInterface model solver interface
             * ( only used to set the time levels )
             */
            void create_field_interpolator_managers( MSI::Model_Solver_Interface* aModelSolverInterface );

            //------------------------------------------------------------------------------
            /**
             * get IWGs
             * @param[ out ] aIWGs cell of IWG pointers
             */
            Vector< std::shared_ptr< IWG > >&
            get_IWGs()
            {
                return mIWGs;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of IWGs
             */
            uint
            get_number_of_IWGs()
            {
                return mIWGs.size();
            }

            //------------------------------------------------------------------------------
            /**
             * get requested IWGs
             * @param[ out ] mRequestedIWGs cell of requested IWG pointers
             */
            const Vector< std::shared_ptr< IWG > >&
            get_requested_IWGs()
            {
                return mRequestedIWGs;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of requested IWGs
             */
            uint
            get_number_of_requested_IWGs()
            {
                return mRequestedIWGs.size();
            }

            //------------------------------------------------------------------------------
            /**
             * building an IQI name to set local index map
             */
            void create_IQI_map();

            //------------------------------------------------------------------------------

            std::unordered_map< moris_index, Vector< real > >
            get_nodal_displacements( const std::unordered_set< moris_index >& aRequestedNodes ) override;

            //------------------------------------------------------------------------------
            /**
             * building an IQI name to set local index map
             */
            std::string get_set_name() override;

            //------------------------------------------------------------------------------
            /**
             * create local requested IQIs for sensitivity analysis
             */
            void create_requested_IQI_list();

            //------------------------------------------------------------------------------
            /**
             * get local requested IQIs for sensitivity analysis
             * @param[ out ] mRequestedIQIs cell of IQIs pointers
             */
            const Vector< std::shared_ptr< IQI > >&
            get_requested_IQIs()
            {
                return mRequestedIQIs;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of local requested IQIs for sensitivity analysis
             * @param[ out ] uint number of local requested IQIs for SA
             */
            uint
            get_number_of_requested_IQIs() override
            {
                return mRequestedIQIs.size();
            }

            //------------------------------------------------------------------------------
            /**
             * set the field interpolator managers for the IWGs
             */
            void set_IWG_field_interpolator_managers();

            //------------------------------------------------------------------------------
            /**
             * set the field interpolator managers for the IQIs
             */
            void set_IQI_field_interpolator_managers();

            //------------------------------------------------------------------------------
            /*
             * set the cluster for the IWG stabilization parameters
             * associated with this set
             * @param[ in ] aFemCluster fem cluster pointer
             */
            void set_IWG_cluster_for_stabilization_parameters( fem::Cluster* aFemCluster );

            //------------------------------------------------------------------------------
            /*
             * set the cluster for the IQI stabilization parameters
             * associated with this set
             * @param[ in ] aFemCluster fem cluster pointer
             */
            void set_IQI_cluster_for_stabilization_parameters( fem::Cluster* aFemCluster );

            //------------------------------------------------------------------------------
            /*
             * build cluster measure specification list and map required on set
             */
            void build_cluster_measure_specifications();

            //------------------------------------------------------------------------------
            /*
             * get cluster measures required on set
             * return cell of tuple with cluster measure specifications
             */
            Vector< std::tuple<
                    fem::Measure_Type,
                    mtk::Primary_Void,
                    mtk::Leader_Follower > >&
            get_cluster_measure_tuples();

            //------------------------------------------------------------------------------
            /*
             * get cluster measures required on set
             * return cell of tuple with cluster measure specifications
             */
            std::set< Cluster_Measure::ClusterMeasureSpecification > const &
            get_cluster_measure_specifications();

            //------------------------------------------------------------------------------
            /**
             * create the dof assembly map for the residual/rows
             */
            void create_residual_dof_assembly_map();

            //------------------------------------------------------------------------------
            /**
             * create the dof assembly map for the jacobian/cols
             * @param[ in ] aIsStaggered bool true if residual evaluation
             */
            void create_dof_assembly_map( const bool aIsStaggered );

            //------------------------------------------------------------------------------
            /**
             * create the dof assembly map for the jacobian/cols
             */
            void create_jacobian_dof_assembly_map();

            //------------------------------------------------------------------------------
            /**
             * create the dof assembly map
             * for the off-diagonal requested jacobian/cols
             * for R = R_0 - A_{01} x_{1}
             */
            void create_staggered_jacobian_dof_assembly_map();

            //------------------------------------------------------------------------------
            /**
             * create the material pdv assembly map
             */
            void create_mat_pdv_assembly_map();

            //------------------------------------------------------------------------------
            /**
             * create the geometry pdv assembly map
             * @param[ in ] aFemCluster a fem cluster pointer
             */
            void create_geo_pdv_assembly_map( const std::shared_ptr< fem::Cluster >& aFemCluster );

            //------------------------------------------------------------------------------
            /**
             * create the adv weights and assembly vector at integration point
             */
            void create_geo_adv_assembly_data( const Matrix< IndexMat >& tVertexMeshIndices );

            //------------------------------------------------------------------------------
            /**
             * create a list of IWGs requested by the solver
             */
            void create_requested_IWG_list( const fem::Time_Continuity_Flag aTimeContinuityOnlyFlag = moris::fem::Time_Continuity_Flag::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * create a dof type list for the list of IWGs requested by the solver
             * @param[ in ] aIsStaggered bool true if residual evaluation
             */
            void build_requested_IWG_dof_type_list( const bool aIsStaggered );

            /**
             * Create a dof type list for the list of IQIs being requested.
             */
            void build_requested_IQI_dof_type_list();

            //------------------------------------------------------------------------------
            /**
             * get the IP geometry type
             * @param[ out ] aGeometryType
             */
            mtk::Geometry_Type
            get_IP_geometry_type()
            {
                return mIPGeometryType;
            }

            //------------------------------------------------------------------------------
            /**
             * get the IG geometry type
             * @param[ out ] aGeometryType
             */
            mtk::Geometry_Type
            get_IG_geometry_type()
            {
                return mIGGeometryType;
            }

            //------------------------------------------------------------------------------
            /**
             * FIXME we should not be able to ask this question
             * as the interpolator order is not unique on the set
             * get the IG space interpolation order
             */
            mtk::Interpolation_Order
            get_IG_space_interpolation_order()
            {
                return mIGSpaceInterpolationOrder;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of integration points
             * @param[ out ] aNumIntegPoints number of integration points
             */
            uint
            get_number_of_integration_points()
            {
                return mIntegWeights.numel();
            }

            //------------------------------------------------------------------------------
            /**
             * get the integration points
             * @param[ out ] aIntegPoints integration points
             */
            const Matrix< DDRMat >&
            get_integration_points()
            {
                return mIntegPoints;
            }

            //------------------------------------------------------------------------------
            /**
             * get the integration weights
             * @param[ out ] aIntegWeights integration weights
             */
            const Matrix< DDRMat >&
            get_integration_weights()
            {
                return mIntegWeights;
            }

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator manager
             * @param[ in ]  aIsLeader an enum for leader or follower
             * @param[ out ] mFIManger a field interpolator manager pointer
             */
            Field_Interpolator_Manager* get_field_interpolator_manager(
                    mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator manager for previous time step
             * @param[ in ]  aIsLeader an enum for leader or follower
             * @param[ out ] mFIManger a field interpolator manager pointer
             */
            Field_Interpolator_Manager* get_field_interpolator_manager_previous_time(
                    mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator manager for eigen vectors
             * @param[ in ]  aIsLeader an enum for leader or follower
             * @param[ out ] mFIManger a field interpolator manager pointer
             */
            Field_Interpolator_Manager* get_field_interpolator_manager_eigen_vectors(
                    mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * get the field interpolator manager for adjoint vectors
             * @param[ in ]  aIsLeader an enum for leader or follower
             * @param[ out ] mFIManger a field interpolator manager pointer
             */
            Field_Interpolator_Manager* get_field_interpolator_manager_adjoint_vectors(
                    mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * auto detect full integration scheme
             * @param[ in ]  aGeometryType       a geometry type
             * @param[ in ]  aInterpolationOrder an interpolation order
             * @param[ out ] aIntegrationOrder   an integration order
             */
            mtk::Integration_Order get_auto_integration_order(
                    const fem::Element_Type        aElementType,
                    const mtk::Geometry_Type       aGeometryType,
                    const mtk::Interpolation_Order aInterpolationOrder );

            //------------------------------------------------------------------------------
            /**
             * auto detect space interpolation scheme
             * @param[ in ] aNumVertices  number of vertices for the considered geometry type
             * @param[ in ] aGeometryType geometry type
             * FIXME: works for Lagrange only
             */
            mtk::Interpolation_Order get_auto_interpolation_order(
                    const moris::uint        aNumVertices,
                    const mtk::Geometry_Type aGeometryType );

            //------------------------------------------------------------------------------
            /**
             * auto detect time interpolation scheme
             * @param[ in ] aNumVertices number of vertices for a line
             */
            mtk::Interpolation_Type get_auto_time_interpolation_type(
                    const moris::uint aNumVertices );

            //------------------------------------------------------------------------------
            /**
             * set size and reset values for Jacobian
             */
            void initialize_mJacobian();

            //------------------------------------------------------------------------------
            /**
             * set size and reset values for residual
             */
            void initialize_mResidual( sint aNumRHS = -1 );

            //------------------------------------------------------------------------------
            /**
             * set size and reset values for QI
             */
            void initialize_mQI();

            //----------------------------------------------------------------------
            /**
             * set size and reset values for dQIdp
             */
            void initialize_mdQIdpMat();

            //----------------------------------------------------------------------
            /**
             * set size and reset values for dQIdp
             */
            void initialize_mdQIdpGeo( const std::shared_ptr< fem::Cluster >& aFemCluster );

            //------------------------------------------------------------------------------
            /**
             * set size and reset values for dRdp mat
             */
            void initialize_mdRdpMat();

            //------------------------------------------------------------------------------
            /**
             * set size and reset values for dRdp geo
             */
            void initialize_mdRdpGeo( const std::shared_ptr< fem::Cluster >& aFemCluster );

            //------------------------------------------------------------------------------
            /*
             * get index from unique dof type map
             *@param[ in ] aDofType a dof type enum
             */
            moris::sint
            get_index_from_unique_dof_type_map( enum MSI::Dof_Type aDofType )
            {
                return mUniqueDofTypeMap( static_cast< int >( aDofType ), 0 );
            }

            //------------------------------------------------------------------------------
            /*
             * get index from unique dv type map
             *@param[ in ] aDofType a dof type enum
             */
            moris::sint
            get_index_from_unique_dv_type_map( enum gen::PDV_Type aDvType )
            {
                return mUniqueDvTypeMap( static_cast< int >( aDvType ), 0 );
            }

            //------------------------------------------------------------------------------
            /*
             * get index from unique field type map
             *@param[ in ] aDofType a dof type enum
             */
            moris::sint
            get_index_from_unique_field_type_map( mtk::Field_Type aFieldType )
            {
                return mUniqueFieldTypeMap( static_cast< int >( aFieldType ), 0 );
            }

            //------------------------------------------------------------------------------
            /*
             * create requested dv assembly map
             */
            void create_requested_dv_assembly_map();

            //------------------------------------------------------------------------------
            /*
             * create dof and dv type map unique
             */
            void create_unique_dof_dv_and_field_type_maps();

            //------------------------------------------------------------------------------
            /**
             * compute a quantity of interest nodal
             * @param[ in ] aVisMeshIndex     vis mesh index to define mesh
             *                                on which values are evaluated
             * @param[ in ] aNodalFieldValues matrix to fill with nodal values
             * @param[ in ] aQINames          list of IQI names to be evaluated
             */
            void compute_quantity_of_interest_nodal(
                    const uint                   aVisMeshIndex,
                    Matrix< DDRMat >*            aNodalFieldValues,
                    const Vector< std::string >& aQINames ) override;

            //------------------------------------------------------------------------------
            /**
             * get list of requested IQIs for nodal evaluation for visualization
             */
            const Vector< std::shared_ptr< fem::IQI > >& get_requested_nodal_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get number of requested IQIs for nodal evaluation for visualization
             */
            uint get_number_of_requested_nodal_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get global indices for the list of requested IQIs
             * for nodal evaluation for visualization
             */
            const Vector< moris_index >& get_requested_nodal_IQIs_global_indices_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * compute a quantity of interest global
             * @param[ in ] aVisMeshIndex      vis mesh index to define mesh
             *                                 on which values are evaluated
             * @param[ in ] aGlobalFieldValues matrix to fill with global values
             * @param[ in ] aQINames           list of IQI names to be evaluated
             */
            void compute_quantity_of_interest_global(
                    const uint                   aVisMeshIndex,
                    Matrix< DDRMat >*            aGlobalFieldValues,
                    const Vector< std::string >& aQINames ) override;

            //------------------------------------------------------------------------------
            /**
             * get list of requested IQIs for global evaluation for visualization
             */
            const Vector< std::shared_ptr< fem::IQI > >& get_requested_global_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get number of requested IQIs for global evaluation for visualization
             */
            uint get_number_of_requested_global_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get global indices for the list of requested IQIs
             * for global evaluation for visualization
             */
            const Vector< moris_index >& get_requested_global_IQIs_global_indices_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * compute a quantity of interest elemental
             * @param[ in ] aVisMeshIndex         vis mesh index to define mesh
             *                                    on which values are evaluated
             * @param[ in ] aElementalFieldValues matrix to fill with elemental values
             * @param[ in ] aQINames              list of IQI names to be evaluated
             * @param[ in ] aOutputAverageValue   whether the value is an average on the element, or the integrated quantity on the element
             */
            void compute_quantity_of_interest_elemental(
                    const uint                   aVisMeshIndex,
                    Matrix< DDRMat >*            aElementalFieldValues,
                    const Vector< std::string >& aQINames,
                    const bool                   aOutputAverageValue = true ) override;

            //------------------------------------------------------------------------------
            /**
             * get list of requested IQIs for elemental evaluation for visualization
             */
            const Vector< std::shared_ptr< fem::IQI > >& get_requested_elemental_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get number of requested IQIs for elemental evaluation for visualization
             */
            uint get_number_of_requested_elemental_IQIs_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get global indices for the list of requested IQIs
             * for elemental evaluation for visualization
             */
            const Vector< moris_index >& get_requested_elemental_IQIs_global_indices_for_visualization();

            //------------------------------------------------------------------------------
            /**
             * get list of requested IQIs for field evaluation
             */
            const Vector< std::shared_ptr< fem::IQI > >& get_requested_field_IQIs();

            //------------------------------------------------------------------------------
            /**
             * get number of requested IQIs for field evaluation
             */
            uint get_number_of_requested_field_IQIs();

            //------------------------------------------------------------------------------
            /**
             * get global indices for the list of requested IQIs
             * for field evaluation
             */
            const Vector< moris_index >& get_requested_field_IQIs_global_indices();

            //------------------------------------------------------------------------------
            /**
             * determine set type from mtk set type
             */
            void determine_set_type();

            //------------------------------------------------------------------------------
            /**
             * set set type
             * only for debug and unit test
             */
            void set_set_type( fem::Element_Type aElementType );

            //------------------------------------------------------------------------------
            /**
             * get ig unique dv types for set from design variable interface
             * @param[ in ] aGeoPdvType list of ig unique pdv types on set
             */
            void get_ig_unique_dv_types_for_set( Vector< enum gen::PDV_Type >& aGeoPdvType );

            //------------------------------------------------------------------------------

            /**
             * get ip dv types for set from design variable interface
             * @param[ in ] aMatPdvType list of group of ip pdv types on set
             * @param[ in ] aIsLeader determine the leader / follower side, only for dbl sided set, is leader by default
             */
            void
            get_ip_dv_types_for_set(
                    Vector< Vector< enum gen::PDV_Type > >& aMatPdvType,
                    mtk::Leader_Follower                    aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------

            void populate_fields(
                    Vector< std::shared_ptr< fem::Field > >& aFieldToPopulate,
                    Vector< std::string > const &            aFieldIQINames ) override;

            //------------------------------------------------------------------------------
            /**
             * get the weights of design extraction operator of the cluster
             *
             * @return a vector of matrices of design extraction operator weights
             */
            const Vector< Matrix< DDRMat > >& get_adv_geo_weights()
            {
                return mAdvGeoWeights;
            }

            //------------------------------------------------------------------------------

            void set_geo_weights_for_cluster_node_index( uint aIndex )
            {
                MORIS_ASSERT( aIndex < mAdvGeoWeights.size(),
                        "Set::set_index_for_geo_weights - Index out of range" );

                mCurrrentAdvGeoWeight = mAdvGeoWeights( aIndex );
            }

            //------------------------------------------------------------------------------

            /**
             * get the local index in cluster of an optimization variable defined by its id
             *
             * @param[ in ] aOptimizationVariableId optimization variable id
             * @return local index
             */
            sint
            get_current_adv_geo_index( uint aOptimizationVariableId )
            {
                return mIGAdvIds.find( aOptimizationVariableId );
            }

            //------------------------------------------------------------------------------

            /**
             * get the current adv geo weight
             *
             * @return current adv geo weight
             */
            const Matrix< DDRMat >& get_current_adv_geo_weight()
            {
                return mCurrrentAdvGeoWeight;
            }
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_SET_HPP_ */
