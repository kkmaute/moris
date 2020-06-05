#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_Geometry_Engine_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_Geometry_Engine_HPP_

// GEN
#include "cl_GEN_Interpolaton.hpp"
#include "cl_GEN_Pending_Node.hpp"
#include "cl_GEN_Phase_Table.hpp"

#include "cl_GEN_Geometry_Object.hpp"
#include "cl_GEN_Geometry_Object_Manager.hpp"

#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

// MTK
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_Mesh_Enums.hpp"

// MRS
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"
#include "fn_trans.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    namespace hmr
    {
        class HMR;
        class Mesh;
    }
    namespace MSI
    {
        class Design_Variable_Interface;
    }

    //------------------------------------------------------------------------------------------------------------------

    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        /**
         * $\frac{\partial{\phi_A}}{\partial{p}}$ (change in phi with respect to a design variable
         * See for more detailed description of this function:
         */
        inline
        void compute_dx_dp_with_linear_basis( moris::Matrix< moris::DDRMat >  & aDPhiADp,
                                              moris::Matrix< moris::DDRMat >  & aDPhiBDp,
                                              moris::Matrix< moris::DDRMat >  & aEdgeCoordinates,
                                              moris::Matrix< moris::DDRMat >  & aEdgeNodePhi,
                                              moris::Matrix< moris::DDRMat >  & aDxDp )
        {
            moris::real const & tPhiA = aEdgeNodePhi(0,0);
            moris::real const & tPhiB = aEdgeNodePhi(1,0);

             // Initialize
            moris::Matrix< moris::DDRMat > tXa = trans(aEdgeCoordinates.get_row(0));
            moris::Matrix< moris::DDRMat > tXb = trans(aEdgeCoordinates.get_row(1));

            // Compute $\frac{\partial x_{\Gamma}}{\partial \phi}$
            moris::DDRMat tDxgammaDphiA = -(tPhiB)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());
            moris::DDRMat tDxgammaDphiB =  (tPhiA)/std::pow((tPhiA-tPhiB),2)*(tXb.matrix_data()-tXa.matrix_data());
            moris::Matrix< moris::DDRMat > tDxgDphiAMat(tDxgammaDphiA);
            moris::Matrix< moris::DDRMat > tDxgDphiBMat(tDxgammaDphiB);

            // Compute dx/dp
            moris::DDRMat tDxDp = tDxgDphiAMat * aDPhiADp + tDxgDphiBMat * aDPhiBDp;
            aDxDp = moris::Matrix< moris::DDRMat >(tDxDp);
        }

        //--------------------------------------------------------------------------------------------------------------
        // GEOMETRY ENGINE
        //--------------------------------------------------------------------------------------------------------------

        class Geometry_Engine
        {
        public:
            moris::uint mSpatialDim;
            moris::real mThresholdValue;
            moris::real mPerturbationValue;
            bool mComputeDxDp = false; // FIXME change this?

        private:

            // HMR refinements
            uint mNumRefinements;

            // ADVs/IQIs
            Matrix<DDRMat> mADVs;
            Matrix<DDRMat> mLowerBounds;
            Matrix<DDRMat> mUpperBounds;
            Cell<std::string> mRequestedIQIs;

            // Geometry
            moris::size_t mActiveGeometryIndex = 0;
            Cell<std::shared_ptr<Geometry>> mGeometry;

            // Property
            Cell<std::shared_ptr<Property>> mProperties;
            Cell<ParameterList> mPropertyParameterLists;

            // Contains all the geometry objects
            Geometry_Object_Manager mGeometryObjectManager;

            // Contains all the pdv hosts
            Pdv_Host_Manager mPdvHostManager;

            // Phase Table
            Phase_Table mPhaseTable;

            // Node Entity Phase Vals - only analytic phase values are stored here to prevent duplicate storage of discrete geometries
            moris::Matrix< moris::DDRMat > mNodePhaseVals;

            // Mesh
            mtk::Mesh_Manager* mMeshManager;
            moris::Cell< std::shared_ptr< moris::hmr::HMR > > mHMRPerformer;
            moris::Cell< std::shared_ptr< moris::hmr::Mesh > > mMesh_HMR; //FIXME needs to be more general to only have a mesh manager as this member

            bool mTypesSet      = false;
            moris::Cell< moris::moris_index > mIntegNodeIndices;

            Matrix<IndexMat> mInterfaceNodeIndices;

        public:

            /**
             * Constructor using cell of cell of parameter lists
             *
             * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
             * @param aLibrary Library used for pulling user-defined functions
             */
            Geometry_Engine(moris::Cell<moris::Cell<ParameterList>> aParameterLists,
                            std::shared_ptr<moris::Library_IO> aLibrary = nullptr);

            /**
             * Constructor using externally created geometries and phase table
             *
             * @param[ in ] aGeometry cell of shared Geometry pointers
             * @param[ in ] aPhaseTable phase table
             * @param[ in ] aSpatialDim spatial dimensions
             */
            Geometry_Engine(Cell< std::shared_ptr<Geometry> >   aGeometry,
                                Phase_Table                     aPhaseTable,
                                uint                            aSpatialDim = 3,
                                real                            aThresholdValue = 0.0,
                                real                            aPerturbationValue = 1E-6);

            /**
             * Destructor
             */
            ~Geometry_Engine();

            /**
             * Sets new advs for the geometry engine
             *
             * @param aNewADVs vector of new advs to use
             */
            void set_advs(Matrix<DDRMat> aNewADVs);

            /**
             * Gets the advs from the geometry engine
             *
             * @return vector of advs
             */
            Matrix<DDRMat>& get_advs();

            /**
             * Gets the lower bounds from the geometry engine
             *
             * @return vector of lower bounds
             */
            Matrix<DDRMat>& get_lower_bounds();

            /**
             * Gets the upper bounds from the geometry engine
             *
             * @return vector of upper bounds
             */
            Matrix<DDRMat>& get_upper_bounds();

            /**
             * Lets MDL know about the stored requested IQIs through the PDV host manager
             */
            void communicate_requested_IQIs();
            void communicate_requested_IQIs(Cell<std::string> aIQINames);

            /**
             * Gets the sensitivities of the critieria with respect to the advs
             *
             * @return Matrix of sensitivities
             */
            Matrix<DDRMat> get_dcriteria_dadv();

            /**
             * Gets the design variable interface from the geometry engine
             *
             * @return member pdv host manager pointer
             */
            MSI::Design_Variable_Interface* get_design_variable_interface();

            /**
             * @brief Initial allocation of geometry objects,
             * this creates a geometry object for each node coordinate.
             * In this case, aNodeCoords needs to be ordered by proc indices.
             * @param[ in ] aNumNodes number of nodes
             */
            void initialize_geometry_objects_for_background_mesh_nodes( moris::size_t const & aNumNodes );

            /**
             * ???
             * @param[ in ] aNodeCoords node coordinates
             */
            void initialize_geometry_object_phase_values( moris::Matrix< moris::DDRMat > const & aNodeCoords );

            /**
             * @brief Creates a geometry object association for pending nodes
             * These nodes have node indices and parent information
             */
            void associate_new_nodes_with_geometry_object( moris::Cell< Pending_Node > & aNewNodes,
                                                           bool                          aInterfaceNodes );

            /**
             * create new node geometry objects
             * @param[ in ] aNodeCoords node coordinates
             */
            void create_new_node_geometry_objects(Cell< moris_index >  const & aNewNodeIndices,
                                                  bool                         aStoreParentTopo,
                                                  Cell<xtk::Topology*> const & aParentTopo,
                                                  Cell<Matrix<DDRMat>> const & aParamCoordRelativeToParent,
                                                  Cell<Matrix<DDRMat>> const & aGlobalNodeCoord);

            /**
             * @brief Links new nodes with an existing geometry object. This is used for unzipped interfaces
             * where more than one node is at the same location
             * @param[in] aNodesIndicesWithGeomObj - Node indices which already have a geometry object
             * @param[in] aNodesIndicesToLink - Node indices to link to the corresponding nodes in aNodesIndicesWithGeomObj
             */
            void link_new_nodes_to_existing_geometry_objects( Matrix< IndexMat > const & aNodesIndicesWithGeomObj,
                                                              Matrix< IndexMat > const & aNodesIndicesToLink );

            /**
             * @brief is_intersected checks to see if an entity provided to it intersects a geometry field. Intersects in this context
             * means a geometry crosses a certain threshold (typically 0). For levelset fields, this can be thought of as a phase change
             *
             * @param[in] aNodeCoords       - Node coordinate
             * @param[in] aNodeToEntityConn - Connectivity between nodes and parent entity
             * @param[in] aCheckType        - Specifies what type of intersection check is to be performed
             *                                   0 - No information on interface required
             *                                   1 - information on interface required
             */
            void is_intersected( moris::Matrix< moris::DDRMat > const &   aNodeCoords,
                                 moris::Matrix< moris::IndexMat > const & aNodetoEntityConn,
                                 moris::size_t                            aCheckType,
                                 Cell<GEN_Geometry_Object> &              aGeometryObjects );

            /*!
             * @brief Computes the interface sensitivity of the provided node indices. After this call,
             * the sensitivity information of these interface nodes can be accessed through the interface
             * nodes respective geometry object.
             * @param[in] aInterfaceNodeIndices - Interface Node Indices (should be interface nodes wrt geometry index provided)
             * @param[in] aNodeCoords -  Node coordinates with location corresponding to indices of aIntefaceNodeIndices.
             * @param[in] aGeomIndex - Geometry Index
             * @param[in] aGlbCoord  - bool to calculate the global coordinate of the intersection point
             */
            void compute_interface_sensitivity( Matrix< IndexMat > const & aInterfaceNodeIndices,
                                                Matrix< DDRMat >   const & aNodeCoords,
                                                moris_index                aGeomIndex,
                                                bool               const   aGlbCoord = false );

            /**
             * Computes the intersection of an isocountour with an entity and returning the local coordinate relative to the parent
             * and the global coordinate if needed
             */
            void get_intersection_location( moris::real const &                      aIsocontourThreshold,
                                            moris::real const &                      aPerturbationThreshold,
                                            moris::Matrix< moris::DDRMat > const &   aGlobalNodeCoordinates,
                                            moris::Matrix< moris::DDRMat > const &   aEntityNodeVars,
                                            moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                            moris::Matrix< moris::DDRMat > &         aIntersectionLocalCoordinates,
                                            moris::Matrix< moris::DDRMat > &         aIntersectionGlobalCoordinates,
                                            bool                                     aCheckLocalCoordinate = true,
                                            bool                                     aComputeGlobalCoordinate = false );

            /**
             * Computes dx/dp using finite differencing
             */
            void compute_dx_dp_finite_difference( moris::real                      const & aPerturbationVal,
                                                  moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                                  moris::Matrix< moris::DDRMat >   const & aEntityNodeCoordinates,
                                                  moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                                  moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                                  moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                                  moris::Matrix< moris::DDRMat >         & aDxDp );

            /**
             * Computes dx/dp of an intersection
             */
            void compute_dx_dp_for_an_intersection( moris::Matrix< moris::IndexMat > const & aEntityNodeIndices,
                                                    moris::Matrix< moris::DDRMat >   const & aGlobalNodeCoordinates,
                                                    moris::Matrix< moris::DDRMat >   const & aIntersectionLclCoordinate,
                                                    moris::Matrix< moris::DDRMat >         & aEntityNodeVars,
                                                    moris::Matrix< moris::DDRMat >         & aDxDp,
                                                    moris::Matrix< moris::IndexMat >       & aADVIndices );

            /**
             * @brief Get the total number of phases in the phase table
             */
            moris::size_t get_num_phases();

            /**
             * @brief Get the 0 or 1 value associated with a given phase and geometry index
             */
            moris::moris_index
            get_phase_sign_of_given_phase_and_geometry( moris::moris_index aPhaseIndex,
                                                        moris::moris_index aGeometryIndex );

            /**
             * @brief Get phase value for a given node and geometry index
             */
            moris::real
            get_entity_phase_val( moris::size_t const & aNodeIndex,
                                  moris::size_t const & aGeomIndex );

            /**
             * @brief Get dxdp for a node
             */
            moris::Matrix< moris::DDRMat > const &
            get_node_dx_dp(moris::size_t const & aNodeIndex) const;

            /**
             * @brief get adv indices for a node
             */
            moris::Matrix< moris::IndexMat > const &
            get_node_adv_indices( moris::size_t const & aNodeIndex ) const;

            /**
             * @brief For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
             */
            void get_phase_index( moris::Matrix< moris::DDSTMat > const & aNodeIndex,
                                  moris::Matrix< moris::DDSTMat > & aNodePhaseIndex );

            /**
              * @brief For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
              */
             void get_phase_index( moris::moris_index const & aNodeIndex,
                                   moris::size_t & aNodePhaseIndex );

            /**
             * @brief Provided the inside and out phase values for an entity, return the phase index
             */
            moris::moris_index
            get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff);

            /**
             * @brief Returns whether a node is inside or outside wrt to a given geometry index
             */
            moris::size_t
            get_node_phase_index_wrt_a_geometry(moris::size_t aNodeIndex,
                                                moris::size_t aGeometryIndex);

            /**
             * @brief Returns the number of geometries
             */
            moris::size_t get_num_geometries();

            /**
             * @brief Returns the number of phases
             */
            moris::size_t get_num_bulk_phase();

            /**
             * @brief Returns the active geometry index
             */
            moris::size_t get_active_geometry_index();

            /**
             * @brief Advance the active geometry index
             */
            void advance_geometry_index();

            /**
             * this function need to be deleted as they are not used in the current PDV_Type interface implementation !!!
             */
            moris::Matrix< moris::IndexMat > get_node_adv_indices_analytic();

            /**
             * @brief Returns the ADV indices of the provided nodes
             * this function need to be deleted as they are not used in the current PDV_Type interface implementation !!!
             */
            moris::Matrix< moris::IndexMat > get_node_adv_indices_discrete
            ( moris::Matrix< moris::IndexMat > const & aEntityNodes );

            /**
             * Returns the number of advs
             */
            moris::size_t get_num_design_variables();

            /**
             * Register an MTK mesh pair to be used for later computation(s)
             */
            void register_mesh( mtk::Mesh_Manager* aMesh );

            /**
             * Register an HMR mesh
             *
             * @warning will be removed in a future update, GE will only be able to register an mtk mesh pair
             */
            moris_index register_mesh( std::shared_ptr< moris::hmr::Mesh > aMesh ); //FIXME: this needs to be deleted and the GE should only be able to register an mtk mesh pair

            /**
             * Register an HMR mesh
             *
             * @warning will be removed in a future update, GE will only be able to register an mtk mesh pair
             */
            moris_index register_mesh( std::shared_ptr< hmr::HMR > aHMR );

            /**
             * Allows GE to become a performer
             */
            void perform( );

            /**
             * Performs refinement on an HMR mesh
             */
            void perform_refinement( );

            /**
             * Fills a cell of MORIS matrices with the level-set values corresponding to each geometry
             */
            void get_field_values_for_all_geometries( moris::Cell< Matrix< DDRMat > > & aAllFieldVals,
                                                      const moris_index                 aWhichMesh = 0 );

            /**
             * @brief assign the pdv type and property for each pdv host in a given set
             */
            void assign_ip_hosts_by_set_name( std::string                     aSetName,
                                              std::shared_ptr< Property > aPropertyPointer,
                                              PDV_Type                     aPdvType,
                                              moris_index                     aWhichMesh = 0 );

            /**
             * @brief assign the pdv type and property for each pdv host in a given set
             */
            void assign_ip_hosts_by_set_index( moris_index                     aSetIndex,
                                               std::shared_ptr< Property > aPropertyPointer,
                                               PDV_Type                     aPdvType,
                                               moris_index                     aWhichMesh = 0 );

            /**
             * Create PDV_Type hosts with the specified PDV_Type types on the interpolation mesh
             *
             * @param aPdvTypes PDV_Type types; set->group->individual
             * @param aMeshIndex Interpolation mesh index
             */
            void create_ip_pdv_hosts(Cell<Cell<Cell<PDV_Type>>> aPdvTypes, moris_index aMeshIndex = 0);

            /**
             * Create PDV_Type hosts with PDVs for each of the spatial dimensions on the integration mesh
             *
             * @param aMeshIndex Integration mesh index
             */
            void create_ig_pdv_hosts(moris_index aMeshIndex = 0);

            /**
             * Assign PDV hosts based on properties constructed through parameter lists
             */
            void assign_pdv_hosts();

        private:

            /**
             * Get the index of analytic geometry
             *
             * @param aGlobalGeometryIndex Global geometry index (active geometry)
             * @return Index for the analytic geometry cell
             */
            size_t analytic_geometry_index(size_t aGlobalGeometryIndex);

            /**
             * Get the index of discrete geometry
             *
             * @param aGlobalGeometryIndex Global geometry index (active geometry)
             * @return Index for the discrete geometry cell
             */
            size_t discrete_geometry_index(size_t aGlobalGeometryIndex);

            /**
             * Compute_intersection_info, calculates the relevant intersection information placed in the geometry object
             *
             * @param[in]  aEntityNodeInds - node to entity connectivity
             * @param[in]  aNodeVars       - node level set values
             * @param[in]  aCheckType      - if a entity local location is necessary 1, else 0.
             * @param[out] Returns an intersection flag and local coordinates if aCheckType 1 in cell 1 and node sensitivity information in cell 2 if intersection point located
             **/
            bool compute_intersection_info( moris::moris_index               const & aEntityIndex,
                                            moris::Matrix< moris::IndexMat > const & aEntityNodeInds,
                                            moris::Matrix< moris::DDRMat >   const & aNodeCoords,
                                            moris::size_t                    const & aCheckType,
                                            moris::Matrix< moris::IndexMat >       & aNodeADVIndices,
                                            GEN_Geometry_Object                    & aGeometryObject );

            /**
             * Interpolate given level set values to a child node location
             *
             * @param aParentTopology
             * @param aGeometryIndex
             * @param aNodeLocalCoordinate
             * @param aLevelSetValues
             */
            void interpolate_level_set_value_to_child_node_location( xtk::Topology                  const & aParentTopology,
                                                                     moris::size_t                  const & aGeometryIndex,
                                                                     moris::Matrix< moris::DDRMat > const & aNodeLocalCoordinate,
                                                                     moris::Matrix< moris::DDRMat >       & aLevelSetValues );

        };
    }
}

#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_Geometry_Engine_HPP_ */
