#ifndef MORIS_CL_Geometry_Engine_HPP_
#define MORIS_CL_Geometry_Engine_HPP_

// WRK
#include "cl_WRK_Performer.hpp"

// GEN
#include "cl_GEN_Pending_Node.hpp"
#include "cl_GEN_Phase_Table.hpp"

#include "cl_GEN_Geometry_Object.hpp"
#include "cl_GEN_Geometry_Object_Manager.hpp"

#include "pdv/cl_GEN_Pdv_Host_Manager.hpp"
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

    namespace MSI
    {
        class Design_Variable_Interface;
    }

    //------------------------------------------------------------------------------------------------------------------

    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------
        // GEOMETRY ENGINE
        //--------------------------------------------------------------------------------------------------------------

        class Geometry_Engine : public wrk::Performer
        {
        private:

            // Level-set isocontour threshold
            real mIsocontourThreshold;
            real mPerturbationValue;

            // Spatial dimensions
            uint mSpatialDim;

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

            // Child nodes
            Cell<std::shared_ptr<Child_Node>> mChildNodes;

            // Contains all the geometry objects
            Geometry_Object_Manager mGeometryObjectManager;

            // Contains all the pdv hosts
            Pdv_Host_Manager mPdvHostManager;

            // Phase Table
            Phase_Table mPhaseTable;

            // Mesh
            std::shared_ptr<mtk::Mesh_Manager> mMeshManager;

            Matrix<IndexMat> mInterfaceNodeIndices;

        public:

            /**
             * Constructor using cell of cell of parameter lists
             *
             * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
             * @param aLibrary Library used for pulling user-defined functions
             */
            Geometry_Engine(Cell<Cell<ParameterList>> aParameterLists,
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
                            real                            aIsocontourThreshold = 0.0,
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
             * Gets all of the geometry field values at the specified coordinates
             *
             * @param aNodeIndices Node indices on the mesh
             * @param aCoordinates Coordinate values for evaluating the geometry fields
             * @param aGeometryIndex Index of the geometry for evaluating the field of
             * @return Field values
             */
            real get_geometry_field_value(      uint            aNodeIndex,
                                          const Matrix<DDRMat>& aCoordinates,
                                                uint            aGeometryIndex = 0);

            /**
             * create new node geometry objects
             * @param[ in ] aNodeCoords node coordinates
             */
            void create_new_node_geometry_objects( const Cell<moris_index>&    aNewNodeIndices,
                                                   bool                        aStoreParentTopo,
                                                   const Cell<xtk::Topology*>& aParentTopo,
                                                   const Cell<Matrix<DDRMat>>& aParamCoordRelativeToParent,
                                                   const Matrix<DDRMat>&       aGlobalNodeCoord );

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

            /**
             * Sets the indices for the nodes on the interface
             *
             * @param aInterfaceNodeIndices Interface node indices
             */
            void set_interface_nodes( Matrix< IndexMat > const & aInterfaceNodeIndices);

            /**
             * Computes the intersection of an isocountour with an entity and returning the local coordinate relative to the parent
             * and the global coordinate if needed
             */
            void get_intersection_location(
                    const Matrix<DDRMat>&   aGlobalNodeCoordinates,
                    const Matrix<DDRMat>&   aEntityNodeVars,
                    const Matrix<IndexMat>& aEntityNodeIndices,
                    Matrix<DDRMat>&         aIntersectionLocalCoordinates,
                    Matrix<DDRMat>&         aIntersectionGlobalCoordinates,
                    bool                    aCheckLocalCoordinate = false,
                    bool                    aComputeGlobalCoordinate = false);

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
              * For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
              */
            size_t get_phase_index(moris_index aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * @brief Provided the inside and out phase values for an entity, return the phase index
             */
            moris_index get_elem_phase_index(moris::Matrix< moris::IndexMat > const & aElemOnOff);

            /**
             * @brief Returns whether a node is inside or outside wrt to a given geometry index
             */
            size_t get_node_phase_index_wrt_a_geometry(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, uint aGeometryIndex);

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
             * Register an MTK mesh pair to the geometry engine
             *
             * @param aMeshManager MTK mesh manager with interpolation and integration meshes
             */
            void register_mesh(std::shared_ptr<mtk::Mesh_Manager> aMeshManager);

            /**
             * Return the number of fields that can be used for refinement
             *
             * @return Number of fields for refinement
             */
            uint get_num_refinement_fields();

            /**
             * Gets a flag to determine if refinement should continue
             *
             * @param aFieldIndex The index of a field
             * @param aRefinementIndex The current refinement step being performed
             * @return If refinement is needed for this field
             */
            bool refinement_needed(uint aFieldIndex,
                                   uint aRefinementIndex);

            /**
             * Returns fields so that HMR can perform refinement based on the data from this performer
             *
             * @param aFieldIndex Index of the field
             * @param aNodeIndex Index of the node
             * @param aCoordinates Coordinates of the node
             */
            real get_field_value(uint aFieldIndex,
                                 uint aNodeIndex,
                                 const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the index of an HMR user-defined refinement function for the given field index
             *
             * @param aFieldIndex Index of the field
             * @param aRefinementIndex The current refinement step being performed
             * @return User-defined function index, or -1 to use default refinement
             */
            sint get_refinement_function_index(uint aFieldIndex,
                                               uint aRefinementIndex);

            /**
             * @brief assign the pdv type and property for each pdv host in a given set
             */
            void assign_ip_hosts_by_set_name( std::string                 aSetName,
                                              std::shared_ptr<Property> aPropertyPointer,
                                              PDV_Type                    aPdvType,
                                              moris_index                 aWhichMesh = 0 );

            /**
             * @brief assign the pdv type and property for each pdv host in a given set
             */
            void assign_ip_hosts_by_set_index( moris_index               aSetIndex,
                                               std::shared_ptr<Property> aPropertyPointer,
                                               PDV_Type                  aPdvType,
                                               moris_index               aWhichMesh = 0 );

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

            /**
             * Resets the stored info specific to meshes
             */
            void reset();

        private:

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
                                            GEN_Geometry_Object                    & aGeometryObject );

        };
    }
}

#endif /* MORIS_CL_Geometry_Engine_HPP_ */
