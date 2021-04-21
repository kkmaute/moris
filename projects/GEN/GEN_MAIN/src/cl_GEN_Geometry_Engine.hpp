#ifndef MORIS_CL_Geometry_Engine_HPP_
#define MORIS_CL_Geometry_Engine_HPP_

// MRS
#include "cl_Param_List.hpp"
#include "cl_Library_IO.hpp"
#include "fn_trans.hpp"

// WRK
#include "cl_WRK_Performer.hpp"

// GEN
#include "st_GEN_Geometry_Engine_Parameters.hpp"
#include "cl_GEN_Phase_Table.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Geometric_Proximity.hpp"

// MTK
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_Mesh_Enums.hpp"
#include <unordered_map>

// SOL FIXME
#include "cl_SOL_Dist_Vector.hpp"

namespace xtk
{
    class Topology;
}

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

        class Geometry_Engine : public wrk::Performer
        {
        // Fields
        protected:
            Cell<std::shared_ptr<Geometry>> mGeometries;
            Cell<std::shared_ptr<Property>> mProperties;

        private:
            // Phase Table
            Phase_Table mPhaseTable;

            // Intersections
            Intersection_Mode mIntersectionMode;
            real mIsocontourThreshold;
            real mIsocontourTolerance;
            real mIntersectionTolerance;

            // Spatial dimensions
            uint mNumSpatialDimensions;

            // ADVs
            Matrix<DDRMat> mInitialPrimitiveADVs;
            Matrix<DDRMat> mADVs;
            Matrix<DDSMat> mFullADVIds;
            sol::Dist_Vector* mOwnedADVs = nullptr;
            sol::Dist_Vector* mPrimitiveADVs = nullptr;

            // Bounds
            Matrix<DDRMat> mLowerBounds;
            Matrix<DDRMat> mUpperBounds;

            // IQIs
            Cell<std::string> mRequestedIQIs;

            // Library
            std::shared_ptr<Library_IO> mLibrary;

            size_t mActiveGeometryIndex = 0;
            std::string mGeometryFieldFile;
            std::string mOutputMeshFile;
            bool mShapeSensitivities = false;
            real mTimeOffset;

            // PDVs
            Pdv_Host_Manager mPdvHostManager;
            std::shared_ptr<Intersection_Node> mQueuedIntersectionNode;

            // Keeps track of a vertex proximity to each geometry ( NumVerts x NumGeometries)
            // 0 - G(x) < threshold
            // 1 - G(x) == threshold
            // 2 - G(x) > threshold
            // Max not set
            Cell<Geometric_Proximity> mVertexGeometricProximity;

        public:

            /**
             * Constructor using cell of cell of parameter lists
             *
             * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
             * @param aLibrary Library used for pulling user-defined functions
              * @param aMesh Mesh for discrete or mesh based geomtries
             */
            Geometry_Engine(
                    Cell< Cell<ParameterList> > aParameterLists,
                    std::shared_ptr<Library_IO> aLibrary = nullptr,
                    mtk::Mesh*                  aMesh = nullptr);

            /**
             * Constructor
             *
             * @param aMesh Mesh for getting B-spline information
             * @param aParameters Optional geometry engine parameters
             */
            Geometry_Engine(
                    mtk::Interpolation_Mesh*   aMesh,
                    Geometry_Engine_Parameters aParameters = {});

            /**
             * Destructor
             */
            ~Geometry_Engine();

            /**
             * Sets new ADVs for the geometry engine.
             *
             * @param aNewADVs vector of new advs to use
             */
            void set_advs(const Matrix<DDRMat> & aNewADVs);

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
             * Gets the sensitivities of the criteria with respect to the advs
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
             * Determines if the element consisting of the given node coordinates is intersected.
             *
             * @param aNodeIndices Node indices
             * @param aNodeCoordinates Node coordinates
             * @return If the element is intersected
             */
            bool is_intersected(
                    const Matrix<IndexMat>& aNodeIndices,
                    const Matrix<DDRMat>  & aNodeCoordinates);

            /**
             * Determines if the given edge is intersected, and queues an intersection node if it is. If an intersection
             * node has been queued, questions can be asked about the queued node:
             *
             * @param aNodeIndex1 First node index
             * @param aNodeIndex2 Second node index
             * @param aFirstNodeLocalCoordinates Local coordinates of the first node inside the background element
             * @param aSecondNodeLocalCoordinates Local coordinates of the second node inside the background element
             * @param aFirstNodeGlobalCoordinates Global coordinates of the first node
             * @param aSecondNodeLocalCoordinates Global coordinates of the second node
             * @param aBackgroundElementNodeIndices Node indices of the background element
             * @param aBackgroundElementNodeCoordinates Node coordinates of the background element
             * @return If the edge is intersected and a node has been queued
             */
            bool queue_intersection(
                    uint                        aFirstNodeIndex,
                    uint                        aSecondNodeIndex,
                    const Matrix<DDRMat>&       aFirstNodeLocalCoordinates,
                    const Matrix<DDRMat>&       aSecondNodeLocalCoordinates,
                    const Matrix<DDRMat>&       aFirstNodeGlobalCoordinates,
                    const Matrix<DDRMat>&       aSecondNodeGlobalCoordinates,
                    const Matrix<DDUMat>&       aBackgroundElementNodeIndices,
                    const Cell<Matrix<DDRMat>>& aBackgroundElementNodeCoordinates);

            /**
             * Returns if the queued intersection has the first parent node on the active geometry interface.
             *
             * @return If the first parent node is on the interface
             */
            bool queued_intersection_first_parent_on_interface();

            /**
             * Returns if the queued intersection has the second parent node on the active geometry interface.
             *
             * @return If the second parent node is on the interface
             */
            bool queued_intersection_second_parent_on_interface();

            /**
             * Gets the local coordinate of the queued intersection node.
             *
             * @return Intersection node local coordinate (between -1 and 1)
             */
            real get_queued_intersection_local_coordinate();

            /**
             * Gets the global coordinates of the queued intersection node.
             *
             * @return Intersection node global coordinates
             */
            Matrix<DDRMat> get_queued_intersection_global_coordinates();

            /**
             * Admit the queued intersection as a unique, permanent node(s) for sensitivity calculations.
             */
            void admit_queued_intersection(uint aNodeIndex);

            /**
             * Update the queued intersection node with its node ID and node owner.
             *
             * @param aNodeIndex Node Index
             * @param aNodeId Node ID
             * @param aNodeOwner Node owner
             */
            void update_queued_intersection(
                    const moris_index & aNodeIndex,
                    const moris_index & aNodeId,
                    const moris_index & aNodeOwner );

            /**
             * create new node geometry objects
             * @param[ in ] aNodeCoords node coordinates
             */
            void create_new_child_nodes(
                    const Cell<moris_index>&    aNewNodeIndices,
                    const Cell<xtk::Topology*>& aParentTopo,
                    const Cell<Matrix<DDRMat>>& aParamCoordRelativeToParent,
                    const Matrix<DDRMat>&       aGlobalNodeCoord );

            /**
             * Get the total number of phases in the phase table
             */
            size_t get_num_phases();

            /**
             * For a given node index, return the phase index relative to each geometry (i.e. inside/outside indicator)
             */
            size_t get_phase_index(
                    moris_index           aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Use the geometric proximity to tell me whether this vertex is on the interface wrt a provided geometry indx
             */
            moris_index 
            is_interface_vertex(moris_index aNodeIndex,
                                moris_index aGeometryIndex);

            /**
             * @brief Provided the inside and out phase values for an entity, return the phase index
             */
            moris_index get_elem_phase_index(Matrix< IndexMat > const & aElemOnOff);

            /**
             * @brief Returns whether a node is inside or outside wrt to a given geometry index
             */
            size_t get_node_phase_index_wrt_a_geometry(
                    uint                  aNodeIndex,
                    uint                  aGeometryIndex);

            /**
             * @brief Returns the number of geometries
             */
            size_t get_num_geometries();

            /**
             * @brief Returns the number of phases
             */
            size_t get_num_bulk_phase();

            /**
             * @brief Returns the active geometry index
             */
            size_t get_active_geometry_index();

            /**
             * @brief Advance the active geometry index
             */
            void advance_geometry_index();

            /**
             * Return the number of fields that can be used for refinement
             *
             * @return Number of fields for refinement
             */
            uint get_num_refinement_fields();

            const Matrix< DDSMat > & get_num_refinements(uint aFieldIndex );

            const Matrix< DDSMat > & get_refinement_mesh_indices(uint aFieldIndex );

            /**
             * Gets all of the MTK fields that the geometry engine is using.
             *
             * @return MTK fields
             */
            Cell<std::shared_ptr<mtk::Field>> get_mtk_fields();

            /**
             * Returns fields so that HMR can perform refinement based on the data from this performer
             *
             * @param aFieldIndex Index of the field
             * @param aNodeIndex Index of the node
             * @param aCoordinates Coordinates of the node
             */
            real get_field_value(
                    uint                  aFieldIndex,
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates);

            /**
             * Gets the index of an HMR user-defined refinement function for the given field index
             *
             * @param aFieldIndex Index of the field
             * @param aRefinementIndex The current refinement step being performed
             * @return User-defined function index, or -1 to use default refinement
             */
            sint get_refinement_function_index(
                    uint aFieldIndex,
                    uint aRefinementIndex);

            /**
             * Discretize GEN fields on the given mesh and distribute parallel ADVs based on these fields.
             *
             * @param aMeshPair Mesh for discretizing fields
             */
            void distribute_advs(
                    mtk::Mesh_Pair aMeshPair,
                    moris::Cell<std::shared_ptr< mtk::Field > > aFields);

            /**
             * Resets the information that the geometry engine stores about a mesh.
             *
             * @param aMesh Mesh for computing level set data
             */
            void reset_mesh_information(mtk::Interpolation_Mesh* aMesh);

            /**
             * Outputs geometry and property fields on the given mesh, and writes level set fields to a text file.
             * Uses output locations and file names stored from a parameter list or previous call to an output function.
             *
             * @param aMesh Mesh to evaluate fields on
             */
            void output_fields(mtk::Mesh* aMesh);

            /**
             * Creates geometry and property fields on the given mesh, and writes the mesh to an exodus file.
             *
             * @param aMesh Mesh to evaluate fields on
             * @param aExodusFileName Name of an Exodus file to write the mesh to
             */
            void output_fields_on_mesh(
                    mtk::Mesh*  aMesh,
                    std::string aExodusFileName);

            /**
             * Writes all geometry fields to separate text files with the given base file name (suffix appended on to
             * identify each individual geometry by index.
             *
             * @param aMesh Mesh to evaluate fields on
             * @param aFileName Base name of text files to write the geometry field data to
             */
            void write_geometry_fields(
                    mtk::Mesh*  aMesh,
                    std::string aBaseFileName);

            /**
             * Assign PDV hosts based on properties constructed through parameter lists
             *
             * @param aMeshManager Mesh manager
             */
            void create_pdvs(mtk::Mesh_Pair aMeshPair);

        private:

            /**
             * Create PDV_Type hosts with the specified PDV_Type types on the interpolation mesh
             *
             * @param aPdvTypes PDV_Type types; set->group->individual
             * @param aMeshIndex Interpolation mesh index
             */
            void create_interpolation_pdv_hosts(
                    mtk::Interpolation_Mesh    * aInterpolationMesh,
                    mtk::Integration_Mesh      * aIntegrationMesh,
                    Cell<Cell<Cell<PDV_Type>>>   aPdvTypes);

            /**
             * Create PDV_Type hosts with PDVs for each of the spatial dimensions on the integration mesh
             *
             * @param aMeshIndex Integration mesh index
             */
            void set_integration_pdv_types(mtk::Integration_Mesh* aIntegrationMesh);

            /**
             * @brief assign the pdv type and property for each pdv host in a given set
             */
            void assign_property_to_pdv_hosts(
                    std::shared_ptr<Property> aPropertyPointer,
                    PDV_Type                  aPdvType,
                    mtk::Integration_Mesh*    aIntegrationMesh,
                    Matrix<DDUMat>            aSetIndices);

            /**
             * Initialize the PDV type list.
             */
            void initialize_pdv_type_list();

            /**
             * Setup initial geometric proximities
             */ 
            void
            setup_initial_geometric_proximities(mtk::Interpolation_Mesh* aMesh);

            /*
            *   Return the geometric proximity index. Converts from a double value to a index
            */
            moris_index
            get_geometric_proximity_index(moris::real const & aGeometricVal);

            /*
            * Determine the intersection vertex proximity
            */
           moris_index
           get_queued_intersection_geometric_proximity_index(moris_index const & aGeomIndex);

            /**
             * Admit the proximity information for an interface vertex
             */
            void
            admit_queued_intersection_geometric_proximity(uint aNodeIndex);

            /**
             * Decides how to construct the phase table based on the given arguments.
             *
             * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
             * @return Phase table
             */
            static Phase_Table create_phase_table(
                    Cell<Cell<ParameterList>> aParameterLists,
                    std::shared_ptr<Library_IO> aLibrary);

            /**
             * Decides how to construct the phase table based on the given parameter lists
             *
             * @param aNumGeometries Number of geometries
             * @param aBulkPhases Phase table vector
             * @param aPhaseFunction Custom phase function
             * @return Phase table
             */
            static Phase_Table create_phase_table(
                    uint                     aNumGeometries,
                    Matrix<DDUMat>           aBulkPhases,
                    PHASE_FUNCTION aPhaseFunction = nullptr,
                    uint                     aNumPhases = 1);
        };
    }
}

#endif /* MORIS_CL_Geometry_Engine_HPP_ */
