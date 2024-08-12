/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Surface_Mesh_Geometry.hpp
 *
 */

#pragma once

#include "cl_SDF_Object.hpp"

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_ADV_Manager.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_Library_IO.hpp"

namespace moris::gen
{
    // User-defined function that determines which indices are fixed or not
    using VERTEX_FACTOR_FUNCTION = real ( * )( const uint aFacetVertexIndex, const Matrix< DDRMat >& aFacetVertexCoordinates, const uint aDimension );

    /**
     * This is a struct used to simplify \ref moris::gen::Surface_Mesh_Geometry constructors. It contains all field and surface mesh parameters.
     */
    struct Surface_Mesh_Parameters : public Field_Parameters
            , public Design_Parameters
    {
        Vector< real > mOffsets;                     // Initial shift of surface mesh coordinates
        Vector< real > mScale;                       // Option to scale each axis of the surface mesh
        std::string    mFilePath;                    // Surface mesh file path
        real           mIntersectionTolerance;       // Interface tolerance based on intersection distance
        Vector< uint > mADVIndices;                  // Indices of the ADVs that the surface mesh depends on
        std::string    mVertexFactorFunctionName;    // Name of the user-defined function that provides a scaling factor for the facet vertex sensitivities
        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Surface_Mesh_Parameters( const Parameter_List& aParameterList = prm::create_surface_mesh_geometry_parameter_list() );
    };

    class Surface_Mesh_Geometry : public Geometry
            , public sdf::Object // FIXME: remove sdf::Object inheritance
    {
      private:
        bool mBasesComputed = false;

        Surface_Mesh_Parameters mParameters;
        Node_Manager*           mNodeManager;
        std::string             mName;

        // Optimization variables
        VERTEX_FACTOR_FUNCTION             mVertexFactorFunction = nullptr;
        mtk::Mesh*                         mMesh                 = nullptr;    // Pointer to lagrange interpolation mesh
        Vector< uint >                     mFixedVertexIndices;                // Indices of surface mesh vertices that are unaffected by ADVs
        Vector< std::shared_ptr< Field > > mPerturbationFields;                // Vector of perturbation fields
        Matrix< DDRMat >                   mVertexBases;                       // Basis function values for each vertex <number of fields> x <number of vertices>
        Vector< mtk::Cell* >               mVertexBackgroundElements;          // Index of the background element the facet vertex was in on construction
        Vector< Vector< real > >           mOriginalVertexCoordinates;         // All vertex coordinates as they were upon construction <dimension> x <number of vertices>


      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        Surface_Mesh_Geometry(
                mtk::Mesh*                    aMesh,
                ADV_Manager&                  aADVManager,
                Surface_Mesh_Parameters       aParameters  = Surface_Mesh_Parameters(),
                Node_Manager&                 aNodeManager = Node_Manager::get_trivial_instance(),
                std::shared_ptr< Library_IO > aLibrary     = nullptr );

        /**
         * Default destructor
         *
         */
        ~Surface_Mesh_Geometry();

        /**
         * Gets the geometric region of a node, based on this geometry.
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Geometric region enum
         */
        Geometric_Region get_geometric_region(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates ) override;

        /**
         * Creates an intersection node based on the given information. The intersection node may or may not represent an intersection;
         * that is, its position may lie outside of the edge definition based on the given nodal coordinates. This information can be
         * requested from the created intersection node.
         *
         * @param aNodeIndex Node index of the new intersection node
         * @param aBackgroundNodes Background nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return New intersection node
         */
        Intersection_Node* create_intersection_node(
                uint                              aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                mtk::Geometry_Type                aBackgroundGeometryType,
                mtk::Interpolation_Order          aBackgroundInterpolationOrder ) override;

        /**
         * Computes the local coordinate along a parent edge of an intersection node created using this geometry.
         *
         * @param aBackgroundNodes Background nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aParentFacetIndex return value. A pointer to the facet that intersected the edge to create this intersection node
         * @return Parent edge local coordinate, between -1 and 1
         */
        real compute_intersection_local_coordinate(
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                uint&                             aParentFacet );

        /**
         *
         * Whether or not the surface mesh has ADVs
         *
         */
        bool
        depends_on_advs() const;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        Vector< std::shared_ptr< mtk::Field > > get_mtk_fields() override;

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs ) override;

        /**
         * Gets the name of this geometry
         *
         * @return File name of the .obj file that this surface mesh was created with
         */
        std::string
        get_name() override;

        /**
         * Gets the names of all the fields associated with this design
         *
         * @return Vector< std::string > The name of each perturbation field
         */
        virtual Vector< std::string > get_field_names() override;

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         *
         * @param aInterpolationMesh Interpolation mesh containing new nodal data
         */
        void reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh ) override;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         * @param aFieldIndex For geometries that have multiple fields, which field to discretize
         */
        void discretize(
                mtk::Mesh_Pair    aMeshPair,
                sol::Dist_Vector* aOwnedADVs ) override;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMTKField Input MTK field to map based on
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         * @param aFieldIndex For geometries that have multiple fields, which field to discretize
         */
        void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs ) override;

        /**
         * Used to print geometry information to exodus files and print debug information.
         *
         *  @param aNodeIndex decides the point at which the surface mesh displacement is printed. If the node is a derived node, the value is interpolated from the parents.
         * @param aCoordinates The field location to get the value from.
         * @return the value of the surface mesh displacement at the requested location
         */
        void get_design_info(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Vector< real >&         aOutputDesignInfo ) override;

        /**
         * Appends this designs ADV IDs, ijklIDs, lower bounds, and upper bounds to the global matrices stored in the geometry engine.
         * Sets mNumCoeff, mOffsetID, and appends to mSharedADVIDs
         *
         * @param aMeshPair     Background mesh pair
         * @param aOwnedADVIds  IDs of the ADVs that are owned by this design
         * @param aOwnedijklIDs
         * @param aOffsetID     Offset of this Design's ADVs from the global ADV vector
         * @param aLowerBounds  ADV lower bounds
         * @param aUpperBounds  ADV upper bounds
         * @return uint The new offset ID after this geometry has appended its information
         */
        virtual sint append_adv_info(
                mtk::Interpolation_Mesh* aMesh,
                Vector< sint >&          aOwnedADVIds,
                Matrix< IdMat >&         aOwnedijklIDs,
                sint                     aOffsetID,
                Vector< real >&          aLowerBounds,
                Vector< real >&          aUpperBounds ) override;

        /**
         * Gets the number of fields the surface mesh has
         */
        uint get_num_fields() override
        {
            return mPerturbationFields.size();
        }

        /**
         * Allows for access to the perturbation fields
         *
         * @return Underlying fields
         */
        Vector< std::shared_ptr< Field > > get_fields() override;

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @param aADVs ADVs
         */
        void set_advs( sol::Dist_Vector* aADVs ) override;

        /**
         * Sets a new node manager (from the geometry engine, if it was created after this geometry)
         *
         * @param aNodeManager Geometry engine node manager
         */
        void set_node_manager( Node_Manager& aNodeManager ) override;

        /**
         * Gets if this field is to be used for seeding a B-spline field.
         *
         * @return Logic for B-spline creation
         */
        bool intended_discretization() override;

        /**
         * Gets a discretization mesh index for a discretized field.
         *
         * @return Mesh index
         */
        moris_index get_discretization_mesh_index() override;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        real get_discretization_lower_bound() override;

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        real get_discretization_upper_bound() override;

        /**
         * Updates the dependencies of this design based on the given designs
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedDesigns All designs (this design will take fields from the ones it needs)
         */
        void update_dependencies( Vector< std::shared_ptr< Design > > aAllUpdatedDesigns ) override;

        /**
         * Gets the intersection tolerance for creating intersection nodes
         *
         */
        real get_intersection_tolerance() override
        {
            return this->Geometry::get_intersection_tolerance();
        }

        /**
         * Gets the center of the facet at the given local index
         *
         * @param aFacetIndex local index of the facet in the surface mesh
         * @return Matrix< DDRMat > center of the facet
         */
       Matrix< DDRMat > get_facet_center( uint aFacetIndex );
        
        /**
         * Gets the basis functions for the specified vertex
         *
         * @param aVertexIndex Index of the vertex
         * @return Matrix< DDRMat > Basis functions for the vertex
         */
        Matrix< DDRMat > get_vertex_bases( uint aVertexIndex )
        {
            return mVertexBases.get_column( aVertexIndex );
        }

        /**
         * Determines if the requested facet vertex depends on ADVs or not
         *
         * @param aFacetVertexIndex the index of the facet vertex that is queried
         * @return true if the vertex's index is not in mParameters.mFixedVertexIndices and the node's position is within the boundaries of the mesh
         * @return false if either of the above conditions are true
         */
        bool facet_vertex_depends_on_advs( uint aFacetVertexIndex );

        /**
         * Computes and returns the sensitivity of a facet vertex with respect to the ADVs
         * NOTE: This function assumes that the facet vertex depends on ADVs. Check this with facet_vertex_depends_on_advs() if unsure
         *
         * @return Matrix< DDRMat > derivative of global vertex location with respect to each ADV. Size is <dimension> x <number of ADVs>
         */
        Matrix< DDRMat > get_dvertex_dadv( uint aFacetVertexIndex );

        /**
         * Gets the ADV IDs that the facet vertex depends on.
         * These are the ADVs that control the bspline field value in the background element that the vertex lies in.
         * NOTE: This function assumes that the facet vertex depends on ADVs. Check this with facet_vertex_depends_on_advs() if unsure
         *
         *
         * @param aFacetVertexIndex Vertex index of the surface mesh
         * @return Matrix< DDSMat > ADV IDs that the vertex depends on
         */
        Vector< sint > get_vertex_adv_ids( uint aFacetVertexIndex );

        /**
         * Gets the IDs of the ADVs that the given node depends on
         *
         * @param aNodeIndex the query node index on the integration mesh for which ADVs to retrieve
         * @param aCoordinates the query node coordinates for which ADVs to retrieve
         * @return Vector< sint > ADV IDs that the query node depends on
         */
        Vector< sint > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

      private:
        /**
         * Finds the background elemenent in aField that contains aCoordinates
         *
         * @param aCoordinate Search global coordinate location
         * @param aBoundingBox Return variable that holds the bounding box of the found cell
         *
         * @return Index of the element in which aCoordinates resides. If no element is found, -1 is returned
         */
        mtk::Cell* find_background_element_from_global_coordinates( const Matrix< DDRMat >& aCoordinate );

        /**
         * Gets the bounding box of a requested mtk::Cell
         *
         * @param aElement mtk::Cell of which to get the bounding box
         * @return Vector< Vector< real > > 2 x dim 2D vector. First index is the minimum second is the maximum for each dimension
         */
        Vector< Vector< real > > determine_mtk_cell_bounding_box( mtk::Cell* aElement );

        /**
         * @brief Computes the basis functions at a given point in the background element.
         *
         * @param aBackgroundElement the background element in which the point resides
         * @param aParametricCoordinates the local coordinate of the point
         * @param aBasis Return value. The basis functions at the point.
         */
        Matrix< DDRMat > compute_vertex_basis(
                mtk::Cell*              aBackgroundElement,
                const Matrix< DDRMat >& aParametricCoordinates );

        /**
         * Updates the values of the basis functions for all of the facet vertices and the background elements they lie in
         *
         */
        void update_vertex_basis_data();

        /**
         * Determines the field value at a given point in the background element.
         *
         * @param aBackgroundElement The background element in which the point resides
         * @param aFieldIndex        The perturbation field index (spatial dimension)
         * @param aBasis  The facet vertex to perturb
         * @return real
         */
        real interpolate_perturbation_from_background_element(
                mtk::Cell* aBackgroundElement,
                uint       aFieldIndex,
                uint       aFacetVertexIndex );
    };
}    // namespace moris::gen