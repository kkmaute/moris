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

#include "fn_SDF_Raycast.hpp"

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::gen
{
    /**
     * This is a struct used to simplify \ref moris::gen::Surface_Mesh_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Surface_Mesh_Parameters : public Field_Parameters
            , public Design_Parameters
    {
        Vector< real > mOffsets;                  // Initial shift of surface mesh coordinates
        Vector< real > mScale;                    // Option to scale each axis of the surface mesh
        std::string    mFilePath;                 // Surface mesh file path
        real           mIntersectionTolerance;    // Interface tolerance based on intersection distance
        Vector< uint > mADVIndices;

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Surface_Mesh_Parameters( const ParameterList& aParameterList = prm::create_surface_mesh_geometry_parameter_list() );
    };

    class Surface_Mesh_Geometry : public Geometry
            , public sdf::Object
    {
      private:
        Surface_Mesh_Parameters mParameters;
        std::string             mName;
        mtk::Mesh*              mMesh;
        Cell< Cell< real > >    mOriginalVertexCoordinates;    // All vertex coordinates as they were upon construction <dimension> x <number of vertices>

        Cell< std::shared_ptr< Field > > mPerturbationFields;    // Vector of perturbation fields

      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         */
        Surface_Mesh_Geometry( mtk::Mesh* aMesh, Matrix< DDRMat > aADVs, Surface_Mesh_Parameters aParameters = Surface_Mesh_Parameters() );

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
         * @return Parent edge local coordinate, between -1 and 1
         */
        virtual real compute_intersection_local_coordinate(
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode ) override;

        /**
         *
         * Whether or not the surface mesh has ADVs
         *
         */
        bool
        depends_on_advs() const
        {
            // TODO BRENDAN
            return false;
        }

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

        std::string
        get_name() override
        {
            return mName;
        }

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
         */
        void discretize(
                mtk::Mesh_Pair          aMeshPair,
                sol::Dist_Vector*       aOwnedADVs,
                const Matrix< DDSMat >& aSharedADVIds,
                uint                    aADVOffsetID ) override;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMTKField Input MTK field to map based on
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         */
        void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs,
                const Matrix< DDSMat >&       aSharedADVIds,
                uint                          aADVOffsetID ) override;

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
         * Gets the number of fields the surface mesh has
         */
        uint get_num_fields() override
        {
            return mPerturbationFields.size();
        }

        /**
         * Allows for access to the GEN field
         *
         * @return Underlying field'
         */
        std::shared_ptr< Field > get_field()
        {
            // TODO BRENDAN
            return nullptr;
        }

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @param aADVs ADVs
         */
        void set_advs( sol::Dist_Vector* aADVs ) override;

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
        virtual moris_index get_discretization_mesh_index() override;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        virtual real get_discretization_lower_bound() override;

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        virtual real get_discretization_upper_bound() override;

        /**
         * Gets the intersection tolerance for creating intersection nodes
         *
         */
        real get_intersection_tolerance() override
        {
            return this->Geometry::get_intersection_tolerance();
        }

      private:
        void transform_surface_mesh_to_local_coordinate(
                const Parent_Node& aFirstParentNode,
                const Parent_Node& aSecondParentNode );

        /**
         * Finds the background elemenent in aField that contains aCoordinates
         *
         * @param aCoordinate Search global coordinate location
         * @param aBoundingBox Return variable that holds the bounding box of the found cell
         *
         * @return Index of the element in which aCoordinates resides. If no element is found, -1 is returned
         */
        moris_index find_background_element_from_global_coordinates(
                const Matrix< DDRMat >& aCoordinate,
                Cell< Cell< real > >&   aBoundingBox );

        real interpolate_perturbation_from_background_element(
                mtk::Cell*              aBackgroundElement,
                uint                    aFieldIndex,
                const Matrix< DDRMat >& aParametricCoordinates );
    };
}    // namespace moris::gen
