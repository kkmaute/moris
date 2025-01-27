/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.hpp
 *
 */

#pragma once

#include "cl_GEN_Design_Field.hpp"
#include "cl_GEN_Geometry.hpp"
#include "GEN_Data_Types.hpp"

namespace moris::gen
{
    /**
     * This is a struct used to simplify \ref moris::gen::Level_Set_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Level_Set_Parameters : public Field_Parameters
            , public Design_Parameters
    {
        real mIsocontourThreshold;      // Level set isocontour level
        real mIsocontourTolerance;      // Interface tolerance based on geometry value
        real mIntersectionTolerance;    // Interface tolerance based on intersection distance

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Level_Set_Parameters( const Parameter_List& aParameterList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::NONE ) );
    };

    class Level_Set_Geometry : public Geometry
            , public Design_Field
    {
      private:
        Level_Set_Parameters mParameters;

      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         * @param aNodeManager Node manager from the geometry engine, if available
         */
        explicit Level_Set_Geometry(
                std::shared_ptr< Field >    aField,
                const Level_Set_Parameters& aParameters  = Level_Set_Parameters(),
                Node_Manager&               aNodeManager = Node_Manager::get_trivial_instance() );

        /**
         * Sets a new node manager (from the geometry engine, if it was created after this geometry)
         *
         * @param aNodeManager Geometry engine node manager
         */
        void set_node_manager( Node_Manager& aNodeManager ) override;

        /**
         * Accesses the isocontour level that determines the interface for this geometry
         *
         * @return the isocontour level that determines the geometry interface
         */
        real get_isocontour_threshold() const;

        /**
         * Gets if this geometry depends on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() const override;

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
         * Creates a floating node based on the given information.
         *
         * @param aNodeIndex Node index to be assigned to the new floating node
         * @param aBackgroundNodes Background nodes of the element where the floating node lies
         * @param aParametricCoordinates Parametric coordinates inside the background element
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return New floating node
         */
        Floating_Node* create_floating_node(
                uint                              aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Matrix< DDRMat >&           aParametricCoordinates,
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
        real compute_intersection_local_coordinate(
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aParentNode Parent node
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                const Basis_Node& aParentNode,
                Matrix< DDRMat >& aSensitivities ) const;

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
         * Gets the name of this design's field
         *
         * @return Name
         */
        std::string get_name() override;

        /**
         * Gets the names of all the fields associated with this design
         *
         * @return Vector< std::string > the geometry name, as this implementation only has one field
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
         */
        void discretize(
                mtk::Mesh_Pair    aMeshPair,
                sol::Dist_Vector* aOwnedADVs ) override;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMTKField Input MTK field to map based on
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         */
        void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs ) override;

        /**
         * Used to print geometry information to exodus files and print debug information.
         *
         *  @param aNodeIndex decides the point at which the field value is printed. If the node is a derived node, the value is interpolated from the parents.
         * @param aCoordinates The field location to get the value from.
         * @return the value of the level set field at the requested location
         */
        void get_design_info(
                const uint              aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Vector< real >&         aOutputDesignInfo ) override;

        /**
         * Gets the number of fields the level set geometry has
         */
        uint get_num_fields() override
        {
            return 1;
        }

        /**
         * Allows for access to the GEN field
         *
         * @return Underlying field
         */
        Vector< std::shared_ptr< Field > > get_fields() override
        {
            return { Design_Field::mField };
        }

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @param aADVs ADVs
         */
        void set_advs( sol::Dist_Vector* aADVs ) override
        {
            Design_Field::mField->set_advs( aADVs );
        }

        /**
         * Gets if this field is to be used for seeding a B-spline field.
         *
         * @return Logic for B-spline creation
         */
        bool intended_discretization() const override;

        /**
         * Gets a discretization mesh index for a discretized field.
         *
         * @return Mesh index
         */
        moris_index get_discretization_mesh_index() const override;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        real get_discretization_lower_bound() const override;

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        real get_discretization_upper_bound() const override;

        /**
         * Updates the dependencies of this design based on the given designs
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedDesigns All designs (this design will take fields from the ones it needs)
         */
        void update_dependencies( const Vector< std::shared_ptr< Design > >& aAllUpdatedDesigns ) override;

      private:
        /**
         * Determines the geometric region of a point based on a level set value
         *
         * @param aLevelSetValue Value of the level set function
         * @return Geometric region enum
         */
        Geometric_Region determine_geometric_region( real aLevelSetValue ) const;
    };
}    // namespace moris::gen
