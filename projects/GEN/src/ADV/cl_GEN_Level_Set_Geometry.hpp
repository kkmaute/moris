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

namespace moris::ge
{
    enum class Int_Interpolation
    {
        LINEAR,
        MULTILINEAR
    };

    /**
     * This is a struct used to simplify \ref moris::ge::Level_Set_Geometry constructors. It contains all field and level-set parameters.
     */
    struct Level_Set_Parameters : public Field_Parameters
    {
        real mIsocontourThreshold;      // Level set isocontour level
        real mIsocontourTolerance;      // Interface tolerance based on geometry value
        real mIntersectionTolerance;    // Interface tolerance based on intersecction distance

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Parameter list with level set geometry parameters
         */
        explicit Level_Set_Parameters( const ParameterList& aParameterList = prm::create_level_set_geometry_parameter_list() );
    };

    class Level_Set_Geometry : public Geometry
            , public Design_Field
            , public std::enable_shared_from_this< Level_Set_Geometry >    // TODO make it so we don't need enable_shared_from_this, should be possible in intersection node
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
                std::shared_ptr< Field > aField,
                Level_Set_Parameters     aParameters  = Level_Set_Parameters(),
                Node_Manager&            aNodeManager = Node_Manager::get_trivial_instance() );

        /**
         * Sets a new node manager (from the geometry engine, if it was created after this geometry)
         *
         * @param aNodeManager Geometry engine node manager
         */
        void set_node_manager( Node_Manager& aNodeManager ) override;

        /**
         * Gets the mode of intersection used for this geometry
         *
         * @return Intersection_Mode enum
         */
        Intersection_Mode get_intersection_mode();

        /**
         * Accesses the isocontour level that determines the interface for this geometry
         *
         * @return the isocontour level that determines the geometry interface
         */
        real get_isocontour_threshold();

        /**
         * Acccesses the isocontour tolerance for this geometry
         *
         * @return isocontour tolerance
         */
        real get_isocontour_tolerance();

        /**
         * Accesses the intersection tolerance for this geometry
         *
         * @return The real value of the intersection tolerance
         */
        real get_intersection_tolerance();

        /**
         * Gets if this geometry depends on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() override;

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
         * @param aBaseNodes Base nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return New intersection node
         */
        Intersection_Node* create_intersection_node(
                uint                     aNodeIndex,
                const Cell< Node* >&     aBaseNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder ) override;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() override;

        /**
         * Determines the geometric region of a point based on a level set value
         *
         * @param aLevelSetValue Value of the level set function
         * @return Geometric region enum
         */
        Geometric_Region determine_geometric_region( real aLevelSetValue );

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
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        void reset_nodal_data() override;

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
        virtual void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs,
                const Matrix< DDSMat >&       aSharedADVIds,
                uint                          aADVOffsetID ) override;

        /**
         * Used to print geometry information to exodus files and print debug information.
         *
         *  @param aNodeIndex decides the point at which the field value is printed. If the node is a derived node, the value is interpolated from the parents.
         * @param aCoordinates The field location to get the value from.
         * @return the value of the level set field at the requested location
         */
        void get_design_info(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Cell< real >&           aOutputDesignInfo ) override;

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
        std::shared_ptr< Field > get_field()
        {
            return Design_Field::mField;
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
    };
}    // namespace moris::ge
