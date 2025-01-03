/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Property.hpp
 *
 */

#pragma once

#include "cl_GEN_Design_Field.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_GEN_Design.hpp"

namespace moris::gen
{
    /**
     * This struct contains additional parameters that are used by properties.
     */
    struct Property_Parameters : public Field_Parameters
            , public Design_Parameters
    {
        PDV_Type              mPDVType;              //! The type of PDV that this property will be assigned to
        bool                  mInterpolationPDV;     //! If the PDV is defined on the interpolation mesh (always true for now)
        Vector< uint >        mPDVMeshSetIndices;    //! Mesh set indices for assigning PDVs
        Vector< std::string > mPDVMeshSetNames;      //! Mesh set names for assigning PDVs

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList
         */
        explicit Property_Parameters( const Parameter_List& aParameterList = prm::create_gen_property_parameter_list( gen::Field_Type::NONE ) );
    };

    class Property : public Design_Field
            , public Design
    {
      private:
        Property_Parameters mParameters;

      public:
        /**
         * Constructor taking in a field pointer and a set of parameters.
         *
         * @param aField Field for computing nodal values
         * @param aParameters Field parameters
         * @param aNodeManager Node manager from the geometry engine, if available
         */
        explicit Property(
                std::shared_ptr< Field >   aField,
                const Property_Parameters& aParameters  = Property_Parameters(),
                Node_Manager&              aNodeManager = Node_Manager::get_trivial_instance() );

        /**
         * Sets a new node manager (from the geometry engine, if it was created after this property)
         *
         * @param aNodeManager Geometry engine node manager
         */
        void set_node_manager( Node_Manager& aNodeManager );

        /**
         * Updates the dependencies of this design based on the given designs
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedDesigns All designs (this design will take fields from the ones it needs)
         */
        void update_dependencies( const Vector< std::shared_ptr< Design > >& aAllUpdatedDesigns ) override;

        /**
         * Gets the PDV type that this property defines.
         *
         * @return PDV type
         */
        PDV_Type get_pdv_type();

        /**
         * Gets if this property's PDV type is defined on the interpolation mesh.
         *
         * @return mesh type switch
         */
        bool is_interpolation_pdv();

        /**
         * Gets the mesh set indices where this property's PDV is defined.
         *
         * @param aMesh Mesh for getting set indices from set names
         * @return Mesh set indices
         */
        Vector< uint > get_pdv_mesh_set_indices( mtk::Integration_Mesh* aMesh );

        /**
         * Used for writing to mtk meshes and printing for debug info
         *
         * @param aNodeIndex decides the point at which the field value is printed. If the node is a derived node, the value is interpolated from the parents.
         * @param aCoordinates The field location to get the value from.
         * @return the value of the property field at the requested location
         */
        void get_design_info(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Vector< real >&         aOutputDesignInfo ) override;

        /**
         * gets the number of fields the property has
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
        std::shared_ptr< Field > get_field() override
        {
            return Design_Field::mField;
        }

        /**
         * Gets the name of the geometry
         *
         */
        std::string get_name() override;

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
    };
}    // namespace moris::gen
