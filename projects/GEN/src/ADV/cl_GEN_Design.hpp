/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Design.hpp
 *
 */

#pragma once

#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_GEN_ADV_Handler.hpp"
#include "cl_GEN_Field.hpp"
namespace moris::gen
{
    /**
     * This is a struct used to simplify \ref moris::gen::Design_Field constructors. It contains additional parameters that
     * are used by all fields, with given defaults.
     */
    struct Design_Parameters
    {
        Vector< uint > mNumberOfRefinements;        // The number of refinement steps to use for this field
        Vector< uint > mRefinementMeshIndices;      // Indices of meshes to perform refinement on
        sint         mRefinementFunctionIndex;    // Index of a user-defined refinement function (-1 = default)

        /**
         * Constructor with a given parameter list
         *
         * @param aParameterList Design field parameter list
         */
        explicit Design_Parameters( const Parameter_List& aParameterList );

        /**
         * Default constructor
         */
        Design_Parameters();
    };

    class Design
    {
      private:
        Design_Parameters mParameters;

      public:
        /**
         * Constructor
         *
         * @param aParameters Parameters relevant to all designs
         */
        explicit Design( Design_Parameters aParameters );

        /**
         * Default destructor
         */
        virtual ~Design() = default;

        /**
         * This function will return true when called less than the number of refinements set for this field,
         * and false otherwise. This is to determine for a given refinement call if this field needs refinement.
         *
         * @return if to perform an additional refinement with this field
         */
        const Vector< uint >& get_num_refinements();

        const Vector< uint >& get_refinement_mesh_indices();

        /**
         * Gets the index of a user-defined refinement function used within HMR.
         *
         * @return User-defined refinement function index
         */
        sint get_refinement_function_index();

        /**
         * Returns the value of the design at the given node location. Often times, this is the field value of the design.
         *
         * @param aNodeIndex
         */
        virtual void get_design_info(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Vector< real >&         aOutputDesignInfo ) = 0;

        /**
         * Gets the number of fields that the design has
         *
         * @return Number of fields
         */
        virtual uint get_num_fields() = 0;

        /**
         * Gets the name of the design
         *
         * @return Design name
         */
        virtual std::string get_name() = 0;

        /**
         * Gets if this field is to be used for seeding a B-spline field.
         *
         * @return Logic for B-spline creation
         */
        virtual bool intended_discretization() = 0;

        /**
         * Gets a discretization mesh index for a discretized field.
         *
         * @return Mesh index
         */
        virtual moris_index get_discretization_mesh_index() = 0;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        virtual real get_discretization_lower_bound() = 0;

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        virtual real get_discretization_upper_bound() = 0;

        /**
         * Allows for access to the GEN field
         *
         * @return Underlying field
         */
        virtual std::shared_ptr< Field > get_field() = 0;

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @param aADVs ADVs
         */
        virtual void set_advs( sol::Dist_Vector* aADVs ) = 0;

        /**
         * Updates the dependencies of this design based on the given designs
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedDesigns All designs (this design will take fields from the ones it needs)
         */
        virtual void update_dependencies( Vector< std::shared_ptr< Design > > aAllUpdatedDesigns ) = 0;
    };
}    // namespace moris::gen
