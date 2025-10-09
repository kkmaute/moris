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
        Vector< uint >        mNumberOfRefinements;        // The number of refinement steps to use for this field
        Vector< uint >        mRefinementMeshIndices;      // Indices of meshes to perform refinement on
        sint                  mRefinementFunctionIndex;    // Index of a user-defined refinement function (default = -1)
        Vector< GQI_Type >    mRequestedGQIs;              // List of GQIs to be computed on this design
        Vector< std::string > mRequestedGQINames;          // Names for the requested GQIs

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

      protected:
        uint                     mOffsetID;        // Offset of the global ADVs to this Design's ADVs
        Vector< Vector< sint > > mSharedADVIDs;    // IDs of the ADVs that this design shares. Size = number of fields

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
                const uint              aNodeIndex,
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
         * Gets the names of all the fields associated with this design
         *
         * @return Vector< std::string >
         */
        virtual Vector< std::string > get_field_names() = 0;

        /**
         * Gets if this field is to be used for seeding a B-spline field.
         *
         * @return Logic for B-spline creation
         */
        virtual bool intended_discretization() const = 0;

        /**
         * Gets a discretization mesh index for a discretized field.
         *
         * @return Mesh index
         */
        virtual moris_index get_discretization_mesh_index() const = 0;

        /**
         * Gets the lower bound for a discretized field.
         *
         * @return Lower bound
         */
        virtual real get_discretization_lower_bound() const = 0;

        /**
         * Get the upper bound for a discretized field.
         *
         * @return Upper bound
         */
        virtual real get_discretization_upper_bound() const = 0;

        /**
         * Allows for access to all the GEN Fields for this design
         *
         * @return Underlying fields
         */
        virtual Vector< std::shared_ptr< Field > > get_fields() = 0;

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
                Vector< real >&          aUpperBounds,
                uint                     aFieldIndex = 0 );

        static void communicate_missing_owned_coefficients(
                mtk::Interpolation_Mesh* aMesh,
                Matrix< IdMat >&         aAllCoefIds,
                Matrix< IdMat >&         aAllCoefOwners,
                Matrix< IdMat >&         aAllCoefijklIds,
                uint                     aNumCoeff,
                uint                     aDiscretizationMeshIndex,
                mtk::MeshType            aMeshType );

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @param aADVs ADVs
         */
        virtual void set_advs( sol::Dist_Vector* aADVs ) = 0;

        /**
         * Check if the design depends on ADVs
         */
        virtual bool depends_on_advs() const = 0;

        /**
         * Updates the dependencies of this design based on the given designs
         * (fields may have been mapped/updated).
         *
         * @param aAllUpdatedDesigns All designs (this design will take fields from the ones it needs)
         */
        virtual void update_dependencies( const Vector< std::shared_ptr< Design > >& aUpdatedFields ) = 0;


        //------------------------------------------------------------------------------
        // Geometry Quantity of Interest (GQI) functions brendan document all
        //------------------------------------------------------------------------------

      public:
        const uint get_num_GQIs() const;

        const Vector< std::string >& get_all_GQI_names() const;

        /**
         * Loops through all GQIs requested on this design, computes their values, and if requested,
         * computes their sensitivities and stores them in the given distributed vector.
         *
         * @param aGQISensitivities Distributed vector to store GQI sensitivities in. Contains all GQI sensitivities for all designs.
         * @param aRequestIndices Vector indices in aGQISensitivities that this design's GQI sensitivities should be stored in. MORIS_UINT_MAX if the GQI is not requested.
         *
         * @return Vector< real > Values of ALL the GQIs for this design. NOTE: This is the size of mParameters.mRequestedGQIs
         */
        Vector< real > compute_GQIs( sol::Dist_Vector* aGQISensitivities, const Vector< uint >& aRequestIndices );

      protected:
        /**
         * Computes the value of a requested geometric quantity of interest (GQI) for this design.
         */
        virtual real compute_GQI( GQI_Type aGQIType ) = 0;

        virtual void compute_GQI_sensitivities( GQI_Type aGQIType, sol::Dist_Vector* aGQISensitivities, uint aRequestIndex ) const = 0;
    };
}    // namespace moris::gen
