/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field.hpp
 *
 */

#pragma once

#include <utility>
#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "cl_MTK_Field.hpp"

/**
 * \def ADV_ARG_TYPES
 * All the ADV arguments to a field
 */
#define ADV_ARG_TYPES Matrix< DDRMat >& aADVs, Matrix< DDUMat > aFieldVariableIndices, Matrix< DDUMat > aADVIndices, Matrix< DDRMat > aConstants, std::string aName = ""

/**
 * \def ADV_ARGS
 * Just the ADV arguments a field receives, to pass to the base constructor
 */
#define ADV_ARGS aADVs, aFieldVariableIndices, aADVIndices, aConstants, aName

/**
 * \def VARIABLE_CHECK( num_variables )
 * Generates a moris error that checks that the number of variables passed to a field is correct. This should be called in most child field constructors.
 */
#define VARIABLE_CHECK( num_variables ) MORIS_ERROR( aFieldVariableIndices.length() + aConstants.length() == num_variables, \
        "A GEN %s must be created with a total of exactly %d variables (ADVs + constants)", __FUNCTION__, num_variables )

namespace moris::mtk
{
    class Mesh;
}

namespace moris::gen
{
    // Forward declare derived node
    class Derived_Node;

    class Field
    {
      protected:
        uint mNumOriginalNodes = 0; // FIXME
        ADV_Manager mADVManager;
        Matrix< DDRMat > mSensitivities;

      private:
        std::string mName;
        inline static uint mCount = 0;

      public:
        // FIXME this is only used for Field_Analytic, this should be removed when refinement is done through a proper interface
        mtk::Mesh_Pair mMeshPairForAnalytic = mtk::Mesh_Pair( nullptr, nullptr );
        static inline uint gDiscretizationIndex = 0;

        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         */
        Field( Matrix< DDRMat >& aADVs,
                Matrix< DDUMat > aFieldVariableIndices,
                Matrix< DDUMat > aADVIndices,
                Matrix< DDRMat > aConstants,
                std::string      aName );

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        Field( Matrix< DDRMat > aConstants,
                std::string     aName );

        /**
         * Constructor that sets all field variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aFieldVariableIndices Field variable indices for assigning the shared ADV IDs
         * @param aSharedADVIds Shared ADV IDs needed for this field
         */
        Field( const Matrix< DDSMat >& aSharedADVIds,
                std::string            aName );

        /**
         * Copy constructor with replacement variables for new constants.
         *
         * @param aCopy Analytic field to copy
         * @param aReplaceVariables Variable indices to replace
         * @param aNewConstants New constants
         */
        Field( const Field& aCopy,
                const Vector< uint >& aReplaceVariables,
                const Vector< real >& aNewConstants );

        /**
         * Destructor
         */
        virtual ~Field() = default;

        /**
         * Copies the current field into a shared pointer with replacement variables for new constants.
         * Right now, the only practical use for this is to translate analytic field, so the default
         * implementation for discrete fields does not copy. This can be changed in the future.
         *
         * @param aReplaceVariables Variable indices to replace
         * @param aNewConstants New constants
         * @return Shared pointer to copied field
         */
        virtual std::shared_ptr< Field > copy(
                const Vector< uint >& aReplaceVariables,
                const Vector< real >& aNewConstants );

        /**
         * Sets the ADVs and grabs the field variables needed from the ADV vector
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADVs
         */
        template< typename Vector_Type >
        void set_advs( Vector_Type& aADVs )
        {
            mADVManager.set_advs( aADVs );
        }

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        virtual void import_advs( sol::Dist_Vector* aOwnedADVs );

        /**
         * Gets the number of reference coordinates this field has.
         *
         * @return Number of reference coordinates
         */
        virtual uint get_number_of_reference_coordinates();

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        virtual real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Gets a field value of a derived node.
         *
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         * @return Field value
         */
        virtual real get_field_value(
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) = 0;

        /**
         * Gets an interpolated field value based on given basis nodes.
         *
         * @param aBasisNodes Basis nodes of a derived node
         * @return Field value
         */
        real get_interpolated_field_value(
                const Vector< Basis_Node >& aBasisNodes,
                const Node_Manager&       aNodeManager );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return d(field value)/d(ADV_j)
         */
        virtual const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) = 0;

        /**
         * Gets a vector of the field derivatives with respect to ADVs of a derived node.
         *
         * @param aSensitivities Sensitivities to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        virtual void get_dfield_dadvs(
                Matrix< DDRMat >&   aSensitivities,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) = 0;

        /**
         * Interpolates the field derivatives with respect to ADVs based on given basis nodes.
         *
         * @param aInterpolatedSensitivities Appended d(field value)/d(ADV_j)
         * @param aBasisNodes Basis nodes of a derived node
         */
        void append_interpolated_dfield_dadvs(
                Matrix< DDRMat >&         aInterpolatedSensitivities,
                const Vector< Basis_Node >& aBasisNodes,
                const Node_Manager&       aNodeManager,
                real                      aBasisFactor = 1.0 );

        /**
         * Gets the IDs of ADVs that this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        virtual Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets the IDs of ADVs that this field depends on for evaluations at a derived node.
         *
         * @param aDeterminingADVIDs Determining ADV IDs to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        virtual void get_determining_adv_ids(
                Matrix< DDSMat >&   aDeterminingADVIDs,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager );

        /**
         * Gets the interpolated IDs of ADVs that this field depends on given basis nodes.
         *
         * @param aInterpolatedADVIDs Appended determining ADV IDs
         * @param aBasisNodes Basis nodes of a derived node
         */
        void append_interpolated_determining_adv_ids(
                Matrix< DDSMat >&         aInterpolatedADVIDs,
                const Vector< Basis_Node >& aBasisNodes,
                const Node_Manager&       aNodeManager );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        virtual void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) = 0;

        /**
         * Gets if this manager has ADVs (at least one non-constant parameter)
         *
         * @return if this manager has ADVs
         */
        bool has_advs();

        /**
         * Sets the dependencies of this field after they have been found by the owning field. By default
         * does nothing.
         *
         * @param aDependencyFields Other fields that this field depends on.
         */
        virtual void set_dependencies( Vector< std::shared_ptr< Field > > aDependencyFields );

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        virtual void reset_nodal_data( mtk::Interpolation_Mesh* aMesh );

        /**
         * Gets the name of this field
         *
         * @return Name
         */
        std::string get_name();

        /**
         * Gets an MTK field, if this field needs to be remapped to a new mesh
         *
         * @return MTK field
         */
        virtual std::shared_ptr< mtk::Field > get_mtk_field() = 0;

      protected:

        /**
         * Creates an MTK field based on this field on a given mesh. Should only need to be called by discrete fields.
         *
         * @param aMesh Mesh pointer
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > create_mtk_field( const mtk::Mesh_Pair& aMeshPair );

      private:

        /**
         * Checks that a unique name has been given, and if not assigns a default name
         */
        void verify_name();

    };
}
