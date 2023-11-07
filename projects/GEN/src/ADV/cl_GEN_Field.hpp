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
#include "cl_GEN_Child_Node.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

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

namespace moris::ge
{
    class Field
    {
      public:
        mtk::Mesh_Pair mMeshPair = mtk::Mesh_Pair( nullptr, nullptr );

      protected:
        uint mNumOriginalNodes = 0;
        ADV_Manager mADVManager;
        Matrix< DDRMat > mSensitivities;

      private:
        std::string mName;

      public:

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
                const Cell< uint >& aReplaceVariables,
                const Cell< real >& aNewConstants );

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
                const Cell< uint >& aReplaceVariables,
                const Cell< real >& aNewConstants );

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
         * Gets the IDs of ADVs that this manager depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        virtual Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * Gets if this manager has ADVs (at least one non-constant parameter)
         *
         * @return if this manager has ADVs
         */
        bool has_advs();

        /**
         * Sets the dependencies of this field after they have been found by the owning property. By default
         * does nothing.
         *
         * @param aDependencyFields Other fields that this field depends on.
         */
        virtual void set_dependencies( Cell< std::shared_ptr< Field > > aDependencyFields );

        /**
         * Add a new child node for evaluation, implemented for discrete integration fields.
         *
         * @param aNodeIndex Index of the child node
         * @param aChildNode Contains information about how the child node was created
         */
        virtual void add_child_node(
                uint                          aNodeIndex,
                std::shared_ptr< Child_Node > aChildNode );

        /**
         * In relevant derived classes, uses additional information from the given interpolation mesh to define
         * potentially new nodes on the field. Implemented for discrete interpolation fields.
         *
         * @param aMesh Interpolation mesh with additional nodes
         */
        virtual void add_nodal_data( mtk::Interpolation_Mesh* aMesh );

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         */
        virtual void reset_nodal_data();

        void set_num_original_nodes( uint aNumOriginalNodes );

        /**
         * Gets the name of this field
         *
         * @return Name
         */
        std::string get_name();

    };
}
