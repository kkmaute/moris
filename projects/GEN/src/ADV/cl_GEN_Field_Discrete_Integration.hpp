/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Discrete_Integration.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

namespace moris::ge
{
    class Field_Discrete_Integration : public Field
    {

      private:
        Cell< std::shared_ptr< Child_Node > > mChildNodes;

      public:

        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         * @param aName Name of this field
         * @param aNumOriginalNodes Number of nodes originally on the integration mesh, before intersections
         */
        Field_Discrete_Integration(
                ADV_ARG_TYPES,
                uint aNumOriginalNodes = 0 )
                : Field( ADV_ARGS )
        {
            mNumOriginalNodes = aNumOriginalNodes;
        }

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        Field_Discrete_Integration(
                Matrix< DDRMat > aConstants,
                uint             aNumOriginalNodes );

        /**
         * Constructor that sets all field variables as consecutive ADVs. Assumes the use of distributed ADVs.
         *
         * @param aSharedADVIds Shared ADV IDs needed for this field
         * @param aNumOriginalNodes Number of original nodes for this field
         * @param aName Name of this field
         */
        Field_Discrete_Integration(
                const Matrix< DDSMat >& aSharedADVIds,
                uint                    aNumOriginalNodes,
                std::string             aName );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node index, returns the field value
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        virtual real get_field_value( uint aNodeIndex ) = 0;

        /**
         * Given a node index or coordinate, returns a matrix all sensitivities.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Given a node index, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @return Vector of sensitivities
         */
        virtual const Matrix< DDRMat >& get_dfield_dadvs( uint aNodeIndex ) = 0;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations, including child nodes.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations for non-child nodes.
         *
         * @param aNodeIndex Node index
         * @return Determining ADV IDs at this node
         */
        virtual Matrix< DDSMat > get_determining_adv_ids( uint aNodeIndex );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;

        /**
         * Add a new child node for evaluation.
         *
         * @param aNodeIndex Index of the child node
         * @param aChildNode Contains information about how the child node was created
         */
        void add_child_node( uint aNodeIndex, std::shared_ptr< Child_Node > aChildNode ) override;

        /**
         * Resets all child nodes, called when a new XTK mesh is being created.
         */
        void reset_nodal_data() override;
    };
}