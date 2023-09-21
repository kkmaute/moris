/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Discrete_Interpolation.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris::ge
{
    class Field_Discrete_Interpolation : public Field
    {
    private:
        uint (Field_Discrete_Interpolation::*get_node_index)(uint) = &Field_Discrete_Interpolation::return_same_index;
        Matrix<DDUMat> tNodeToBaseIndex;

    public:

        /**
         * Constructor using pointers to ADVs for variable evaluations.
         *
         * @tparam Vector_Type Type of vector where ADVs are stored
         * @param aADVs ADV vector
         * @param aFieldVariableIndices Indices of field variables to be filled by the ADVs
         * @param aADVIndices The indices of the ADV vector to fill in the field variables
         * @param aConstants The constant field variables not filled by ADVs
         */
        template< typename Vector_Type >
        Field_Discrete_Interpolation(
              Vector_Type&     aADVs,
              Matrix< DDUMat > aFieldVariableIndices,
              Matrix< DDUMat > aADVIndices,
              Matrix< DDRMat > aConstants,
              mtk::Interpolation_Mesh* aMesh )
              : Field( aADVs, aFieldVariableIndices, aADVIndices, aConstants )
        {
            if (aMesh)
            {
                this->add_nodal_data( aMesh );
            }
        }

        /**
         * Constructor using only constants (no ADVs).
         *
         * @param aConstants The parameters that define this field
         */
        Field_Discrete_Interpolation(
                Matrix< DDRMat > aConstants,
                mtk::Interpolation_Mesh* aMesh );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Field value
         */
        real get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinate, returns a matrix all sensitivities.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix<DDSMat> get_determining_adv_ids(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities);

        /**
         * Uses additional information from the given interpolation mesh to define
         * potentially new nodes on the field.
         *
         * @param aMesh Interpolation mesh with additional nodes
         */
        void add_nodal_data(mtk::Interpolation_Mesh* aMesh);

        /**
         * Resets all nodal information. This should be called when a new XTK mesh is being created.
         */
        void reset_nodal_data();

    private:

        /**
         * Given a node index, returns the field value at the base node index
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        virtual real get_base_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates) = 0;

        /**
         * Given a node index, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Vector of sensitivities
         */
        virtual const Matrix<DDRMat>& get_base_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates) = 0;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        virtual Matrix<DDSMat> get_base_determining_adv_ids(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates);

        /**
         * Returns the same node index, for the case if there is no mesh information to go off of.
         *
         * @return Input node index
         */
        uint return_same_index(uint aNodeIndex);

        /**
         * Gets a base node index based on the information obtained from the mesh with nodal data.
         *
         * @return Base node index
         */
        uint get_base_node_index(uint aNodeIndex);

    };
}
