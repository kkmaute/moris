/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Stored_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris::ge
{
    class Stored_Field : public Field_Discrete_Integration
    {

    private:
        std::shared_ptr< Field > mField;
        mtk::Mesh* mMesh;
        Cell< real > mFieldValues;

    public:

        /**
         * Constructor
         *
         * @param aMesh The mesh pointer where node information can be obtained
         * @param aField Field for obtaining values to store
         */
      Stored_Field(
                mtk::Mesh*               aMesh,
                std::shared_ptr< Field > aField );

        /**
         * Given a node index or coordinate, returns the field value.
         *
         * @param aNodeIndex Node index
         * @return Field value
         */
        real get_field_value(uint aNodeIndex);

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to its ADVs.
         *
         * @param aNodeIndex Node index
         * @return Vector of sensitivities
         */
        const Matrix<DDRMat>& get_dfield_dadvs(uint aNodeIndex);

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @return Determining ADV IDs at this node
         */
        Matrix<DDSMat> get_determining_adv_ids(uint aNodeIndex);

        /**
         * Resets all nodal information about field values.
         */
        void reset_nodal_data();

    private:

        /**
         * Evaluates and stores the nodal values of this field for use later.
         */
        void evaluate_nodal_values();

    };
}
