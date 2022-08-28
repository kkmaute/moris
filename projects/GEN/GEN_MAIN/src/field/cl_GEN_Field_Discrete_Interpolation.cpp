/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Discrete_Interpolation.cpp
 *
 */

#include "cl_GEN_Field_Discrete_Interpolation.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Discrete_Interpolation::Field_Discrete_Interpolation(mtk::Interpolation_Mesh* aMesh)
        {
            if (aMesh)
            {
                this->add_nodal_data(aMesh);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Discrete_Interpolation::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            // FIXME: map onto base node not correct however coordinate is consistent with incorrect node
            return this->get_base_field_value( (this->*get_node_index)(aNodeIndex), aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Field_Discrete_Interpolation::get_dfield_dadvs(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_base_dfield_dadvs( (this->*get_node_index)(aNodeIndex), aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_determining_adv_ids(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_base_determining_adv_ids( (this->*get_node_index)(aNodeIndex), aCoordinates );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete_Interpolation::get_dfield_dcoordinates(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            MORIS_ERROR(false, "A discrete interpolation field does not have spatial derivatives.");
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete_Interpolation::add_nodal_data(mtk::Interpolation_Mesh* aMesh)
        {
            // Obtain all data needed
            tNodeToBaseIndex.resize(aMesh->get_num_nodes(), 1);
            for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
            {
                tNodeToBaseIndex(tNodeIndex) = aMesh->get_base_node_index(tNodeIndex);
            }

            // Make sure that we get node indices using this data
            get_node_index = &Field_Discrete_Interpolation::get_base_node_index;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete_Interpolation::reset_nodal_data()
        {
            get_node_index = &Field_Discrete_Interpolation::return_same_index;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_base_determining_adv_ids(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return Field::get_determining_adv_ids(
                    (this->*get_node_index)(aNodeIndex),
                    aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Field_Discrete_Interpolation::return_same_index(uint aNodeIndex)
        {
            return aNodeIndex;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Field_Discrete_Interpolation::get_base_node_index(uint aNodeIndex)
        {
            // until Property PDVs are defined on nodes of un-zipped mesh return trivial relation
            return aNodeIndex; //tNodeToBaseIndex(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

