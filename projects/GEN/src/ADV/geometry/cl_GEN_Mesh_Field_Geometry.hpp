/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Mesh_Field_Geometry.hpp
 *
 */

#pragma once

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace moris::ge
{
    class Mesh_Field_Geometry: public Field_Discrete_Integration
    {
        private:

            mtk::Mesh*     mMesh;
            std::string    mFieldName;
            real           mOffset;
            mtk::EntityRank     mEntityRank;

            Matrix<DDRMat> mFieldData;
            bool           mUseOwnData;

        public:
            /**
             * Constructor
             *
             * @param aMesh       Mesh with the level set fields
             * @param aFieldName  Name of the field
             */
            Mesh_Field_Geometry(mtk::Mesh*  aMesh,
                    std::string aFieldName,
                    mtk::EntityRank  aEntityRank = mtk::EntityRank::NODE,
                    Level_Set_Parameters    aParameters = Level_Set_Parameters() );

            /**
              * Constructor
              *
              * @param aMesh In-core mesh
              * @param aFileName    Name of the file with field data
              * @param aFieldName   Name of the field
              * @param aFileFormat  Name of file format (e.g., exodus)
              *
              */
            Mesh_Field_Geometry(
                            mtk::Mesh*  aMesh,
                            std::string aFileName,
                            std::string aFieldName,
                            std::string aFileFormat,
                            real        aOffset,
                            mtk::EntityRank  aEntityRank = mtk::EntityRank::NODE,
                    Level_Set_Parameters aParameters = Level_Set_Parameters());

            /**
             * Given a node index, returns the field value.
             *
             * @param aNodeIndex Node index
             * @return Distance to this geometry
             */
            real get_field_value(uint aNodeIndex);

        private:

            /**
             * Given a node index, evaluates the sensitivity of the geometry field with respect to all of the
             * geometry variables. This is currently not implemented for a mesh field geometry.
             *
             * @param aNodeIndex Node index
             * @return Vector of sensitivities
             */
            const Matrix<DDRMat>& get_dfield_dadvs(uint aNodeIndex);

    };
}