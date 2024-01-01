/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Data_STK.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_DATA_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_DATA_STK_HPP_

// Third-party header files.
#include <stk_io/StkMeshIoBroker.hpp>             // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>          // for count_entities
#include <stk_mesh/base/MetaData.hpp>             // for MetaData
#include <stk_mesh/base/BulkData.hpp>             // for BulkData
#include <stk_mesh/base/Selector.hpp>             // for Selector
#include <stk_mesh/base/FEMHelpers.hpp>           // for Selector
#include "stk_io/DatabasePurpose.hpp"             // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/CoordinateSystems.hpp"    // for Cartesian
#include "stk_mesh/base/CreateFaces.hpp"          // for handling faces
#include "stk_mesh/base/CreateEdges.hpp"          // for handling faces
#include "stk_mesh/base/Bucket.hpp"               // for buckets
#include "stk_mesh/base/Field.hpp"                // for coordinates
#include "stk_mesh/base/GetEntities.hpp"          // for coordinates
#include "stk_mesh/base/FieldParallel.hpp"        // for handling parallel fields

#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Block_Set.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Matrix_Field_Info.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"

#include "cl_Vector.hpp"

// For Vector and vertex APi
#include "cl_MTK_Cell_STK.hpp"
#include "cl_MTK_Vertex_STK.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh_Data_STK
        {
          public:
            Mesh_Data_STK()
            {
            }

            ~Mesh_Data_STK()
            {
                delete mMeshReader;
                delete mMtkMeshBulkData;
                delete mMtkMeshMetaData;
            }

            // STK specific Member variables
            stk::io::StkMeshIoBroker* mMeshReader      = nullptr;
            stk::mesh::MetaData*      mMtkMeshMetaData = nullptr;
            stk::mesh::BulkData*      mMtkMeshBulkData = nullptr;

            // General mesh trait member variables
            bool mDataGeneratedMesh = false;
            bool mCreatedFaces      = false;
            bool mCreatedEdges      = false;

            // Local to Global c
            // Vector(0) - Node Local to Global
            // Vector(1) - Edge Local to Global
            // Vector(2) - Face Local to Global
            // Vector(3) - Elem Local to Global
            Vector< moris::Matrix< IdMat > > mEntityLocaltoGlobalMap;

            // Local to Global c
            // Vector(0) - Node Global to Local
            // Vector(1) - Edge Global to Local
            // Vector(2) - Face Global to Local
            // Vector(3) - Elem Global to Local
            // Keep in mind not all of these are created
            Vector< std::unordered_map< moris_id, moris_index > > mEntityGlobalToLocalMap;

            // Exterior Vector Entity Rank (Same structure as local to global
            //  Interior Vector (processor rank)
            Vector< Vector< moris::Matrix< IndexMat > > > mEntitySendList;
            Vector< Vector< moris::Matrix< IndexMat > > > mEntityReceiveList;

            //    // Vector and Vertex
            Vector< mtk::Cell_STK >                 mMtkCells;
            Vector< mtk::Vertex_Core_STK >          mMtkVertices;
            Vector< mtk::Vertex_Interpolation_STK > mMtkVerticeInterpolation;

            // cell connectivity
            Vector< std::shared_ptr< mtk::Cell_Info > > mCellInfo;

            uint mMaxNumFields = 20;
            uint mNumDims      = 0;

            typedef stk::mesh::Field< real >                          Field1CompReal;
            typedef stk::mesh::Field< real, stk::mesh::Cartesian2d >  Field2CompReal;
            typedef stk::mesh::Field< real, stk::mesh::Cartesian3d >  Field3CompReal;
            typedef stk::mesh::Field< real, stk::mesh::FullTensor22 > Field4CompReal;
            typedef stk::mesh::Field< real, stk::mesh::FullTensor >   Field9CompReal;

            std::vector< Field1CompReal* > mField1CompVecsReal;
            std::vector< Field2CompReal* > mField2CompVecsReal;
            std::vector< Field3CompReal* > mField3CompVecsReal;
            std::vector< Field4CompReal* > mField4CompVecsReal;
            std::vector< Field9CompReal* > mField9CompVecsReal;

            typedef stk::mesh::Field< sint >                          Field1CompInt;
            typedef stk::mesh::Field< sint, stk::mesh::Cartesian2d >  Field2CompInt;
            typedef stk::mesh::Field< sint, stk::mesh::Cartesian3d >  Field3CompInt;
            typedef stk::mesh::Field< sint, stk::mesh::FullTensor22 > Field4CompInt;
            typedef stk::mesh::Field< sint, stk::mesh::FullTensor >   Field9CompInt;

            std::vector< Field1CompInt* > mField1CompVecsInt;
            std::vector< Field2CompInt* > mField2CompVecsInt;
            std::vector< Field3CompInt* > mField3CompVecsInt;
            std::vector< Field4CompInt* > mField4CompVecsInt;
            std::vector< Field9CompInt* > mField9CompVecsInt;

            std::vector< bool > mSetRankFlags;    // Flags for user-defined node [0], side [1], and block [2] sets.

            // TODO: ADD edge sets (not sure why these were neglected in previous implementation
            std::vector< std::vector< std::string > > mSetNames;    // User-defined names for node [0], side [1], and block [2] sets.

            std::map< uint, uint > mProcsSharedToIndex;

            // Fields to Declare on Output (note this is needed for supplementary fields
            // provided when mesh is loaded from a file only
            Vector< Field1CompReal* > mRealNodeScalarFieldsToAddToOutput;

            // Shared Vertex Information
            std::unordered_map< moris_index, moris_index > mVertexSharingProcsMap;    // maps between processor rank and location in cell of mVertexSharingData
            moris::Matrix< moris::IdMat >                  mProcsWithSharedVertex;
            Vector< moris::Matrix< moris::IdMat > >   mVertexSharingData;    // for a processor (single cell in vector) , Vertex Ids and the vertex index on the other processor. Sorted by vertex ids ascending order
            Vector< moris::Matrix< moris::IdMat > >   mCellSharingData;
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_DATA_STK_HPP_ */
