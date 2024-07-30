/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_WRK_DataBase_Performer.cpp  
 * 
 */

#include "cl_WRK_DataBase_Performer.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_MTK_Interpolation_Mesh_Editor.hpp"
#include "cl_MTK_Mesh_DataBase_IP.hpp"

#include "cl_MTK_Integration_Mesh_Editor.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"

#include "cl_MIG_Mesh_Editor.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace wrk
    {
        //------------------------------------------------------------------------------

        DataBase_Performer::DataBase_Performer( std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer) :
            mMTKInputPerformer( std::move( tMTKPerformer ) )
        {
        }

        //------------------------------------------------------------------------------

        DataBase_Performer::~DataBase_Performer()
        {
            delete mIGMeshEditor;
        }

        //------------------------------------------------------------------------------

        void
        DataBase_Performer::perform()
        {
            Tracer tTracer("DataBase", "No Type", "Perform");
            
            // get the ip and ig meshes from the mesh manager
            mtk::Interpolation_Mesh* tEnrInterpMesh = mMTKInputPerformer->get_mesh_pair( 0 ).get_interpolation_mesh();
            mtk::Integration_Mesh*   tEnrIntegMesh  = mMTKInputPerformer->get_mesh_pair( 0 ).get_integration_mesh();

            mtk::Interpolation_Mesh_Editor tIPMeshEditor = mtk::Interpolation_Mesh_Editor( *tEnrInterpMesh );

            mInterpolationMesh = tIPMeshEditor.perform();

            mIGMeshEditor = new mtk::Integration_Mesh_Editor( tEnrIntegMesh, mInterpolationMesh, mCheckMesh );

            mIntegrationMesh = mIGMeshEditor->perform();

            // assign a name to the mesh
            std::string tXTKMeshName = "DataBase_Meshes";

            // register the mesh pair and grant ownership of the pointers created
            mMTKOutputPerformer->register_mesh_pair( mInterpolationMesh, mIntegrationMesh, true, tXTKMeshName );
        }

        //------------------------------------------------------------------------------

        void
        DataBase_Performer::set_output_performer( std::shared_ptr< mtk::Mesh_Manager > tMTKOutPerformer )
        {
            mMTKOutputPerformer = std::move( tMTKOutPerformer );
        }

        //------------------------------------------------------------------------------

        void
        DataBase_Performer::free_memory()
        {
            mIGMeshEditor->free_memory();
        }

    }// namespace wrk
}// namespace moris
