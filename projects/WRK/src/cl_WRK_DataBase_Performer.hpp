/*
 * cl_WRK_DataBase_Performer.hpp
 *
 *  Created on: Dec  29, 2021
 *      Author: momo
 */
#ifndef SRC_cl_WRK_DataBase_Performer
#define SRC_cl_WRK_DataBase_Performer
#include "cl_WRK_Performer.hpp"

#include <memory>

namespace moris
{

    namespace mtk
    {
        class Mesh_Manager;
        class Interpolation_Mesh_DataBase;
        class Integration_Mesh_DataBase;
        class Interpolation_Mesh_DataBase_IP;
        class Integration_Mesh_DataBase_IG;
        class Integration_Mesh_Editor;
    }// namespace mtk

    namespace wrk
    {
        class DataBase_Performer
        {
          private:
            std::shared_ptr< mtk::Mesh_Manager > mMTKInputPerformer; /*!< Input performer that holds the meshes database is going to work on */

            std::shared_ptr< mtk::Mesh_Manager > mMTKOutputPerformer; /*!< Output performer that will hold the meshes database creates */

            mtk::Interpolation_Mesh_DataBase_IP* mInterpolationMesh;
            mtk::Integration_Mesh_DataBase_IG*   mIntegrationMesh;
            mtk::Integration_Mesh_Editor*        mIGMeshEditor;

            // flag allowing the mesh check to be turned on/off for debugging
            bool mCheckMesh = true;

          public:
            //------------------------------------------------------------------------------

            /**
             * @brief Construct a new DataBase_Performer object
             *
             * @param tMTKPerformer
             */
            DataBase_Performer( std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer);

            //------------------------------------------------------------------------------

            /**
             * @brief Destroy the DataBase_Performer object
             *
             */

            ~DataBase_Performer();

            //------------------------------------------------------------------------------

            /**
             * @brief perform function that constrcuts the data base
             *
             */
            void
            perform();

            //------------------------------------------------------------------------------

            /**
             * @brief Set the output performer
             *
             * @param tMTKOutPerformer
             */

            void
            set_output_performer( std::shared_ptr< mtk::Mesh_Manager > tMTKOutPerformer );

            //------------------------------------------------------------------------------

            /**
             * @brief
             *
             */

            void
            free_memory();

            //------------------------------------------------------------------------------

            /**
             * @brief Set whether a mesh check should be performed or not after building the mesh data base
             * 
             * @param aCheckMesh 
             */
            void
            set_mesh_check( bool aCheckMesh )
            {
                mCheckMesh = aCheckMesh;
            }

            //------------------------------------------------------------------------------
        };
    }// namespace wrk
}// namespace moris


#endif /* cl_WRK_DataBase_Performer.hpp */