/*
 * cl_mig_Mesh_Editor.hpp
 *
 *  Created on: Jan  25, 2022
 *      Author: momo
 */
#ifndef SRC_cl_MIG_Mesh_Editor
#define SRC_cl_MIG_Mesh_Editor

#include "typedefs.hpp"
#include <memory>
#include "cl_MTK_Integration_Mesh_Editor.hpp"

namespace moris::ge
{
    class Geometry_Engine;
}

namespace moris::mtk
{
    class Integration_Mesh_DataBase_IG;
}// namespace moris::mtk


namespace moris::mig
{
    class Periodic_2D;
    class Periodic_3D;

    class Periodic_Mesh_Editor : public mtk::Integration_Mesh_Editor
    {
      private:
        mig::Periodic_2D* mPeriodicData2D = nullptr;

        mig::Periodic_3D* mPeriodicData3D = nullptr;

        moris::ge::Geometry_Engine* mGeometryEngine;

      public:
        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief Construct a new Periodic_Mesh_Editor object
         * 
         */

        Periodic_Mesh_Editor() = default;

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief Construct a new Periodic_Mesh_Editor object
         * 
         * @param aIGMesh 
         * @param aPeriodicData2D 
         */

        Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_2D* aPeriodicData2D );

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief Construct a new Periodic_Mesh_Editor object
         * 
         * @param aIGMesh 
         * @param mPeriodicData3D 
         */

        Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_3D* mPeriodicData3D );

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief Destroy the Periodic_Mesh_Editor object
         * 
         */

        ~Periodic_Mesh_Editor();
        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief added newly genrated data to the data base
         * 
         * @param aSideClusterToVertexIndices 
         * @param aVerticesCoords 
         * @param aSideClusterToCells 
         * @param aCellToVertexIndices 
         * @param aSideClusterToIPCell 
         * @param aVertexParametricCoords 
         * @param aDoubleSidedClustersIndex 
         * @param mNumDblSideCluster 
         * @param aNumGeometry 
         */

        virtual void
        construct_periodic_data_base(
            moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices,
            Matrix< DDRMat >                           aVerticesCoords,
            moris::Cell< moris::Cell< moris_index > >& aSideClusterToCells,
            moris::Cell< moris::Cell< moris_index > >& aCellToVertexIndices,
            moris::Cell< moris_index >&                aSideClusterToIPCell,
            Matrix< DDRMat >&                          aVertexParametricCoords,
            moris::Cell< moris_index >&                aDoubleSidedClustersIndex,
            uint                                       mNumDblSideCluster,
            uint                                       aNumGeometry ) override;

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief adds data to the ge and also to the mtk data base
         * 
         */

        void  
        perform();

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief Set the geometry engine object
         * 
         * @param aGeometryEngine 
         */

        void
        set_geometry_engine( moris::ge::Geometry_Engine* aGeometryEngine );

        //------------------------------------------------------------------------------------------------------------

        /**
         * @brief adds vertex data to the geomtry
         * 
         */

        void
        link_nodes_to_geomtry_engine();

        //------------------------------------------------------------------------------------------------------------

        //############################################################################################################/
        //These function are not used, added here for later development
        //############################################################################################################/

        void
        create_vertices();

        //------------------------------------------------------------------------------------------------------------


        void
        create_cells();

        //------------------------------------------------------------------------------------------------------------

        void
        create_side_clusters();

        //------------------------------------------------------------------------------------------------------------

        void
        create_side_sets();

        //------------------------------------------------------------------------------------------------------------

        void
        create_double_sided_clusters();

        //------------------------------------------------------------------------------------------------------------

        void
        create_double_sided_set();

        //------------------------------------------------------------------------------------------------------------

        void
        create_parallel_consistnet_new_vertex_ids();

        //------------------------------------------------------------------------------------------------------------

        void
        reconstrcut_connectivity();

        //------------------------------------------------------------------------------------------------------------

        void
        merge_meshes();
    };
}// namespace moris::mig


#endif /* cl_mig_Mesh_Editor.hpp */