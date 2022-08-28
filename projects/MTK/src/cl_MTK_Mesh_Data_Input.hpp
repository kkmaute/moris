/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Data_Input.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Fields_Info.hpp"
#include "cl_MTK_Cell_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"

#include <iomanip>      // std::setw

namespace moris
{
namespace mtk
{
    //////////////////////////
    // STRUC FOR BLOCK SET  //
    //////////////////////////
    struct MtkMeshData
    {
        // Number of spatial dimensions for the mesh
        uint*                           SpatialDim ;

        // element to node connectivity
        // Cell - topology of element, inner matrix  - element to node connectivity
        moris::Cell<Matrix < IdMat >*>  ElemConn;

        // Owning processor of a node
        Matrix < IdMat >*               NodeProcOwner;

        // Processor Ids that share a node (row - node index, col - procs)
        // Should contain MORIS_ID_MAX as a dummy value
        Matrix < IdMat >*               NodeProcsShared;

        // Clustering Inputs
        // Cell clusters
        Cell_Cluster_Input * CellClusterInput = nullptr;

        // Side clusters
        Side_Cluster_Input * SideClusterInput = nullptr;

        // Double Side clusters
        Double_Side_Cluster_Input * DoubleSideClusterInput = nullptr;

        // Node coordinates (row - node index, col - (x,y,z)
        Matrix < DDRMat >*              NodeCoords;
        moris::Cell<Matrix < IdMat >*>  LocaltoGlobalElemMap;
        moris::Cell<enum CellTopology>  CellTopology;
        Matrix < IdMat >*               LocaltoGlobalNodeMap;
        bool                            CreateAllEdgesAndFaces;
        MtkFieldsInfo*                  FieldsInfo;
        MtkSetsInfo*                    SetsInfo;
        real                            TimeStamp = 0.0;
        bool                            AutoAuraOptionInSTK = true;
        bool                            SupplementaryToFile = false;
        bool                            Verbose = false;
        bool                            MarkNoBlockForIO = true;

        MtkMeshData(uint aNumElementTypes):
            SpatialDim(),
            ElemConn(aNumElementTypes),
            NodeProcOwner(),
            NodeProcsShared(nullptr),
            NodeCoords(),
            LocaltoGlobalElemMap(aNumElementTypes),
            CellTopology(aNumElementTypes,CellTopology::INVALID),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(true),
            FieldsInfo(),
            SetsInfo(nullptr)
        {

        }

        MtkMeshData():
            SpatialDim(),
            ElemConn(1),
            NodeProcOwner(),
            NodeProcsShared(nullptr),
            NodeCoords(),
            LocaltoGlobalElemMap(1),
            CellTopology(1,CellTopology::INVALID),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(true),
            FieldsInfo(),
            SetsInfo(nullptr)
        {

        }

        bool
        has_mesh_sets()
        {
            bool tAnswer = true;
            if ( SetsInfo == NULL )
            {
                return false;
            }

            return tAnswer;
        }

        /*
         * count the number of elements provided
         */
        uint
        get_num_elements()
        {
            uint tNumElements = 0;
            for(uint i = 0; i<ElemConn.size(); i++)
            {
                tNumElements += ElemConn(i)->n_rows();
            }
            return tNumElements;
        }

        /*
         * Get the number of nodes
         */
        uint
        get_num_nodes()
        {
            return LocaltoGlobalNodeMap->numel();
        }

        /*
         * Number of elements in the element map
         */
        uint
        size_local_to_global_elem_map()
        {
            uint tSizeLocalToGlobal = 0;
            for(uint i = 0; i <LocaltoGlobalElemMap.size(); i++)
            {
                tSizeLocalToGlobal += LocaltoGlobalElemMap(i)->numel();
            }
            return tSizeLocalToGlobal;
        }

        /*
         * Collapse the element cell to a continuous vector
         */
        Matrix< IdMat >
        collapse_element_map()
        {
            Matrix< IdMat > tCollapsedMap(this->get_num_elements(),1);

            uint tCount = 0;
            for(uint i = 0; i <LocaltoGlobalElemMap.size(); i++)
            {
                for(uint j = 0; j<LocaltoGlobalElemMap(i)->numel(); j++)
                {
                    tCollapsedMap(tCount) = (*LocaltoGlobalElemMap(i))(j);
                    tCount++;
                }
            }

            return tCollapsedMap;
        }

        bool
        has_node_sharing_info()
        {
            if( NodeProcsShared == nullptr)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        void
        print_summary()
        {
            std::cout<<"--------------------------------------------------------"<<std::endl;
            std::cout<<"Mesh Data Output Summary: "<<std::endl;
            std::cout<<"    Processor Rank: "<<par_rank()<<std::endl;
            std::cout<<"    Spatial Dimension: "<<*SpatialDim<<std::endl;
            std::cout<<"    Number of cells: "<<this->get_num_elements()<<std::endl;
            std::cout<<"    Number of vertices: "<<this->get_num_nodes()<<std::endl;
            std::cout<<"--------------------------------------------------------"<<std::endl;
        }

        void
        print_vertex_coordinates()
        {
            std::cout<<"Vertex Coordinates:"<< std::endl;
            for(moris::uint i = 0; i < NodeCoords->n_rows(); i++)
            {
                std::cout<<"    Vert Id: "<<std::right<<std::setw(5)<<(*LocaltoGlobalNodeMap)(i)<<" | ";
                for(moris::uint j = 0; j < *SpatialDim; j++)
                {
                    std::cout<<std::scientific<<std::setw(14)<<(*NodeCoords)(i,j)<< "   ";
                }

                std::cout<<std::endl;
            }
        }

        void
        print_vertex_sharing()
        {
            std::cout<<"Vertex Sharing:"<<std::endl;
            if(NodeProcsShared != nullptr)
            {
                for(moris::uint i = 0; i < NodeProcsShared->n_rows(); i++)
                {
                    std::cout<<"    Vert Id: "<<std::right<<std::setw(5)<<(*LocaltoGlobalNodeMap)(i)<<" | ";
                    for(moris::uint j = 0; j < NodeProcsShared->n_cols(); j++)
                    {
                        if((*NodeProcsShared)(i,j) == MORIS_ID_MAX)
                        {
                            break;
                        }
                        std::cout<<std::scientific<<(*NodeProcsShared)(i,j)<< "   ";
                    }

                    std::cout<<std::endl;
                }
            }
            else
            {
                std::cout<<"No Vertex Sharing provided"<<std::endl;
            }
        }

        void
        print_cell_to_node()
        {
            for(moris::uint  iType = 0; iType < ElemConn.size(); iType++)
            {
                std::cout<<"\nNumber of elements: "<<ElemConn(iType)->n_rows()<<" | Cell Topology: "<<get_enum_str(CellTopology(iType))<<std::endl;

                for(moris::uint i = 0; i < ElemConn(iType)->n_rows(); i++)
                {
                    std::cout<<"    Cell Id: "<<std::right<<std::setw(5)<<(*LocaltoGlobalElemMap(iType))(i)<<" | ";
                    for(moris::uint j = 0; j < ElemConn(iType)->n_cols(); j++)
                    {
                        std::cout<<std::right<<std::setw(10)<<(*ElemConn(iType))(i,j)<< "   ";
                    }

                    std::cout<<std::endl;
                }
                std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            }
        }

        void
        print_sets()
        {
            if(SetsInfo !=nullptr)
            {
                SetsInfo->print();
            }
            else{
                std::cout<<"Cell sets: 0"<<std::endl;
                std::cout<<"Side sets: 0"<<std::endl;
                std::cout<<"Vertex sets: 0"<<std::endl;
            }
        }

        void
        print_fields()
        {
            if(FieldsInfo != nullptr)
            {
                FieldsInfo->print();
            }
            else{
                std::cout<<"No Fields Information Provided"<<std::endl;
            }
        }

        void
        print_cell_clusters()
        {
            if(CellClusterInput != nullptr)
            {
                CellClusterInput->print();
            }
            else
            {
                std::cout<<"No Cell Cluster Information Provided"<<std::endl;
            }
        }

        void
        print_side_cluster()
        {
            if(SideClusterInput != nullptr)
            {
                SideClusterInput->print();
            }
            else
            {
                std::cout<<"No Side Cluster Information Provided"<<std::endl;
            }
        }

        void
        print_full_details()
        {
            print_summary();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_vertex_coordinates();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_vertex_sharing();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_cell_to_node();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_sets();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_fields();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_cell_clusters();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
            print_side_cluster();
            std::cout<<"\n-------------------------------------------------------------------------"<<std::endl;
        }

    };
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_ */

