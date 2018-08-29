/*
 * cl_Mesh_Tools.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_TOOLS_HPP_
#define SRC_MESH_CL_MESH_TOOLS_HPP_

#include"mesh/cl_Mesh_Data.hpp"
#include"mesh/cl_Mesh_Enums.hpp"

#include"tools/fn_tet_volume.hpp"
#include"linalg/cl_XTK_Matrix.hpp"




using namespace xtk;

namespace mesh
{
class Mesh_Helper
{
public:
    template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
    static moris::Mat_New<Integer, Integer_Matrix>
    get_glb_entity_id_from_entity_loc_index_range(Mesh_Data<Real,Integer, Real_Matrix, Integer_Matrix> const & aMeshData,
                                                  moris::Mat_New<Integer, Integer_Matrix> const & tEntityIndices,
                                                  enum EntityRank aEntityRank)
    {
        Integer tNumEntities = tEntityIndices.n_cols();
        moris::Mat_New<Integer, Integer_Matrix> tEntityIds(1,tNumEntities);

        for(Integer i =0; i<tNumEntities; i++)
        {
            tEntityIds(0,i) = aMeshData.get_glb_entity_id_from_entity_loc_index(tEntityIndices(0,i),aEntityRank);
        }
        return tEntityIds;
    }


    template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
    static void write_mesh_to_template(Mesh_Data<Real,Integer, Real_Matrix, Integer_Matrix> const & aMeshData,
                                       std::string const & aIndiceFile,
                                       std::string const & aOffsetFile)
    {
        std::ofstream tIndiceOutput;
            tIndiceOutput.open(aIndiceFile);

            std::ofstream tOffsetOutput;
            tOffsetOutput.open(aOffsetFile);

            // Get numbers of entities
            //                Integer tNumNode = get_num_entities(EntityRank::NODE);
            //                Integer tNumEdge = get_num_entities(EntityRank::EDGE);
            //                Integer tNumFace = get_num_entities(EntityRank::FACE);

            moris::Mat_New<Integer,Integer_Matrix> tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT, EntityRank::NODE);
            // 0d1d-------------------------------------------------------
            Integer tNumEnt = aMeshData.get_num_entities(EntityRank::NODE);
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices0d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets0d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            Integer j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::NODE,EntityRank::EDGE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 0d2d-------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices0d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets0d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::NODE,EntityRank::FACE);

                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 0d3d-------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices0d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets0d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::NODE,EntityRank::ELEMENT);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 1d0d ---------------------------------------------------------
            tNumEnt = aMeshData.get_num_entities(EntityRank::EDGE);
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices1d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets1d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::EDGE,EntityRank::NODE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 1d2d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices1d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets1d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::EDGE,EntityRank::FACE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;


            // 1d3d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices1d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets1d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::EDGE,EntityRank::ELEMENT);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 2d0d ---------------------------------------------------------
            tNumEnt = aMeshData.get_num_entities(EntityRank::FACE);
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices2d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets2d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::FACE,EntityRank::NODE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 2d1d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices2d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets2d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::FACE,EntityRank::EDGE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            // 2d3d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices2d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets2d3d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::FACE,EntityRank::ELEMENT);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;


            //3d to 0d ---------------------------------------------------------
            tNumEnt = aMeshData.get_num_entities(EntityRank::ELEMENT);
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices3d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets3d0d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::NODE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            //3d to 1d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices3d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets3d1d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::EDGE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;                       //3d to 2d ---------------------------------------------------------
            tIndiceOutput<< "moris::Mat_New<Integer,Integer_Matrix> tIndices3d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            tOffsetOutput<< "moris::Mat_New<Integer,Integer_Matrix> tOffsets3d2d = mMatrixFactory->create_integer_type_matrix_base( ";
            j = 0;
            tIndiceOutput<<"{{ ";
            tOffsetOutput<<"{{0, ";
            for(Integer i = 0; i<tNumEnt;i++)
            {
                tConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::FACE);
                for(Integer c = 0; c<tConnectivity->n_cols();c++)
                {
                    tIndiceOutput<<(tConnectivity)(0,c);
                    j++;
                    if((c==tConnectivity->n_cols()-1) && (i==tNumEnt-1)) continue;
                    else tIndiceOutput<<", ";
                }
                tOffsetOutput<<j;
                if(i==tNumEnt-1) continue;
                else tOffsetOutput<<", ";
            }

            tIndiceOutput<<"}});"<<std::endl;
            tOffsetOutput<<"}});"<<std::endl;

            tIndiceOutput.close();
            tOffsetOutput.close();
    }

    /*
     * Only works on Tets and Cube Hex 8s (for now) volumes seperated by buckets
     */
    template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
    static void compute_mesh_volume_by_buckets(Mesh_Data<Real,Integer, Real_Matrix, Integer_Matrix> const & aMeshData,
                                               xtk::Cell<Real>  & aVolumes)
    {

        enum EntityRank tElementRank = EntityRank::ELEMENT;
        Integer tNumBuckets = aMeshData.get_num_buckets(EntityRank::ELEMENT);
        aVolumes = xtk::Cell<Real>(tNumBuckets,0);

        moris::Mat_New<Real,Real_Matrix> tNodeCoordinates(0,0);
        moris::Mat_New<Integer, Integer_Matrix>  tElementsInBucket(0,0);
        moris::Mat_New<Integer, Integer_Matrix>  tElementToNodeConnectivity(0,0);

        Real tVol;
        for(Integer iBucket = 0; iBucket <tNumBuckets; iBucket++)
        {
            tElementsInBucket = aMeshData.get_entities_in_bucket_loc_index(iBucket,tElementRank);

            if(aMeshData.get_entity_connected_to_entity_loc_inds((tElementsInBucket)(0,0),EntityRank::ELEMENT,EntityRank::NODE).n_cols() == 4 )
            {
                for(Integer iElem = 0; iElem< tElementsInBucket.n_cols(); iElem++)
                {

                    tElementToNodeConnectivity = aMeshData.get_entity_connected_to_entity_loc_inds((tElementsInBucket)(0,iElem),EntityRank::ELEMENT,EntityRank::NODE);
                    tNodeCoordinates = aMeshData.get_selected_node_coordinates_loc_inds(tElementToNodeConnectivity);
                    tVol = xtk::vol_tetrahedron(tNodeCoordinates);
                    if(tVol<0)
                    {
                        std::cout<<"Warning inverted element! Volume = "<<tVol<<std::endl;
                        tVol = -tVol;
                    }
                    aVolumes(iBucket) += tVol;
                }
            }

            else if(aMeshData.get_entity_connected_to_entity_loc_inds((tElementsInBucket)(0,0),EntityRank::ELEMENT,EntityRank::NODE).n_cols() == 8 )
            {
                for(Integer iElem = 0; iElem< tElementsInBucket.n_cols(); iElem++)
                {
                    Integer tLx = 1;
                    Integer tLy = 1;
                    Integer tLz = 1;
                    aVolumes(iBucket) +=1;
                }
            }
        }


    }
};

}

#endif /* SRC_MESH_CL_MESH_TOOLS_HPP_ */
