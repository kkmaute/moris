/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Contact_Sandbox.hpp
 *
 */

// #ifndef SRC_XTK_cl_XTK_Contact_Sandbox_
// #define SRC_XTK_cl_XTK_Contact_Sandbox_

// #include "cl_MTK_Integration_Mesh.hpp"
// #include "cl_MTK_Set.hpp"
// #include "cl_Vector.hpp"
// #include "moris_typedefs.hpp"
// #include "cl_Matrix.hpp"
// #include "linalg_typedefs.hpp"
// #include "fn_Parsing_Tools.hpp"
// #include <ArborX.hpp>
// #include <Kokkos_Core.hpp>
// #include "cl_Stopwatch.hpp" //CHR/src

// namespace moris::xtk
// {
//     using ExecSpace = Kokkos::DefaultExecutionSpace;
//     using MemSpace = typename ExecSpace::memory_space;
//     using DeviceType = Kokkos::Device<ExecSpace, MemSpace>;

//     class Bounding_Box
//     {
//         private:
//         ArborX::Box mBox;
//         moris::uint mSpatialDim;

//         public:
//         Bounding_Box(moris::uint const & aSpatialDim):
//         mSpatialDim(aSpatialDim)
//         {

//         };

//         Bounding_Box():
//         mSpatialDim(3)
//         {

//         };

//         void
//         inflate(moris::real const & aInflateValue)
//         {
//             ArborX::Point & tMinPoint = mBox.minCorner();
//             ArborX::Point & tMaxPoint = mBox.maxCorner();

//             for(moris::uint i = 0; i < mSpatialDim; i++)
//             {
//                 tMinPoint[i] += -aInflateValue;
//                 tMaxPoint[i] +=  aInflateValue;
//             }
//         }

//         void
//         set_bounds(Matrix<DDRMat> const & aBounds)
//         {
//             ArborX::Point & tMinPoint = mBox.minCorner();
//             tMinPoint[0] = aBounds(0,0);
//             tMinPoint[1] = aBounds(1,0);

//             ArborX::Point & tMaxPoint = mBox.maxCorner();
//             tMaxPoint[0] = aBounds(0,1);
//             tMaxPoint[1] = aBounds(1,1);

//             if(mSpatialDim == 2)
//             {
//                 tMinPoint[2] = 0.0;
//                 tMaxPoint[2] = 0.0;
//             }
//             else
//             {
//                 tMinPoint[2] = aBounds(2,0);
//                 tMaxPoint[2] = aBounds(2,1);
//             }

//         }

//         Matrix<DDRMat>
//         get_bounding_box_as_cell() const
//         {
//             ArborX::Point const & tMinPoint = mBox.minCorner();
//             ArborX::Point const & tMaxPoint = mBox.maxCorner();

//             if(mSpatialDim == 2)
//             {
//                 return { { tMinPoint[0], tMinPoint[1]},
//                          { tMaxPoint[0], tMinPoint[1]},
//                          { tMaxPoint[0], tMaxPoint[1]},
//                          { tMinPoint[0], tMaxPoint[1]}};
//             }
//             else
//             {
//                 return {
//                          { tMinPoint[0], tMinPoint[1], tMinPoint[2]},
//                          { tMaxPoint[0], tMinPoint[1], tMinPoint[2]},
//                          { tMaxPoint[0], tMaxPoint[1], tMinPoint[2]},
//                          { tMinPoint[0], tMaxPoint[1], tMinPoint[2]},
//                          { tMinPoint[0], tMinPoint[1], tMaxPoint[2]},
//                          { tMaxPoint[0], tMinPoint[1], tMaxPoint[2]},
//                          { tMaxPoint[0], tMaxPoint[1], tMaxPoint[2]},
//                          { tMinPoint[0], tMaxPoint[1], tMaxPoint[2]} };
//             }
//         }

//         ArborX::Box const &
//         get_arborx_box() const
//         {
//             return  mBox;
//         }

//         ArborX::Box &
//         get_arborx_box()
//         {
//             return  mBox;
//         }

//         ArborX::Point &
//         get_max_point()
//         {
//             return mBox.maxCorner();
//         }

//         ArborX::Point &
//         get_min_point()
//         {
//             return mBox.minCorner();
//         }

//     };

//     // All the sets related to a given mtk side set
//     struct Mesh_Set_AABBs
//     {
//         // B_n = {b_n^1, b_n^2,..., b_n^k }
//         Bounding_Box mAllVertsAABB;
//         Vector<Bounding_Box> mIGVertexAABBs;

//         // B_f = {b_f^1, b_f^2,..., b_f^k }
//         Bounding_Box mAllFacetClusterAABBs;
//         Vector<Bounding_Box> mFacetClusterAABBs;

//     };

// }

// namespace ArborX
// {
// namespace Traits
// {
// template <>
// struct Access<Vector<xtk::Bounding_Box>, PrimitivesTag>
// {
//   inline static std::size_t size(Vector<xtk::Bounding_Box> const &boxes) { return boxes.size(); }
//   KOKKOS_FUNCTION static ArborX::Box const &  get(Vector<xtk::Bounding_Box>  const & boxes, std::size_t i)
//   {
//       return boxes(i).get_arborx_box();
//   }
//   using memory_space = Kokkos::HostSpace;
// };

// template <>
// struct Access<Vector<xtk::Bounding_Box>, PredicatesTag>
// {
//   inline static std::size_t size(Vector<xtk::Bounding_Box> const &boxes) { return boxes.size(); }

//   KOKKOS_INLINE_FUNCTION static auto get(Vector<xtk::Bounding_Box> const &boxes, std::size_t i)
//   {
//     return intersects(boxes(i).get_arborx_box());
//   }
//   using memory_space = Kokkos::HostSpace;
// };
// }
// }
// namespace moris::xtk
// {

//     class Contact_Sandbox
//     {
//     private:
//         mtk::Integration_Mesh * mIntegrationMesh;
//         Matrix<IndexMat> mContactPhasesPairs;
//         moris::real mBBEpsilon;
//         std::string mContactSet0;
//         std::string mContactSet1;
//     public:
//         /*!
//          *
//          */
//         Contact_Sandbox(mtk::Integration_Mesh  * aIntegrationMesh,
//                         std::string              aContactSet0,
//                         std::string              aContactSet1,
//                         moris::real      const & aBBEpsilon):
//                         mIntegrationMesh(aIntegrationMesh),
//                         mContactSet0(aContactSet0),
//                         mContactSet1(aContactSet1),
//                         mBBEpsilon(aBBEpsilon)
//                         {};

//         ~Contact_Sandbox(){};

//         void
//         perform_global_contact_search( Matrix<DDRMat> const & aCurrentDispVec,
//                                        Matrix<DDRMat> const & aPredictedDispVec)
//         {
//                         // This global search detects which bounding boxes of set 1 intersect bounding boxes of set 2
//             // Based on Ghosting for contact global search (hansen et al)
//             // Bounding boxes on proc

//             // get the sets from the mesh
//             moris::mtk::Set* tSet0 = mIntegrationMesh->get_set_by_name(mContactSet0);
//             moris::mtk::Set* tSet1 = mIntegrationMesh->get_set_by_name(mContactSet1);

//             // bounding boxes for mtk sets
//             Mesh_Set_AABBs tMeshSet0AABBs;
//             this->construct_set_AABB(tSet0,aCurrentDispVec,aPredictedDispVec,tMeshSet0AABBs);

//             // inflate the bounding boxes
//             this->inflate_set_AABBS(mBBEpsilon,tMeshSet0AABBs);

//             // visualize
//             this->viz_mtk_set_AABBs(mContactSet0,tMeshSet0AABBs);

//             Mesh_Set_AABBs tMeshSet1AABBs;
//             this->construct_set_AABB(tSet1,aCurrentDispVec,aPredictedDispVec,tMeshSet1AABBs);

//             // inflate the bounding boxes
//             this->inflate_set_AABBS(mBBEpsilon,tMeshSet1AABBs);

//             // visualize
//             this->viz_mtk_set_AABBs(mContactSet1,tMeshSet1AABBs);

//             // Compute communication partners and ghosting candidates
//             // TODO:

//             tic tTimer2;
//             // construct the BVH for vertices
//             ArborX::BVH<DeviceType> tBVH0(tMeshSet0AABBs.mFacetClusterAABBs);
//             ArborX::BVH<DeviceType> tBVH1(tMeshSet1AABBs.mIGVertexAABBs);

//             Kokkos::View<int*, DeviceType> tIndices("indices",0);
//             Kokkos::View<int*, DeviceType> tOffset("offset",0);

//             tBVH0.query(tMeshSet1AABBs.mIGVertexAABBs, tIndices, tOffset);

//             std::string tOverlapFile = "overlap_viz.vtk";
//             this->viz_overlap(tOverlapFile, tMeshSet0AABBs.mFacetClusterAABBs, tMeshSet1AABBs.mIGVertexAABBs, tIndices, tOffset);

//                             // print output
//                             // stop timer
//             real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;
//             std::cout<<"Contact Search Time: "<<( double ) tElapsedTime / 1000<<std::endl;
//             // MORIS_LOG_INFO( "Global contact search:",           ( double ) tElapsedTime / 1000 );

//         }

//         /*!
//         * Construct axis aligned bounding box (AABB) for the set provided
//         * @param[in] aSet Set to construct AABB for
//         * @param[in] aCurrentDispVec Current displacement vector( u_n )
//         * @param[in] aPredictedDispVec Current displacement vector( u_n+1 )
//         * @param[out] aSetAABBs MTK set relevant bounding boxes
//         */
//             void
//             construct_set_AABB(moris::mtk::Set *const &aSet,
//                                Matrix<DDRMat> const &aCurrentDispVec,
//                                Matrix<DDRMat> const &aPredictedDispVec,
//                                Mesh_Set_AABBs &aSetAABBs)
//             {
//                 // construct the vertex bounding boxes
//                 this->construct_set_ig_vertex_AABBs(aSet, aCurrentDispVec, aPredictedDispVec, aSetAABBs.mIGVertexAABBs);

//                 // construct the side cluster bounding boxes
//                 this->construct_side_cluster_AABB(aSet, aCurrentDispVec, aPredictedDispVec, aSetAABBs.mFacetClusterAABBs);

//                 // all vertex bounding box is just the bounding box of all the vertex bounding boxes
//                 aSetAABBs.mAllVertsAABB = this->bounding_box_of_bounding_box(aSetAABBs.mIGVertexAABBs);

//                 // all facet cluster bounding box is just the bounding box of all the facet cluster bounding boxes
//                 aSetAABBs.mAllFacetClusterAABBs = this->bounding_box_of_bounding_box(aSetAABBs.mFacetClusterAABBs);
//         }

//         void
//         construct_side_cluster_AABB(
//                     moris::mtk::Set*   const & aSet,
//                     Matrix<DDRMat>     const & aCurrentDispVec,
//                     Matrix<DDRMat>     const & aPredictedDispVec,
//                     Vector<Bounding_Box>       & aSideClusterBBs)
//         {
//             // number of clusters
//             moris::uint tNumClusters = aSet->get_num_clusters_on_set();

//             // size the number of clusters
//             aSideClusterBBs.resize(tNumClusters, mIntegrationMesh->get_spatial_dim());

//             // get the set clusters
//             Vector<mtk::Cluster const *> tClusters = aSet->get_clusters_on_set();

//             for(moris::uint i = 0; i < tNumClusters; i++)
//             {
//                   Matrix<IndexMat> tVerticesInCluster = tClusters(i)->get_primary_vertices_inds_in_cluster();\

//                   this->construct_AABB_from_vertices(tVerticesInCluster,aCurrentDispVec,aPredictedDispVec,aSideClusterBBs(i));
//             }
//         }

//         void
//         construct_set_ig_vertex_AABBs(
//                     moris::mtk::Set*   const & aSet,
//                     Matrix<DDRMat>     const & aCurrentDispVec,
//                     Matrix<DDRMat>     const & aPredictedDispVec,
//                     Vector<Bounding_Box>       & aVertexBBs)
//         {
//             // get the set clusters
//             Matrix< IndexMat > tSetIgVerts = aSet->get_ig_vertices_inds_on_block(true);

//             // number of clusters
//             moris::uint tNumVerts  = tSetIgVerts.numel();

//             // size the number of clusters
//             aVertexBBs.resize(tNumVerts, mIntegrationMesh->get_spatial_dim());

//             for(moris::uint i = 0; i < tNumVerts; i++)
//             {
//                   this->construct_single_vertex_AABB(tSetIgVerts(i),aCurrentDispVec,aPredictedDispVec,aVertexBBs(i));
//             }
//         }

//         /*!
//         * Construct axis aligned bounding box (AABB) for the point cloud provided
//         * @param[in] aVerticesForBB Vertices to construct BB
//         * @param[in] aCurrentDispVec Current displacement vector( u_n )
//         * @param[in] aPredictedDispVec Current displacement vector( u_n+1 )
//         * @param[out] aVertexSetBB Bounding box for the provided vertices
//         */
//         void
//         construct_AABB_from_vertices(Matrix<IndexMat> const & aVerticesForBB,
//                                      Matrix<DDRMat> const & aCurrentDispVec,
//                                      Matrix<DDRMat> const & aPredictedDispVec,
//                                      Bounding_Box         & aVertexSetBB)
//         {
//             moris::uint tSpatialDim = mIntegrationMesh->get_spatial_dim();

//             moris::uint tNumNodes = aVerticesForBB.numel();

//             // figure out bounding box
//             Vector<real> tMins(tSpatialDim,MORIS_REAL_MAX);
//             Vector<real> tMaxs(tSpatialDim,-MORIS_REAL_MAX);

//             moris::real tCurrentLoc = 0.0;
//             moris::real tPredictedLoc = 0.0;

//             for(moris::moris_index i = 0; i <(moris::moris_index)tNumNodes; i++)
//             {
//                 Matrix< DDRMat > tNodeCoord =  mIntegrationMesh->get_node_coordinate(aVerticesForBB(i));

//                 for(moris::uint iSpatial = 0; iSpatial<tSpatialDim; iSpatial++)
//                 {
//                     tCurrentLoc   = tNodeCoord(iSpatial) + aCurrentDispVec(aVerticesForBB(i),iSpatial);
//                     tPredictedLoc = tNodeCoord(iSpatial) + aPredictedDispVec(aVerticesForBB(i),iSpatial);

//                     if(tCurrentLoc >= tPredictedLoc)
//                     {
//                         if(tCurrentLoc > tMaxs(iSpatial))
//                         {
//                             tMaxs(iSpatial) = tCurrentLoc;
//                         }
//                         if(tPredictedLoc < tMins(iSpatial))
//                         {
//                             tMins(iSpatial) = tPredictedLoc;
//                         }
//                     }
//                     // case where predicted is higher than the current
//                     else
//                     {
//                         if(tPredictedLoc > tMaxs(iSpatial))
//                         {
//                             tMaxs(iSpatial) = tPredictedLoc;
//                         }
//                         if(tCurrentLoc < tMins(iSpatial))
//                         {
//                             tMins(iSpatial) = tCurrentLoc;
//                         }
//                     }
//                 }
//             }

//             // Place in matrix and set bounding box
//             Matrix<DDRMat> tBounds(tSpatialDim,2);
//             for(moris::uint i = 0; i < tSpatialDim; i++)
//             {
//                 tBounds(i,0) = tMins(i);
//                 tBounds(i,1) = tMaxs(i);
//             }

//             aVertexSetBB.set_bounds(tBounds);

//         }

//         void
//         construct_single_vertex_AABB(moris_index    const & aVertexIndex,
//                                      Matrix<DDRMat> const & aCurrentDispVec,
//                                      Matrix<DDRMat> const & aPredictedDispVec,
//                                      Bounding_Box         & aVertexSetBB)
//             {
//                 moris::uint tSpatialDim = mIntegrationMesh->get_spatial_dim();
//                 Matrix< DDRMat > tNodeCoord =  mIntegrationMesh->get_node_coordinate(aVertexIndex);

//                 // allocate bounds
//                 Matrix<DDRMat> tBounds(tSpatialDim,2);

//                 moris::real tCurrentLoc   = 0;
//                 moris::real tPredictedLoc = 0;

//                 for(moris::uint iSpatial = 0; iSpatial<tSpatialDim; iSpatial++)
//                 {
//                     tCurrentLoc   = tNodeCoord(iSpatial) + aCurrentDispVec(aVertexIndex,iSpatial);
//                     tPredictedLoc = tNodeCoord(iSpatial) + aPredictedDispVec(aVertexIndex,iSpatial);

//                     if(tCurrentLoc >= tPredictedLoc)
//                     {
//                         tBounds(iSpatial,0) = tPredictedLoc;
//                         tBounds(iSpatial,1) = tCurrentLoc;
//                     }
//                     else
//                     {
//                         tBounds(iSpatial,0) = tCurrentLoc;
//                         tBounds(iSpatial,1) = tPredictedLoc;
//                     }
//                 }

//                 aVertexSetBB.set_bounds(tBounds);
//         }
//         void
//         inflate_set_AABBS(moris::real const & aInflateVal,
//                           Mesh_Set_AABBs & aMeshSetAABBs)
//         {

//             aMeshSetAABBs.mAllVertsAABB.inflate(aInflateVal);
//             this->inflate_aabbs(aInflateVal,aMeshSetAABBs.mIGVertexAABBs);

//             aMeshSetAABBs.mAllFacetClusterAABBs.inflate(aInflateVal);
//             this->inflate_aabbs(aInflateVal,aMeshSetAABBs.mFacetClusterAABBs);
//         }

//         void
//         inflate_aabbs( moris::real const & aInflateVal,
//                        Vector<Bounding_Box> & aAABBs)
//         {
//             for(moris::uint i = 0; i < aAABBs.size(); i++)
//             {
//                 aAABBs(i).inflate(aInflateVal);
//             }
//         }

//         Bounding_Box
//         bounding_box_of_bounding_box(Vector<Bounding_Box> & aBoundingBoxes)
//         {

//             moris::uint tSpatialDim = mIntegrationMesh->get_spatial_dim();
//             Bounding_Box tBoundingBox( tSpatialDim );
//             ArborX::Point & tMinPoint = tBoundingBox.get_min_point();
//             ArborX::Point & tMaxPoint = tBoundingBox.get_max_point();

//             for(moris::uint iSpatial = 0; iSpatial<tSpatialDim; iSpatial++)
//             {
//                 tMaxPoint[iSpatial] = -MORIS_REAL_MAX;
//                 tMinPoint[iSpatial] =  MORIS_REAL_MAX;
//             }

//             if(tSpatialDim == 2)
//             {
//                 tMinPoint[2] = 0.0;
//                 tMaxPoint[2] = 0.0;
//             }

//             for(moris::uint iBB =0; iBB< aBoundingBoxes.size(); iBB++)
//             {
//                 for(moris::uint iSpatial = 0; iSpatial<tSpatialDim; iSpatial++)
//                 {
//                     if(aBoundingBoxes(iBB).get_min_point()[iSpatial] < tMinPoint[iSpatial])
//                     {
//                         tMinPoint[iSpatial] = aBoundingBoxes(iBB).get_min_point()[iSpatial];
//                     }
//                     if(aBoundingBoxes(iBB).get_max_point()[iSpatial] > tMaxPoint[iSpatial] )
//                     {
//                         tMaxPoint[iSpatial] = aBoundingBoxes(iBB).get_max_point()[iSpatial];
//                     }
//                 }
//             }

//             return tBoundingBox;
//         }

//         void
//         viz_mtk_set_AABBs(const std::string & aFileBase,
//                           Mesh_Set_AABBs & aSetAABBs)
//         {
//             std::string tAllVertexName  = aFileBase + "_all_vert_AABB.vtk";
//             std::string tVertexName     = aFileBase + "_vert_AABB.vtk";
//             std::string tAllClusterName = aFileBase + "_all_clust_AABB.vtk";
//             std::string tClusterName    = aFileBase + "_clust_AABB.vtk";

//             Vector<Bounding_Box const *> tAllVertsBB(1);
//             tAllVertsBB(0) = &aSetAABBs.mAllVertsAABB;

//             Vector<Bounding_Box const *> tAllClusterBB(1);
//             tAllClusterBB(0) =  &aSetAABBs.mAllFacetClusterAABBs;

//             Vector<Bounding_Box const *> tIGVertexAABBs(aSetAABBs.mIGVertexAABBs.size());
//             for(moris::uint i = 0; i < aSetAABBs.mIGVertexAABBs.size() ; i++)
//             {
//                 tIGVertexAABBs(i) = &aSetAABBs.mIGVertexAABBs(i);
//             }

//             Vector<Bounding_Box const *> tFacetClusterAABBs(aSetAABBs.mFacetClusterAABBs.size());
//             for(moris::uint i = 0; i < aSetAABBs.mFacetClusterAABBs.size() ; i++)
//             {
//                 tFacetClusterAABBs(i) = &aSetAABBs.mFacetClusterAABBs(i);
//             }

//             this->save_to_vtk( tAllVertexName,tAllVertsBB );
//             this->save_to_vtk( tVertexName, tIGVertexAABBs);
//             this->save_to_vtk( tAllClusterName, tAllClusterBB);
//             this->save_to_vtk( tClusterName, tFacetClusterAABBs);
//         }

//         void
//         viz_overlap(
//                 const std::string  & aFilePath,
//                 Vector<Bounding_Box> & aBVHBoxes,
//                 Vector<Bounding_Box> & aQueryBoxes,
//                 Kokkos::View<int*, DeviceType> & aIndices,
//                 Kokkos::View<int*, DeviceType> & aOffset)
//         {
//             Vector<Bounding_Box const *> tBoxToViz;

//             // create kokkos::view around std::vector
//             auto tHostViewIndices = Kokkos::create_mirror_view(aIndices);
//             auto tHostViewOffsets = Kokkos::create_mirror_view(aOffset);

//             for(uint iQB = 0; iQB < aQueryBoxes.size(); iQB++)
//             {
//                 if(tHostViewOffsets(iQB+1) - tHostViewOffsets(iQB) > 0)
//                 {
//                     std::cout<<"Found one"<<std::endl;
//                     tBoxToViz.push_back(&aQueryBoxes(iQB));
//                 }
//              }

//             this->save_to_vtk(aFilePath,tBoxToViz);
//         }
//         /*
//         * Save the cell of bounding boxes to a VTK file
//         */
//         void
//         save_to_vtk(const std::string   & aFilePath,
//                     Vector<Bounding_Box const *> & aBoxesToViz)
//         {
//             // modify filename
//             std::string tFilePath = parallelize_path( aFilePath );

//             // open the file
//             std::ofstream tFile(tFilePath, std::ios::binary);

//             // containers
//             float tFChar = 0;
//             int   tIChar = 0;

//             tFile << "# vtk DataFile Version 3.0" << std::endl;
//             tFile << "GO BUFFS!" << std::endl;
//             tFile << "BINARY" << std::endl;

//             // initialize element counter
//             luint tNumberOfElements = aBoxesToViz.size();

//             // number of nodes per element
//             uint tNumberOfNodesPerElement = std::pow( 2, mIntegrationMesh->get_spatial_dim() );

//             // count number of nodes
//             uint tNumberOfNodes = tNumberOfNodesPerElement*tNumberOfElements;

//             // write node data
//             tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

//             tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;

//             // temporary matrix containing corder nodes
//             Matrix< DDRMat > tNodes( mIntegrationMesh->get_spatial_dim(), tNumberOfNodesPerElement );

//             std::cout<<"tNumberOfNodesPerElement = "<<tNumberOfNodesPerElement<<std::endl;

//             // VTK cell type
//             int tCellType = 0;

//             if ( mIntegrationMesh->get_spatial_dim() == 2 )
//             {
//                 // loop over all elements
//                 for( auto tElement: aBoxesToViz )
//                 {
//                         // ask background mesh for corner nodes
//                         Matrix<DDRMat> tBoundCoord = tElement->get_bounding_box_as_cell();

//                         // write node coordinates to file
//                         for( uint k=0; k<tNumberOfNodesPerElement; ++k )
//                         {
//                             tFChar = swap_byte_endian( (float) tBoundCoord( k, 0) );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                             tFChar = swap_byte_endian( (float) tBoundCoord( k, 1 )  );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                             tFChar = swap_byte_endian( (float) 0 );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                         }
//                 }

//                 // VTK celltype for Quad4
//                 tCellType = 9;
//             }
//             else if ( mIntegrationMesh->get_spatial_dim() == 3 )
//             {
//                // loop over all elements
//                 for( auto tElement: aBoxesToViz )
//                 {
//                         // ask background mesh for corner nodes
//                         Matrix<DDRMat> tBoundCoord = tElement->get_bounding_box_as_cell();

//                         // write node coordinates to file
//                         for( uint k=0; k<tNumberOfNodesPerElement; ++k )
//                         {
//                              tFChar = swap_byte_endian( (float) tBoundCoord( k, 0 ) );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                             tFChar = swap_byte_endian( (float) tBoundCoord( k, 1 ) );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                             tFChar = swap_byte_endian( (float)tBoundCoord( k, 2 ) );
//                             tFile.write( (char*) &tFChar, sizeof(float));
//                         }
//                 }

//                 // VTK celltype for Hex8
//                 tCellType = 12;
//             }

//             // create new line
//             tFile << std::endl;

//             // can only write element data if vtk map exists
//             if ( tCellType != 0 )
//             {
//                 // value to write in VTK file
//                 int tNumberOfNodesVTK = swap_byte_endian( (int) tNumberOfNodesPerElement );

//                 // reset node counter
//                 int tCount = 0;

//                 // write header for cells
//                 tFile << "CELLS " << tNumberOfElements << " "
//                         << ( tNumberOfNodesPerElement + 1 )*tNumberOfElements  << std::endl;

//                 // loop over all elements
//                 for( uint iBB = 0 ; iBB<aBoxesToViz.size(); iBB++ )
//                 {

//                         tFile.write( (char*) &tNumberOfNodesVTK, sizeof(int) );

//                         // loop over all nodes of this element
//                         for( uint k=0; k<tNumberOfNodesPerElement; ++k )
//                         {
//                             // write node to mesh file
//                             tIChar = swap_byte_endian( tCount++ );
//                             tFile.write((char *) &tIChar, sizeof(int));
//                         }

//                 }

//                 // write cell types
//                 tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
//                 tIChar = swap_byte_endian( tCellType );
//                 for ( luint k = 0; k < tNumberOfElements; ++k)
//                 {
//                     tFile.write( (char*) &tIChar, sizeof(int));
//                 }

//                 // close the output file
//                 tFile.close();

//             }

//         }

//     };

// }

// #endif
