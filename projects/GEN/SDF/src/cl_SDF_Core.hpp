/*
 * cl_SDF_Core.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_


#include "typedefs.hpp" // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Triangle_Mesh.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Data.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class Core
        {
                  Mesh          & mMesh;
                  Data          & mData;

                  uint            mCandidateSearchDepth = 1;
                  real            mCandidateSearchDepthEpsilon = 0.01;
                  bool            mVerbose;

//-------------------------------------------------------------------------------
        public :
//-------------------------------------------------------------------------------

            Core(       Mesh & aMesh,
                        Data & aData,
                        bool   aVerbose=false );

//-------------------------------------------------------------------------------

            ~Core(){};

//-------------------------------------------------------------------------------

            void
            calculate_raycast();

//-------------------------------------------------------------------------------

            void
            calculate_raycast(
                    Matrix< IndexMat > & aElementsAtSurface );

//-------------------------------------------------------------------------------
            void
            calculate_raycast(
                    Matrix< IndexMat > & aElementsAtSurface,
                    Matrix< IndexMat > & aElementsInVolume );

//-------------------------------------------------------------------------------

            void
            calculate_raycast_and_sdf( Matrix< DDRMat> & aSDF );

//-------------------------------------------------------------------------------

            void
            save_to_vtk( const std::string & aFilePath );

//-------------------------------------------------------------------------------
        private :
//-------------------------------------------------------------------------------


//-------------------------------------------------------------------------------

            void
            voxelize( const uint aAxis );


//-------------------------------------------------------------------------------

            void
            calculate_udf();

//-------------------------------------------------------------------------------

            void
            sweep();

//-------------------------------------------------------------------------------


            void
             fill_sdf_with_values( Matrix< DDRMat > & aSDF );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_x( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_y( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_z( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            intersect_triangles(
                    const uint aAxis,
                    const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            intersect_ray_with_triangles(
                    const uint aAxis,
                    const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            check_if_node_is_inside(
                    const uint aAxis,
                    const uint aNodeIndex );
//-------------------------------------------------------------------------------

            void
            calculate_candidate_points_and_buffer_diagonal();

//-------------------------------------------------------------------------------

            void
            get_nodes_withing_bounding_box_of_triangle(
                            Triangle * aTriangle, moris::Cell< Vertex* > & aNodes );

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_ */
