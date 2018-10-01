/*
 * cl_SDF_Core.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_

#include <GEN/SDF/src/cl_SDF_Triangle_Mesh.hpp>
#include "typedefs.hpp" // COR/src
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Data.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class Core
        {
            const mtk::Mesh          * mMesh;
                  Data               & mData;

//-------------------------------------------------------------------------------
        public :
//-------------------------------------------------------------------------------

            Core( const mtk::Mesh * aMesh,
                        Data      & aData );

//-------------------------------------------------------------------------------

            ~Core(){};

//-------------------------------------------------------------------------------
        private :
//-------------------------------------------------------------------------------

            void
            voxelize( const uint aAxis );

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
                            const uint                aAxis,
                            const Matrix< F31RMat > & aPoint );

//-------------------------------------------------------------------------------

            void
            check_if_node_is_inside(
                    const uint aAxis,
                    const uint aLocalNodeInd,
                    const Matrix< F31RMat > & aPoint );

//-------------------------------------------------------------------------------

            void
            calculate_candidate_points_and_buffer_diagonal();

//-------------------------------------------------------------------------------
        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_ */
