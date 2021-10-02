/*
 * cl_XTK_Diagnostics.hpp  
 * 
 *  Created on: Oct  01, 2021 
 *      Author: Keenan Doble
 */
#ifndef SRC_cl_XTK_Diagnostics
#define SRC_cl_XTK_Diagnostics

#include "cl_Matrix.hpp"

namespace xtk
{
    class Cut_Integration_Mesh;
    class Model;

     /**
    * @brief Collection of diagnostic free functions for XTK
    */

    /**
    * @brief Checks if the interpolation to a coordinate is as expected
    * 
    * @param aCutMesh 
    * @return true - all coordinates of integration vertices evaluate to their expected value
    * @return false 
    */
    bool
    interpolated_coordinate_check(Cut_Integration_Mesh* aCutMesh);
    
    /**
     * @brief Checks that interface vertices are sufficiently close to the interface
     * 
     * @param aModel - XTK model
     * @param aIsocontourThreshold - Isocontour threshold
     * @param aIsocontourTolerance - Isocontour tolerance
     * @return true - all interface vertices are on the interface
     * @return true - all interface vertices are NOT on the interface
     */
    bool
    verify_interface_vertices(
        xtk::Model* aModel, 
        moris::Matrix<moris::DDRMat> const & aIsocontourThreshold,
        moris::Matrix<moris::DDRMat> const & aIsocontourTolerance);
}








#endif /* cl_XTK_Diagnostics.hpp */