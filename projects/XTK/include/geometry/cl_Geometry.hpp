/*
 * cl_Functon.hpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#ifndef INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_
#define INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_


#include "assert/fn_xtk_assert.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

namespace xtk
{



/*
 * Any geometry must inherit either the discrete or analytic functions interfaces provided in include/geometry
 * For examples on creating analytic or discrete geometries see examples/cl_Geometry_Examples.cpp.
 */

template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Geometry
{
public:

    ~Geometry()
    {
    }


    /*
     * **************************************************************************************
     *              REQUIRED INTERFACE FOR BOTH ANALYTIC AND DISCRETE GEOMETRY
     * **************************************************************************************
     */

    /**
     * This tells the user whether the given geometry is an analytical expression
     * @param[out] True if analytic geometry (i.e. Sphere from an equation)
     */
    virtual bool is_analytic() const = 0;

    /**
     *
     */
    virtual void get_dphi_dp_size(Integer & aNumRows, Integer & aNumCols) const  = 0;


    /*
     * **************************************************************************************
     *                             ANALYTIC GEOMETRY INTERFACE
     * **************************************************************************************
     */

    /*
     * Given a node coordinate, the geometry needs to return the distance to the nearest function.
     */
    virtual Real evaluate_field_value_with_coordinate(Integer const & aRowIndex,
                                                      moris::Matrix< Real_Matrix > const & aCoordinates) const
    {
        XTK_ERROR<<"evaluate_field_value_with_coordinate not implemented. This could be due to a geometry not being based on an analytic expression.";
        return 0;
    }

    /*
     * Given a node coordinate @param[in] aCoordinates, the function returns a matrix of relevant node coordinates
     * Where each cell represents a design variable x,y,z vector
     */
    virtual moris::Matrix< Real_Matrix > evaluate_sensitivity_dphi_dp_with_coordinate(Integer const & aRowIndex,
                                                                                moris::Matrix< Real_Matrix > const & aCoordinates) const
    {
        XTK_ERROR<<"evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on an analytic expression.";
        return moris::Matrix< Real_Matrix >();
    }


    /*
     * **************************************************************************************
     *                             DISCRETE GEOMETRY INTERFACE
     * **************************************************************************************
     */

    /*
     * Given a node index and rank, provide the field value or minimum distance to a geometry feature
     * NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
     */
    virtual Real access_field_value_with_entity_index(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        XTK_ERROR<<" access_field_value_with_entity_id not implemented. This could be due to a geometry not being based on a mesh.";
        return 0;
    }


    /*
     * Given a local coordinate, node index and entity rank, the function returns a matrix of relevant node coordinates
     *  NOTE: THIS CURRENTLY REQUIRES THE BACKGROUND MESH AND GEOMETRY ARE COINCIDENT (NEEDS AN OBJECT THAT MAPS BETWEEN THE TWO)
     */
    virtual moris::Matrix< Real_Matrix > evaluate_sensitivity_dphi_dp(moris::Matrix< Real_Matrix > const & aLocalCoordinate, Integer aEntityIndex, enum EntityRank aEntityRank)
    {
//        XTK_ERROR<<" evaluate_sensitivity_dx_dp not implemented. This could be due to a geometry not being based on a mesh.";
        return moris::Matrix< Real_Matrix >(1,1,0);
    }

    /*
     * Returns the advs that the provided node indices are dependent on
     */
     virtual moris::Matrix< Integer_Matrix > get_node_adv_indices(moris::Matrix< Integer_Matrix > const & aNodeIndices)
    {
         XTK_ERROR<<" get_node_adv_indices not implemented, This could be due to a geometry not being based on a mesh.";
         return moris::Matrix< Integer_Matrix >(1,1,0);
    }



};
}




#endif /* INCLUDE_GEOMETRY_CL_GEOMETRY_HPP_ */
