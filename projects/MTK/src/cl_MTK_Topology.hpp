/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Topology.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_TOPOLOGY_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_TOPOLOGY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
namespace mtk
{
class Interpolation
{
    virtual
    void evaluate_interpolation_function(Matrix< DDRMat > const & aLocalCoordinate,
                                         Matrix< DDRMat > &      aBasisFunctionValues) const = 0;
}
}
}

namespace moris
{
namespace mtk
{
class Interpolation_Hex8: public Interpolation
{
    void evaluate_interpolation_function(moris::Matrix< DDRMat > const & aLocalCoordinate,
                                 moris::Matrix< DDRMat > &      aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,8);

        real tXi = aLocalCoordinate(0,0);
        real tEta = aLocalCoordinate(0,1);
        real tZeta = aLocalCoordinate(0,2);

        aBasisFunctionValues(0,0) = 0.125*(1.0-tXi)*(1.0-tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,1) = 0.125*(1.0+tXi)*(1.0-tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,2) = 0.125*(1.0+tXi)*(1.0+tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,3) = 0.125*(1.0-tXi)*(1.0+tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,4) = 0.125*(1.0-tXi)*(1.0-tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,5) = 0.125*(1.0+tXi)*(1.0-tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,6) = 0.125*(1.0+tXi)*(1.0+tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,7) = 0.125*(1.0-tXi)*(1.0+tEta)*(1.0+tZeta);
    }
};

}
}

namespace moris
{
namespace mtk
{

    Interpolation*
    get_interpolation(enum CellTopology aCellTopology) const
    {

        Interpolation* tBasis = NULL;
        switch(aCellTopology)
        {
            case(CellTopology::HEX8):
            {
                tBasis = new Interpoli
                break;
            }
            default:
                MORIS_ERROR(0,"Invalide Cell Topology");
                break;
        }
    }
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_TOPOLOGY_HPP_ */

