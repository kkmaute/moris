/*
 * cl_XTK_Interpolation_Vertex_Unzipped.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_


#include "typedefs.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
namespace mtk
{
class Vertex_Interpolation;
}
}


using namespace moris;

namespace xtk
{
class Interpolation_Vertex_Unzipped: public mtk::Vertex
{
public:
    Interpolation_Vertex_Unzipped(){};

    Interpolation_Vertex_Unzipped(mtk::Vertex*               aBaseInterpVertex,
                                  moris_id                   aVertexId,
                                  moris_index                aVertexIndex,
                                  moris_index                aVertexOwner,
                                  uint                       aInterpolationOrder,
                                  mtk::Vertex_Interpolation* aVertexInterp);
    //------------------------------------------------------------------------------
    Matrix< DDRMat >
    get_coords() const;
    //------------------------------------------------------------------------------
    moris_id
    get_id() const;
    //------------------------------------------------------------------------------
    moris_index
    get_index() const;
    //------------------------------------------------------------------------------
    moris_index
    get_owner() const;
    //------------------------------------------------------------------------------
    mtk::Vertex_Interpolation *
    get_interpolation( const uint aOrder );
    //------------------------------------------------------------------------------
    const mtk::Vertex_Interpolation *
    get_interpolation( const uint aOrder ) const;
private:
    // the parent vertex
    mtk::Vertex* mBaseInterpVertex;

    moris_id    mVertexId;
    moris_index mVertexIndex;
    moris_index mVertexOwner;

    // interpolation of vertex
    uint mOrder; /* Order of the Vertex interpolation*/
    mtk::Vertex_Interpolation* mInterpolation;

};
}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_ */
