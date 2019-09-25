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

#include "cl_XTK_Vertex_Enrichment.hpp"

namespace xtk
{
class Vertex_Enrichment;
}


using namespace moris;

namespace xtk
{
class Interpolation_Vertex_Unzipped: public mtk::Vertex
{
public:
    Interpolation_Vertex_Unzipped(){};

    Interpolation_Vertex_Unzipped(mtk::Vertex*       aBaseInterpVertex,
                                  moris_id           aVertexId,
                                  moris_index        aVertexIndex,
                                  moris_index        aVertexOwner,
                                  uint               aInterpolationOrder,
                                  Vertex_Enrichment* aVertexInterp);
    //------------------------------------------------------------------------------
    Matrix< DDRMat >                  get_coords() const;
    moris_id                          get_id() const;
    moris_index                       get_index() const;
    moris_index                       get_owner() const;
    mtk::Vertex_Interpolation *       get_interpolation( const uint aOrder );
    const mtk::Vertex_Interpolation * get_interpolation( const uint aOrder ) const;
    //------------------------------------------------------------------------------

    Vertex_Enrichment *       get_xtk_interpolation( const uint aOrder );
    Vertex_Enrichment const * get_xtk_interpolation( const uint aOrder ) const ;
    mtk::Vertex const *       get_base_vertex(  ) const;

private:
    // the parent vertex
    mtk::Vertex* mBaseInterpVertex;

    moris_id    mVertexId;
    moris_index mVertexIndex;
    moris_index mVertexOwner;

    // interpolation of vertex
    uint mOrder; /* Order of the Vertex interpolation*/
    Vertex_Enrichment* mInterpolation;

};

std::ostream &
operator<<(std::ostream & os, const xtk::Interpolation_Vertex_Unzipped & dt)
{
    os<<"Vertex Id: "<< dt.get_id()
      << " | Vertex Index: "<<dt.get_index()
      << " | Base Vertex Id: "<<dt.get_base_vertex()->get_id()
      << " | Vertex Interpolation: "<<*(dt.get_xtk_interpolation(0));
    return os;
}
}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_VERTEX_UNZIPPED_HPP_ */
