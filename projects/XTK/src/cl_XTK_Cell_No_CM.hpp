/*
 * cl_MTK_Cell_XTK_Impl.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_NO_CM_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_NO_CM_HPP_

#include "cl_MTK_Cell.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Cell_Info.hpp"


using namespace moris;


//------------------------------------------------------------------------------
namespace xtk
{
//------------------------------------------------------------------------------
/**
 * \brief This is the XTK cell implementation when there is a child mesh to use.
 */

class Cell_XTK_No_CM: public moris::mtk::Cell
{
    //------------------------------------------------------------------------------
public:
    //------------------------------------------------------------------------------

    /**
     * trivial constructor
     */
    Cell_XTK_No_CM(){};

    Cell_XTK_No_CM(moris::moris_id                 aElementId,
                   moris::moris_index              aElementIndex,
                   moris::moris_index              aElementOwner,
                   std::shared_ptr<mtk::Cell_Info> aCellInfo,
                   moris::Cell< mtk::Vertex* >     aVertices);
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_XTK_No_CM(){};

    //------------------------------------------------------------------------------
    /**
     * tells how many vertices are connected to this cell
     */
    uint
    get_number_of_vertices() const
    {
        return mCellVertices.size();
    }

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    //FIXME: SDF's Triangle_Vertex causes this to not be able to return a reference.
    moris::Cell< mtk::Vertex* >
    get_vertex_pointers() const
    {
        return mCellVertices;
    }

    void
    set_vertex_pointers(moris::Cell< mtk::Vertex* > & aVertexPointers)
    {
        mCellVertices = aVertexPointers;
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Mat of dimension
     * < number of vertices * number of dimensions >
     */
    Matrix< DDRMat >
    get_vertex_coords() const;

private:
    mtk::Cell_Info const *    mCellInfo = nullptr;
    moris::Cell<mtk::Vertex*> mCellVertices;

    //------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} /* namespace mtk */
//------------------------------------------------------------------------------



#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */
