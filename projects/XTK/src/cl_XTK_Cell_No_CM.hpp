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

namespace xtk
{
class Background_Mesh;
}

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

    Cell_XTK_No_CM(moris::moris_id             aElementId,
                   moris::moris_index          aElementIndex,
                   moris::moris_index          aElementOwner,
                   mtk::Cell_Info const *      aCellInfo,
                   moris::Cell< mtk::Vertex* > aVertices);
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_XTK_No_CM(){};

    //------------------------------------------------------------------------------

    /**
     * returns the domain wide id of the cell
     *
     * @return moris_id ID
     */
    moris_id
    get_id() const
    {
        return mCellId;
    }

    //------------------------------------------------------------------------------

    /**
     * returns the local index of the cell
     *
     * @return moris_index ID
     */
    moris_index
    get_index() const
    {
        return mCellInd;
    }

    //------------------------------------------------------------------------------

    /**
     * tells how many vertices are connected to this cell
     */
    uint
    get_number_of_vertices() const
    {
        return mCellVertices.size();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the proc id of the owner of this cell
     * ( this information is needed for STK )
     */
    moris_id
    get_owner() const
    {
        return mCellOwner;
    }

    //------------------------------------------------------------------------------

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    //FIXME: SDF's Triangle_Vertex causes this to not be able to return a reference.
    moris::Cell< mtk::Vertex* >
    get_vertex_pointers() const
    {
        return mCellVertices;
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Mat of dimension
     * < number of vertices * number of dimensions >
     */
    Matrix< DDRMat >
    get_vertex_coords() const;

    //------------------------------------------------------------------------------

    /**
     * returns an enum that defines the geometry type of the element
     */
    mtk::Geometry_Type
    get_geometry_type() const
    {
        return mCellInfo->get_cell_geometry();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the order of the element
     */
    mtk::Interpolation_Order
    get_interpolation_order() const
    {
        return mCellInfo->get_cell_interpolation_order();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the integration order of the element
     */
    mtk::Integration_Order
    get_integration_order() const
    {
        return mCellInfo->get_cell_integration_order();
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const;


    //------------------------------------------------------------------------------
    moris::real
    compute_cell_measure() const;

    //------------------------------------------------------------------------------
    moris::real
    compute_cell_measure_deriv(uint aLocalVertexID, uint aDirection) const;

    moris::real
    compute_cell_side_measure(moris_index const & aSideOrdinal) const;

private:
    mtk::Cell_Info const *    mCellInfo;
    moris_id                  mCellId;
    moris_index               mCellInd;
    moris_index               mCellOwner;
    moris::Cell<mtk::Vertex*> mCellVertices;

    //------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} /* namespace mtk */
//------------------------------------------------------------------------------



#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */