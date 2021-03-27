/*
 * cl_XTK_Interpolation_Cell.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
using namespace moris;

namespace moris
{
    namespace mtk
    {
    class Vertex;
    }
}


namespace xtk
{

class Interpolation_Cell: public mtk::Cell
{
public:
    Interpolation_Cell(moris_id                  aCellId,
                       moris_index               aCellIndex,
                       moris_id                  aCellOwner,
                       moris::mtk::Cell_Info* aConnectivity):
                           mCellId(aCellId),
                           mCellIndex(aCellIndex),
                           mCellOwner(aCellOwner),
                           mCellInfo(aConnectivity)
    {
    }

    Interpolation_Cell(){};

    moris_id
    get_id() const {return mCellId;};

    moris_index
    get_index() const {return mCellIndex;};

    moris_id
    get_owner() const {return mCellOwner;};

    virtual moris::Cell< mtk::Vertex* >
    get_vertex_pointers() const = 0;

    virtual Matrix< DDRMat >
    get_vertex_coords() const = 0;

    mtk::Geometry_Type
    get_geometry_type() const { return mCellInfo->get_cell_geometry(); };


    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const;


    moris::Cell<moris::mtk::Vertex const *>
    get_geometric_vertices_on_side_ordinal( moris::moris_index aSideOrdinal )  const;

    moris::Matrix<moris::DDRMat>
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const;

    mtk::Interpolation_Order
    get_interpolation_order() const{ return mCellInfo->get_cell_interpolation_order();};

    mtk::Integration_Order
    get_integration_order() const{ return mCellInfo->get_cell_integration_order();};

    void set_id(moris_id aId);

    moris::mtk::Cell_Info const *
    get_connectivity() const {return mCellInfo;}

    moris::mtk::Cell_Info *
    get_connectivity() {return mCellInfo;}

    moris::real
    compute_cell_measure() const { return mCellInfo->compute_cell_size(this);}

    moris::real
    compute_cell_measure_deriv(uint aLocalVertexID, uint aDirection) const
    {
        return mCellInfo->compute_cell_size_deriv(this, aLocalVertexID, aDirection);
    }

    moris::real
    compute_cell_side_measure(moris_index const & aCellSideOrd) const
    {
        return mCellInfo->compute_cell_side_size(this,aCellSideOrd);
    }

    moris::real
    compute_cell_side_measure_deriv(moris_index const & aCellSideOrd, uint aLocalVertexID, uint aDirection) const
    {
        return mCellInfo->compute_cell_side_size_deriv(this, aCellSideOrd, aLocalVertexID, aDirection);
    }

private:
    moris_id               mCellId;
    moris_id               mCellIndex;
    moris_id               mCellOwner;
    moris::mtk::Cell_Info* mCellInfo;
};

}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_HPP_ */
