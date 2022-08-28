/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Proxy.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Logger.hpp"
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Vertex_Proxy.hpp" //MTK/src
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
//------------------------------------------------------------------------------
namespace moris
{
namespace mtk
{

/**
 * \brief the mtk::Cell class provides the cell information that is
 * provided by the mesh.
 */

class Cell_Proxy: public moris::mtk::Cell
{
    Cell_Info*           mCellInfo;
    moris_id             mCellId;
    moris_index          mCellInd;
    moris::Cell<Vertex*> mCellVertices;
    Mesh*                mSTKMeshData;

    //------------------------------------------------------------------------------
public:
    //------------------------------------------------------------------------------

    /**
     * trivial constructor
     */
    Cell_Proxy():mCellInfo(nullptr){};

    //------------------------------------------------------------------------------

    /**
     *  constructor
     */
    Cell_Proxy(       moris_id               aCellId,
                const moris::Cell<Vertex*> & aCellVertices ) : mCellId(aCellId),
                                                               mCellVertices(aCellVertices)
    { };
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_Proxy()
    {};

    //------------------------------------------------------------------------------

    /**
     * returns the domain wide id of the cell
     *
     * @return luint ID
     */
    moris_id
    get_id() const
    {
        return mCellId;
    };

    //------------------------------------------------------------------------------

    /**
     * returns the domain wide id of the cell
     *
     * @return luint ID
     */
    moris_index
    get_index() const
    {
    	MORIS_ERROR(false,"");
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
    };

    //------------------------------------------------------------------------------

    /**
     * returns the proc id of the owner of this cell
     * ( this information is needed for STK )
     */
    moris_id
    get_owner() const
    {
    	MORIS_ERROR(false,"");
        return 0;
    }

    //------------------------------------------------------------------------------

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    moris::Cell< Vertex* >
    get_vertex_pointers() const
    {
        return mCellVertices;
    }

    //------------------------------------------------------------------------------
    /**
     * returns a Matrix with IDs of connected vertices
     */
    Matrix< IdMat >
    get_vertex_ids() const
    {
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< IdMat > tVertexIds(1, tNumVertices);
        for(size_t i = 0; i<tNumVertices; i++)
        {
            tVertexIds(i) = mCellVertices(i)->get_id();
        }
        return tVertexIds;
    }

    /**
     * returns a Mat with indices of connected vertices
     */
    Matrix< IndexMat >
    get_vertex_inds() const
    {
    	MORIS_ERROR(false,"");
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< IndexMat > tVertexInds(1, tNumVertices);
        for(size_t i = 0; i<tNumVertices; i++)
        {
            tVertexInds(i) = mCellVertices(i)->get_index();
        }
        return tVertexInds;
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Matrix of dimension
     * < number of vertices * number of dimensions >
     */
    Matrix< DDRMat >
    get_vertex_coords() const
    {
    	MORIS_ERROR(false,"");
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< DDRMat > tVertexCoords(tNumVertices, mSTKMeshData->get_spatial_dim());
        for(size_t i = 0; i<tNumVertices; i++)
        {
            tVertexCoords.set_row(i,mCellVertices(i)->get_coords());
        }
        return tVertexCoords;
    }

    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {
    	MORIS_ERROR(false,"");
        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

        moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mCellInfo->get_node_to_facet_map(aSideOrdinal);

        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
        for(moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
        {
            tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
        }

        return tVerticesOnSide;
    }

    moris::Matrix<moris::DDRMat>
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const
    {
    	MORIS_ERROR(false,"");
        // get the vertex coordinates
        moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

        // Get vector along these edges
        moris::Matrix<moris::DDRMat> tEdge0Vector(tVertexCoords.numel(),1);
        moris::Matrix<moris::DDRMat> tEdge1Vector(tVertexCoords.numel(),1);

        // Get the nodes which need to be used to compute normal
        moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = mCellInfo->get_node_map_outward_normal(aSideOrdinal);

        // Get vector along these edges
        tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
        tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));

        // Take the cross product to get the normal
        Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector,tEdge1Vector);

        // Normalize
        Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);

        return tUnitOutwardNormal;

    }

    //------------------------------------------------------------------------------

    /**
     * returns an enum that defines the geometry type of the element
     */
    Geometry_Type
    get_geometry_type() const
    {
    	MORIS_ERROR(false,"");
        return mCellInfo->get_cell_geometry();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the order of the element
     */
    Interpolation_Order
    get_interpolation_order() const
    {
    	MORIS_ERROR(false,"");
        return mCellInfo->get_cell_interpolation_order();
    }

};

//------------------------------------------------------------------------------
} /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_ */

