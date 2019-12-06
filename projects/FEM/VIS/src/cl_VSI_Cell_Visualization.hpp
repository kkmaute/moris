/*
 * cl_VIS_Cell_Visualization.hpp
 *
 *  Created on: Dec 02, 2019
 *      Author: schmidt
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Logger.hpp"
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MTK_Cell.hpp" //MRS/CON/src
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"
//------------------------------------------------------------------------------
namespace moris
{
namespace vis
{


/**
 * \brief the mtk::Cell class provides the cell information that is
 * provided by the mesh.
 */

class Cell_Visualization: public moris::mtk::Cell
{
//    Cell_Info*           mCellInfo;
    moris_id             mCellId;
    moris_index          mCellInd;
    moris::Cell< mtk::Vertex* > mCellVertices;
    const moris::mtk::Cell *   mIntegrationcell = nullptr;

    //------------------------------------------------------------------------------
public:
    //------------------------------------------------------------------------------

    /**
     * trivial constructor
     */
    Cell_Visualization() //: mCellInfo(nullptr)
{};

    //------------------------------------------------------------------------------

    /**
     *  constructor
     */
    Cell_Visualization(       moris_id             aCellId,
                              moris_index          aCellInd,
                        const moris::Cell< mtk::Vertex* > aCellVertices,
                        const moris::mtk::Cell *   aIntegrationCell ) : mCellId( aCellId ),
                                                                        mCellInd( aCellInd ),
                                                                        mCellVertices( aCellVertices ),
                                                                        mIntegrationcell( aIntegrationCell )
    { };
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_Visualization()
    {
    };

    //------------------------------------------------------------------------------
    /**
     * returns the domain wide id of the cell
     *
     * @return luint ID
     */
    moris_id get_id() const
    {
        return mCellId;
    };

    //------------------------------------------------------------------------------
    /**
     * returns the domain wide id of the cell
     *
     * @return luint ID
     */
    moris_index get_index() const
    {
        return mCellInd;
    }

    //------------------------------------------------------------------------------
    /**
     * tells how many vertices are connected to this cell
     */
    uint get_number_of_vertices() const
    {
        return mCellVertices.size();
    };

    //------------------------------------------------------------------------------

    /**
     * returns the proc id of the owner of this cell
     * ( this information is needed for STK )
     */
    moris_id get_owner() const
    {
        return mIntegrationcell->get_owner();
    }

    //------------------------------------------------------------------------------
    /**
     * fills a moris::cell with pointers to connected vertices
     */
    moris::Cell< mtk::Vertex* > get_vertex_pointers() const
    {
        return mCellVertices;
    }

    //------------------------------------------------------------------------------
    /**
     * returns a Matrix with IDs of connected vertices
     */
    Matrix< IdMat > get_vertex_ids() const
    {
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< IdMat > tVertexIds( 1, tNumVertices );
        for(size_t i = 0; i<tNumVertices; i++)
        {
            tVertexIds(i) = mCellVertices(i)->get_id();
        }
        return tVertexIds;
    }

    //------------------------------------------------------------------------------
    /**
     * returns a Mat with indices of connected vertices
     */
    Matrix< IndexMat > get_vertex_inds() const
    {
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< IndexMat > tVertexInds( 1, tNumVertices );
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
    Matrix< DDRMat > get_vertex_coords() const
    {
        size_t tNumVertices = this->get_number_of_vertices();
        Matrix< DDRMat > tVertexCoords( tNumVertices, mCellVertices( 0 )->get_coords().numel() );
        for(size_t i = 0; i<tNumVertices; i++)
        {
            tVertexCoords.set_row(i,mCellVertices(i)->get_coords() );
        }
        return tVertexCoords;
    }

    //------------------------------------------------------------------------------
    moris::Cell< moris::mtk::Vertex const * > get_vertices_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
//        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();
//
//        moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mCellInfo->get_node_to_facet_map(aSideOrdinal);
//
//        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
//        for(moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
//        {
//            tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
//        }
//
//        return tVerticesOnSide;
        return moris::Cell< moris::mtk::Vertex const * >( 0 );
    }

    //------------------------------------------------------------------------------
    moris::Matrix< moris::DDRMat > compute_outward_side_normal( moris::moris_index aSideOrdinal ) const
    {
//        // get the vertex coordinates
//        moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();
//
//        // Get vector along these edges
//        moris::Matrix<moris::DDRMat> tEdge0Vector(tVertexCoords.numel(),1);
//        moris::Matrix<moris::DDRMat> tEdge1Vector(tVertexCoords.numel(),1);
//
//        // Get the nodes which need to be used to compute normal
//        moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = mCellInfo->get_node_map_outward_normal(aSideOrdinal);
//
//        // Get vector along these edges
//        tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
//        tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));
//
//        // Take the cross product to get the normal
//        Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector,tEdge1Vector);
//
//        // Normalize
//        Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);
//
//        return tUnitOutwardNormal;
        return moris::Matrix<moris::DDRMat>( 0, 0 );
    }

    //------------------------------------------------------------------------------
    /**
     * returns an enum that defines the geometry type of the element
     */
    mtk::Geometry_Type get_geometry_type() const
    {
        MORIS_ERROR( false, "get_geometry_type(), not implemented for visualization mesh");
        return mtk::Geometry_Type::UNDEFINED;
    }

    //------------------------------------------------------------------------------
    /**
     * returns the order of the element
     */
    mtk::Interpolation_Order get_interpolation_order() const
    {
        MORIS_ERROR( false, "get_interpolation_order(), not implemented for visualization mesh");
        return mtk:: Interpolation_Order::UNDEFINED;
    }

};

//------------------------------------------------------------------------------
} /* namespace vis */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_VIS_CELL_VISULIZATION_HPP_ */
