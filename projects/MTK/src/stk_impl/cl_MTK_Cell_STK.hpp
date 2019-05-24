/*
 * cl_MTK_Cell_STK.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Logger.hpp"
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MTK_Vertex_STK.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Tetra4_Connectivity.hpp"
#include "cl_MTK_Hex8_Connectivity.hpp"
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

class Cell_STK: public moris::mtk::Cell
{
    enum CellTopology    mCellType;
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
    Cell_STK(){};

    //------------------------------------------------------------------------------

    /**
     *  constructor
     */
    Cell_STK( enum CellTopology               aCellType,
                   moris_id               aCellId,
                   moris_index            aCellInd,
                   const moris::Cell<Vertex*> & aCellVertices,
                   Mesh* aStkImplementation):
                       mCellId(aCellId),
                       mCellInd(aCellInd),
                       mCellVertices(aCellVertices),
                       mSTKMeshData(aStkImplementation)
    {

    };
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_STK()
    {
        //                delete mInterpolation;
    };

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
        return mSTKMeshData->get_entity_owner( mCellInd, EntityRank::ELEMENT);
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
        Geometry_Type tGeomType = this->get_geometry_type();

        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

        switch(tGeomType)
        {
            case(Geometry_Type::TET):
            {
                if(this->get_number_of_vertices() == 4)
                {
                    MORIS_ASSERT(aSideOrdinal<4,"Side ordinal out of bounds for cell type tet");
                    moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = moris::Tetra4_Connectivity::get_node_to_face_map(aSideOrdinal);

                    moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(3);
                    tVerticesOnSide(0) = tVertices(tNodeOrdsOnSide(0));
                    tVerticesOnSide(1) = tVertices(tNodeOrdsOnSide(1));
                    tVerticesOnSide(2) = tVertices(tNodeOrdsOnSide(2));
                    return tVerticesOnSide;
                    break;
                }
                else
                {
                    MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
                    return moris::Cell<moris::mtk::Vertex const *>(0);
                    break;
                }

            }
            case(Geometry_Type::HEX):
            {

                if(this->get_number_of_vertices() == 8)
                {
                    MORIS_ASSERT(aSideOrdinal<6,"Side ordinal out of bounds for cell type hex");
                    moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = moris::Hex8::get_node_to_face_map(aSideOrdinal);

                    moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(4);
                    tVerticesOnSide(0) = tVertices(tNodeOrdsOnSide(0));
                    tVerticesOnSide(1) = tVertices(tNodeOrdsOnSide(1));
                    tVerticesOnSide(2) = tVertices(tNodeOrdsOnSide(2));
                    tVerticesOnSide(3) = tVertices(tNodeOrdsOnSide(3));
                    return tVerticesOnSide;
                    break;
                }
                else
                {
                    MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
                    return moris::Cell<moris::mtk::Vertex const *>(0);
                    break;
                }

            }

            default:
                MORIS_ERROR(0,"Invalid geometry type, currently only supports hex8 and tet4");
                return moris::Cell<moris::mtk::Vertex const *>(0);
                break;
        }
    }

    moris::Matrix<moris::DDRMat>
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const
    {

        enum Geometry_Type tGeomType = this->get_geometry_type();

        // Get vector along these edges
        moris::Matrix<moris::DDRMat> tEdge0Vector(3,1);
        moris::Matrix<moris::DDRMat> tEdge1Vector(3,1);


        switch(tGeomType)
        {
            case(Geometry_Type::HEX):
            {
                MORIS_ERROR(aSideOrdinal<6,"Side ordinal out of bounds.");

#ifdef DEBUG
                if(this->get_vertex_pointers().size() > 8)
                {
                 MORIS_LOG_DEBUG("Warning: this normal computation only valid for flat facets. Ensure your higher order element has flat facets");
                }
#endif

                // get the vertex coordinates
                moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

                // Get the nodes which need to be used to compute normal
                moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = Hex8::get_node_map_outward_normal(aSideOrdinal);

                // Get vector along these edges
                tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
                tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));
                break;
            }
            case(Geometry_Type::TET):
            {
                MORIS_ERROR(aSideOrdinal<4,"Side ordinal out of bounds.");

#ifdef DEBUG
                if(this->get_vertex_pointers().size() > 4)
                {
                 MORIS_LOG_DEBUG("Warning: this normal computation only valid for flat facets. Ensure your higher order element has flat facets");
                }
#endif

                // get the vertex coordinates
                moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

                // Get the nodes which need to be used to compute normal
                moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = Tetra4_Connectivity::get_node_map_outward_normal(aSideOrdinal);

                // Get vector along these edges
                tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
                tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));

                break;
            }

            default:
                MORIS_ERROR(0,"Only implemented for hex8 compute_outward_side_normal");
                break;
        }

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
        //MORIS_ASSERT(false," Cell_STK::get_geometry_type - Not implemented");
        //return Geometry_Type::UNDEFINED;

        //FIXME: only for Lagrange LINE, QUAD, HEX, TET, TRI
        switch ( mSTKMeshData->get_spatial_dim() )
        {
            case ( 1 ) :
                                return Geometry_Type::LINE;
            break;

            case ( 2 ) :
                        {
                uint tNumEdges = mSTKMeshData->get_edges_connected_to_element_loc_inds( this->get_index() ).numel();
                switch ( tNumEdges )
                {
                    case ( 3 ):
                                        return Geometry_Type::TRI;
                    break;

                    case ( 4 ):
                                        return Geometry_Type::QUAD;
                    break;

                    default:
                        return Geometry_Type::UNDEFINED;
                        break;
                }
                break;
                        }

            case ( 3 ) :
                        {
                uint tNumFaces = mSTKMeshData->get_faces_connected_to_element_loc_inds( this->get_index() ).numel();
                switch ( tNumFaces )
                {
                    case ( 4 ):
                                       return Geometry_Type::TET;
                    break;

                    case ( 6 ):
                                       return Geometry_Type::HEX;
                    break;

                    default:
                        return Geometry_Type::UNDEFINED;
                        break;
                }
                break;
                        }

            default :
                return Geometry_Type::UNDEFINED;
                break;
        }
    }

    //------------------------------------------------------------------------------

    /**
     * returns the order of the element
     */
    Interpolation_Order
    get_interpolation_order() const
    {
        MORIS_ASSERT(false,"Cell_STK::get_interpolation_order - Not implemented");
        return Interpolation_Order::UNDEFINED;
    }

};

//------------------------------------------------------------------------------
} /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_ */
