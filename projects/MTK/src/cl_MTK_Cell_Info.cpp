#include "cl_MTK_Cell_Info.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"
#include "cl_MTK_Enums.hpp"
namespace moris
{
    namespace mtk
    {
        // ---------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info::get_geometric_node_to_facet_map() const
        {
            MORIS_ERROR(0,"ERROR");
            return moris::Matrix<moris::IndexMat>();
        }

        // ---------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info::get_geometric_node_to_facet_map(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"ERROR");
            return moris::Matrix<moris::IndexMat>();
        }

        // ---------------------------------------------------------------------------------

        moris::uint
        Cell_Info::get_adjacent_side_ordinal(moris::uint aSideOrdinal) const
        {
            MORIS_ERROR(0,"get_adjacent_side_ordinal only makes sense for hex/quad type cells");
            return MORIS_UINT_MAX;
        }

        // ---------------------------------------------------------------------------------
        Matrix<DDRMat>
        Cell_Info::get_vertex_loc_coord(moris_index aVertexOrdinal) const
        {
            MORIS_ERROR(0,"get_loc_coord_on_side_ordinal not implemented for given cell info type");
            return Matrix<DDRMat>(0,0);
        }

        // ---------------------------------------------------------------------------------

        void
        Cell_Info::get_loc_coord_on_side_ordinal(
                moris::uint const & aSideOrdinal,
                Matrix<DDRMat>    & aXi ) const
        {
            // number of vertices on facet side
            uint tNumVertsPerFacet = this->get_num_verts_per_facet();

            // dimension of local coordinates
            uint tDimXi = this->get_loc_coord_dim();

            // allocate space
            aXi.resize(tNumVertsPerFacet,tDimXi);

            // get the vertex ordinals on facet
            Matrix<IndexMat> tVerticesOnSide = this->get_node_to_facet_map(aSideOrdinal);

            MORIS_ASSERT(tNumVertsPerFacet == tVerticesOnSide.numel(),"Dimension mismatch between get_num_verts_per_facet() and get_node_to_facet_map()");

            // iterate through vertices and get their local coordinate
            for(moris::uint i = 0; i < tNumVertsPerFacet; i++)
            {
                aXi.get_row(i) = this->get_vertex_loc_coord(tVerticesOnSide(i)).get_row(0);
            }

        }

        // ---------------------------------------------------------------------------------

        void
        Cell_Info::get_loc_coords_of_cell(Matrix<DDRMat> & aXi) const
        {
            // number of vertices
            uint tNumVerts = this->get_num_verts();

            // dimension of local coordinates
            uint tDimXi = this->get_loc_coord_dim();

            // allocate space
            aXi.resize(tNumVerts,tDimXi);

            // iterate through vertices and get their local coordinate
            for(moris::uint i = 0; i < tNumVerts; i++)
            {
                aXi.get_row(i) = this->get_vertex_loc_coord((moris_index) i).get_row(0);
            }
        }

        // ---------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Cell_Info::get_edge_to_face_map() const
        {
            MORIS_ERROR(0,"get edge to face map not implemented for this CELL_INFO, this CELL_INFO currently used in XTK only for tets/tris");
            return moris::Matrix<moris::IndexMat>(0,0);
        }

        // ---------------------------------------------------------------------------------

        void
        Cell_Info::eval_N( const Matrix< DDRMat > & aXi,
                Matrix< DDRMat > & aNXi ) const
        {
            MORIS_ERROR(0,"eval_N not implemented for this type of cell info");
        }


    }

}

