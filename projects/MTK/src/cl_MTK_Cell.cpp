/*
 * cl_MTK_Cell.cpp
 *
 *  Created on: Aug 17, 2020
 *      Author: kedo3694
 */
#include "cl_MTK_Cell.hpp"

namespace moris
{

    namespace mtk
    {


        //------------------------------------------------------------------------------

        uint
        Cell::get_level() const
        {
            return 0;
        }

        //------------------------------------------------------------------------------

        uint
        Cell::get_number_of_vertices() const
        {
            return this->get_vertex_pointers().size();
        }

        //------------------------------------------------------------------------------

        Matrix< IdMat >
        Cell::get_vertex_ids() const
        {
            uint tNumVertices = this->get_number_of_vertices();
            moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

            Matrix< IdMat > tVertexIds(1, tNumVertices);
            for(uint i = 0; i<tNumVertices; i++)
            {
                tVertexIds(i) = tVertices(i)->get_id();
            }
            return tVertexIds;
        }

        //------------------------------------------------------------------------------

        Matrix< IndexMat >
        Cell::get_vertex_inds() const
        {
            uint tNumVertices = this->get_number_of_vertices();
            moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

            Matrix< IdMat > tVertexInds(1, tNumVertices);
            for(uint i = 0; i<tNumVertices; i++)
            {
                tVertexInds(i) = tVertices(i)->get_index();
            }
            return tVertexInds;
        }

        //------------------------------------------------------------------------------

        Matrix< IndexMat >
        Cell::get_vertex_owners() const
        {
            uint tNumVertices = this->get_number_of_vertices();
            moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

            Matrix< IdMat > tVertexOwners(1, tNumVertices);
            for(uint i = 0; i<tNumVertices; i++)
            {
                tVertexOwners(i) = tVertices(i)->get_owner();
            }
            return tVertexOwners;
        }

        //------------------------------------------------------------------------------

        moris::Cell<mtk::Vertex_Interpolation*>
        Cell::get_vertex_interpolations( const uint aOrder ) const
        {
            uint tNumVerts = this->get_number_of_vertices();
            moris::Cell< mtk::Vertex* > tVertexPointers = this->get_vertex_pointers();
            moris::Cell<mtk::Vertex_Interpolation*> tVertexInterp(tNumVerts);

            for(moris::uint i = 0; i < tNumVerts; i++)
            {
                tVertexInterp(i) =tVertexPointers(i)->get_interpolation(aOrder);
            }

            return tVertexInterp;
        }

        //------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        Cell::get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
        {
            MORIS_ERROR(0,"get_vertices_on_side_ordinal has no default implementation");
            return  moris::Cell<moris::mtk::Vertex const *>(0);
        }

        //------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        Cell::get_geometric_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
        {
            MORIS_ERROR(0,"get_geometric_vertices_on_side_ordinal has no default implementation");
            return  moris::Cell<moris::mtk::Vertex const *>(0);
        }

        //------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Cell::get_cell_physical_coords_on_side_ordinal(moris::moris_index aSideOrdinal) const
        {

            // FIXME: Add assert to check side ordinal

            // get the vertex pointers on the side
            moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = this->get_vertices_on_side_ordinal(aSideOrdinal);

            // allocate output coords (note we do not know the spatial dimension at this time)
            moris::Matrix<moris::DDRMat> tVertexPhysCoords(0,0);

            // iterate through vertices and collect local coordinates
            for(moris::uint i = 0; i < tVerticesOnSide.size(); i++)
            {
                moris::Matrix<moris::DDRMat> tVertexCoord = tVerticesOnSide(i)->get_coords();

                if( i == 0 )
                {
                    MORIS_ASSERT(isrow(tVertexCoord),"Default implementation assumes row based coordinates");
                    tVertexPhysCoords.resize(tVerticesOnSide.size(), tVertexCoord.numel());
                }

                tVertexPhysCoords.get_row(i) = tVertexCoord.get_row(0);
            }

            return tVertexPhysCoords;
        }

        //------------------------------------------------------------------------------

        moris::Matrix< IndexMat >
        Cell::get_vertices_ind_on_side_ordinal(moris::moris_index aSideOrdinal) const
        {
            moris::Cell<moris::mtk::Vertex const *> tVertices = this->get_vertices_on_side_ordinal(aSideOrdinal);

            uint tNumVertices = tVertices.size();

            Matrix< IndexMat > tVertexInd( 1, tNumVertices );

            for(uint i = 0; i < tNumVertices; i++ )
            {
                tVertexInd( 0, i ) = tVertices( i )->get_index();
            }
            return  tVertexInd;
        }

        //------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Cell::compute_outward_side_normal(moris::moris_index aSideOrdinal) const
        {
            MORIS_ERROR(0,"compute_outward_side_normal has no default implementation");
            return  moris::Matrix<moris::DDRMat>(0,0);
        }

        //------------------------------------------------------------------------------

        moris::real
        Cell::compute_cell_measure() const
        {
            MORIS_ERROR(0,"Compute cell measure not implemented");
            return 0;
        }

        //------------------------------------------------------------------------------

        moris::real
        Cell::compute_cell_side_measure(moris_index const & aCellSideOrd) const
        {
            MORIS_ERROR(0,"Compute cell side measure not implemented");
            return 0;
        }

        //------------------------------------------------------------------------------

    }
}
