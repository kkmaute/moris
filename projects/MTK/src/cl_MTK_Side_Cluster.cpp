
#include "cl_MTK_Side_Cluster.hpp"
#include "typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell.hpp"

#include "cl_MTK_Cluster.hpp"


namespace moris
{
namespace mtk
{

    // ----------------------------------------------------------------------------------

    Side_Cluster::Side_Cluster(){}

    // ----------------------------------------------------------------------------------

    moris::Cell<mtk::Cell const *> const &
    Side_Cluster::get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster ) const
    {
        return this->get_cells_in_side_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_measure(
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Master_Slave aIsMaster     ) const
    {
        moris::real tVolume = 0.0;
        moris::Cell<moris::mtk::Cell const *> const* tCells = nullptr;
        if(aPrimaryOrVoid == mtk::Primary_Void::PRIMARY)
        {
            tCells = &this->get_primary_cells_in_cluster();
        }
        else
        {
            tCells = & this->get_void_cells_in_cluster();
        }

        for(auto iC = tCells->cbegin(); iC < tCells->cend(); iC++)
        {
            tVolume = tVolume+(*iC)->compute_cell_measure();
        }

        return tVolume;
    }

    // ----------------------------------------------------------------------------------

    moris::real
    Side_Cluster::compute_cluster_cell_side_measure(
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Master_Slave aIsMaster     ) const
    {
        MORIS_ASSERT(aPrimaryOrVoid == mtk::Primary_Void::PRIMARY,"Side cluster only operates on primary cells.");

        moris::real tMeasure = 0.0;
        moris::Cell<mtk::Cell const *> const & tCells = this->get_primary_cells_in_cluster();
        moris::Matrix<IndexMat> tSideOrds = this->get_cell_side_ordinals(aIsMaster);

        for(moris::uint iC = 0 ; iC < tCells.size(); iC++)
        {
            tMeasure = tMeasure+tCells(iC)->compute_cell_side_measure(tSideOrds(iC));
        }

        return tMeasure;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Side_Cluster::get_cell_indices_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        moris::Cell<moris::mtk::Cell const *> const & tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix<moris::IndexMat> tCellIndices(1,tNumCells);

        // get cell indices and store
        for(moris::uint i = 0 ; i < tNumCells; i++)
        {
            tCellIndices(i) = tCells(i)->get_index();
        }

        return tCellIndices;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Side_Cluster::get_primary_cell_indices_in_cluster() const
    {
        return this->get_cell_indices_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index
    Side_Cluster::get_interpolation_cell_index() const
    {
        return this->get_interpolation_cell().get_index();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Side_Cluster::get_vertex_indices_in_cluster() const
    {
        // number of cells in cluster
         moris::uint tNumVertices = this->get_num_vertices_in_cluster();

         // cell access
         moris::Cell<moris::mtk::Vertex const *> tVertices = this->get_vertices_in_cluster();

         // initialize output
         moris::Matrix<moris::IndexMat> tVertexIndices(1,tNumVertices);

         // get cell indices and store
         for(moris::uint i = 0 ; i < tNumVertices; i++)
         {
             tVertexIndices(i) = tVertices(i)->get_index();
         }

         return tVertexIndices;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IdMat>
    Side_Cluster::get_cell_ids_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        moris::Cell<moris::mtk::Cell const *> const & tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix<moris::IdMat> tCellIds(1,tNumCells);

        // get cell indices and store
        for(moris::uint i = 0 ; i < tNumCells; i++)
        {
            tCellIds(i) = tCells(i)->get_id();
        }

        return tCellIds;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IdMat>
    Side_Cluster::get_primary_cell_ids_in_cluster() const
    {
        return this->get_cell_ids_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::IndexMat>
    Side_Cluster::get_vertex_ids_in_cluster() const
    {
        // number of cells in cluster
         moris::uint tNumVertices = this->get_num_vertices_in_cluster();

         // cell access
         moris::Cell<moris::mtk::Vertex const *> tVertices = this->get_vertices_in_cluster();

         // initialize output
         moris::Matrix<moris::IndexMat> tVertexIds(1,tNumVertices);

         // get cell indices and store
         for(moris::uint i = 0 ; i < tNumVertices; i++)
         {
             tVertexIds(i) = tVertices(i)->get_id();
         }

         return tVertexIds;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    Side_Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index      aClusterLocalIndex,
            const mtk::Master_Slave aIsMaster ) const
    {
        MORIS_ASSERT(aClusterLocalIndex < (moris_index)this->get_num_sides_in_cluster(),"Integration Cell Cluster index out of bounds");

        // get side ordinal of interest
        moris_index tSideOrdinal = this->get_cell_side_ordinal(aClusterLocalIndex);

        // get the integration cell of interest
        moris::mtk::Cell const * tIntegrationCell = this->get_cells_in_side_cluster()(aClusterLocalIndex);

        // get the vertex pointers on the side
        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tIntegrationCell->get_vertices_on_side_ordinal(tSideOrdinal);

        // allocate output (nnode x dim_xsi)
        moris::Matrix<moris::DDRMat> tVertexParamCoords( tVerticesOnSide.size(), this->get_dim_of_param_coord());

        // iterate through vertices and collect local coordinates
        for(moris::uint i = 0; i < tVerticesOnSide.size(); i++)
        {
            tVertexParamCoords.get_row(i) = this->get_vertex_local_coordinate_wrt_interp_cell(tVerticesOnSide(i)).get_row(0);
        }

        return tVertexParamCoords;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    Side_Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(
            moris::moris_index      aPrimaryCellClusterIndex,
            const mtk::Master_Slave aIsMaster ) const
    {
        return this ->get_cell_local_coords_on_side_wrt_interp_cell( aPrimaryCellClusterIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_sides_in_cluster() const
    {
        return this->get_cells_in_side_cluster().size();
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_primary_cells() const
    {
        return this->get_num_sides_in_cluster();
    }

    // ----------------------------------------------------------------------------------

    moris::uint
    Side_Cluster::get_num_vertices_in_cluster() const
    {
        return this->get_vertices_in_cluster().size();
    }

    // ----------------------------------------------------------------------------------







}
}
