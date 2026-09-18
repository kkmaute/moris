/*
 * fn_MTK_Gather_Surface_Mesh_Arrays.hpp
 *
 * Utility that takes the *local* (per-processor) mtk::Surface_Mesh view of an
 * IG mesh's side set(s) and turns it into three global arrays that end up
 * IDENTICAL on every processor (gather + broadcast, i.e. an "allgather"):
 *
 *   1) aGlobalCells         - cell -> vertex connectivity (global vertex ids,
 *                              renumbered/compacted so shared vertices on
 *                              partition boundaries only appear once)
 *   2) aGlobalVertexIds     - the original (mesh-wide) global vertex index
 *                              for each row of aGlobalVertexCoords
 *   3) aGlobalVertexCoords  - (spatial_dim x n_unique_vertices) coordinates
 *
 * Implementation: gather everything onto one root processor with moris'
 * existing Communication_Tools::gatherv_mats (handles processors
 * contributing different-sized matrices), de-duplicate/assemble there, then
 * broadcast the assembled result back out with Communication_Tools::broadcast_mat
 * so every rank ends up with the same three arrays.
 */

#pragma once

#include "cl_MTK_Integration_Surface_Mesh.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::mtk
{
    /**
     * @brief Gathers a per-processor Integration_Surface_Mesh into global cell/vertex/coordinate
     *        arrays and makes them available, identically, on every processor.
     *
     * @param aLocalSurfaceMesh   the local surface mesh on the calling processor
     *                            (may be empty, i.e. own zero cells/vertices, on some ranks)
     * @param aGlobalCells        [out] on every rank: one entry per global cell,
     *                            each entry holds the (compacted) indices into
     *                            aGlobalVertexIds / aGlobalVertexCoords for that cell's vertices
     * @param aGlobalCellIds      [out] on every rank: (n_cells x 1) original mesh-wide
     *                            global cell index for each entry of aGlobalCells (same order)
     * @param aGlobalCellOwners   [out] on every rank: (n_cells x 1) owning processor rank
     *                            (0-based) for each entry of aGlobalCells (same order) --
     *                            i.e. which processor's local Surface_Mesh this cell came from
     * @param aGlobalVertexIds    [out] on every rank: (n x 1) original mesh-wide
     *                            global vertex index for each compacted vertex
     * @param aGlobalVertexCoords [out] on every rank: (d x n) vertex coordinates,
     * @param aGlobalVertexDisplacements [out] on every rank: (d x n) vertex displacements,
     *                            column i corresponds to aGlobalVertexIds( i )
     * @param aRootProc           rank used as the internal staging processor for assembly
     *                            before the result is broadcast back out to everyone;
     *                            since the output is identical on all ranks it does not
     *                            matter which valid rank you pick, it just needs to be
     *                            the SAME value on every rank (default 0 is fine)
     */
    void
    gather_surface_mesh_arrays(
            Integration_Surface_Mesh const & aLocalSurfaceMesh,
            Vector< Matrix< IndexMat > >&    aGlobalCells,
            Matrix< IndexMat >&              aGlobalCellIds,
            Matrix< IndexMat >&              aGlobalCellOwners,
            Matrix< IndexMat >&              aGlobalVertexIds,
            Matrix< DDRMat >&                aGlobalVertexCoords,
            Matrix< DDRMat >&                aGlobalVertexDisplacements,
            moris_index                      aRootProc = 0 );

}    // namespace moris::mtk