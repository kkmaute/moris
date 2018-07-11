
#include "cl_Database.hpp" // MTK/src
#include "cl_Hierarchical_Mesh_Main.hpp" // STK/src/Heirarchical
namespace moris
{

// -----------------------------------------------------------------------------------

    /**
     * This function takes adds an additional SDF to the Hierarchical Mesh.
     * There must be at least 3 SDFs in existence
     *
     * @param[in]   Hierarchical Mesh object that is to be processed
     * @param[in]   A Lagrange or B-Spline Mesh
     */
    void
    create_holes(
            Hierarchical_Mesh_Main & aHMR,
            const moris::database  & aMesh );

// -----------------------------------------------------------------------------------
}
