/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_MTK_Load_External_Surface_Mesh.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_assert.hpp"

namespace moris::mtk
{
    //-------------------------------------------------------------------------------
    // These functions provide tools to read surface meshes from various file formats
    // and convert them to the MORIS surface mesh format. The MORIS surface mesh format
    // is a Matrix of dimension <spatial_dimension> x <number_of_vertices> where each
    // column is a vertex in the mesh. The mesh is represented by a list of facets,
    // where each facet is a list of vertex indices. The vertex indices are the column
    // indices of the vertex matrix. The facets are stored in a Vector< Vector< moris_index > >.
    // Surface meshes are assumed to be water tight without any hanging vertices or edges.
    // They can however contain holes formed by a closed loop of facets.


    // Example: 2D triangular surface mesh with 3 vertices and 3 facets would have the following format:

    // Matrix< DDRMat > tCoords = { { 0, 1.25, 1.5 }, { 0, 0, 1.0 } };
    // Vector< Vector< moris_index > > tFacets = { { 0, 1 }, { 1, 2 }, { 2, 0 } };

    //-------------------------------------------------------------------------------

    /**
     * loads an ASCII file into a buffer of strings.
     */
    void
    load_ascii_to_buffer( const std::string& aFilePath,
            Vector< std::string >&           aBuffer );

    //-------------------------------------------------------------------------------

    /**
     * loads an ascii file and creates vertex and facet objects
     * Facets are either lines in 2D or triangles in 3D
     */
    Matrix< DDRMat >
    load_vertices_from_object_file( const std::string& aFilePath, const Vector< real >& aOffsets, const Vector< real >& aScale );

    //-------------------------------------------------------------------------------

    /**
     * loads an ascii file and creates vertex and facet objects
     * Facets are either lines in 2D or triangles in 3D
     */
    Vector< Vector< moris_index > >
    load_facets_from_object_file( const std::string& aFilePath );


}    // namespace moris::mtk
