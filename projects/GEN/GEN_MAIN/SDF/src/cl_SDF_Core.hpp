/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.hpp
 *
 */

#pragma once

#include "typedefs.hpp"    // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Object.hpp"
#include "cl_SDF_Raycast.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Parameters.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    class Core
    {
        Mesh& mMesh;
        Object& mObject;

        uint mSurfaceElements = 0;
        uint mVolumeElements = 0;
        real mBufferDiagonal;

        uint mCandidateSearchDepth        = 1;
        real mCandidateSearchDepthEpsilon = 0.01;
        bool mVerbose;

        //-------------------------------------------------------------------------------

      public:
        //-------------------------------------------------------------------------------

        Core( Mesh&   aMesh,
                Object& aObject,
                bool  aVerbose = false );

        //-------------------------------------------------------------------------------

        ~Core(){};

        //-------------------------------------------------------------------------------

        void
        set_candidate_search_depth( const uint aCandidateSearchDepth )
        {
            mCandidateSearchDepth = aCandidateSearchDepth;
        }

        //-------------------------------------------------------------------------------

        void
        set_candidate_search_epsilon( const real aCandidateSearchEpsilon )
        {
            mCandidateSearchDepthEpsilon = aCandidateSearchEpsilon;
        }

        //-------------------------------------------------------------------------------

        void
        raycast_mesh();

        //-------------------------------------------------------------------------------

        void
        calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface );

        //-------------------------------------------------------------------------------
        void
        calculate_raycast(
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume );

        //-------------------------------------------------------------------------------

        void
        calculate_raycast_and_sdf( Matrix< DDRMat >& aSDF );

        //-------------------------------------------------------------------------------

        void
        calculate_raycast_and_sdf(
                Matrix< DDRMat >&   aSDF,
                Matrix< IndexMat >& aElementsAtSurface,
                Matrix< IndexMat >& aElementsInVolume );

        //-------------------------------------------------------------------------------

        void
        save_to_vtk( const std::string& aFilePath );

        //-------------------------------------------------------------------------------

        void
        save_unsure_to_vtk( const std::string& aFilePath );

        //-------------------------------------------------------------------------------

      private:
        //-------------------------------------------------------------------------------

        // void
        // voxelize_3D( const uint aAxis );

        //-------------------------------------------------------------------------------

        moris::Cell< Vertex* >
        set_candidate_list();

        //-------------------------------------------------------------------------------

        void
        calculate_udf( moris::Cell< Vertex* >& aCandidateList );

        //-------------------------------------------------------------------------------

        /**
         * Kehrwoche :
         * make sure that each vertex is really associated
         * to its closest triangle
         */
        void
        sweep();

        //-------------------------------------------------------------------------------

        void
        fill_sdf_with_values( Matrix< DDRMat >& aSDF );
                
        //-------------------------------------------------------------------------------

        void
        calculate_candidate_points_and_buffer_diagonal();

        //-------------------------------------------------------------------------------

        void
        get_nodes_withing_bounding_box_of_triangle(
                Facet*                  aFacet,
                moris::Cell< Vertex* >& aNodes,
                moris::Cell< Vertex* >& aCandList );

        //--------------------------------------------------------------------------------
    };

    //-------------------------------------------------------------------------------
} /* namespace moris::sdf */
