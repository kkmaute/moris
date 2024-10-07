/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Generator.cpp
 *
 */

#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Core.hpp"
#include "cl_Tracer.hpp"

#include "cl_SDF_Generator.hpp"

#include <utility>

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    SDF_Generator::SDF_Generator(
            const std::string& aObjectPath,
            const bool         aVerboseFlag )
            : mObject( aObjectPath )
            , mVerboseFlag( aVerboseFlag )
    {
    }

    //-------------------------------------------------------------------------------

        SDF_Generator::SDF_Generator(
                const std::string& aObjectPath,
                Vector< real >&      aObjectOffset,
                const bool         aVerboseFlag )
                : mObject( aObjectPath, 1e-8, aObjectOffset )
                , mVerboseFlag( aVerboseFlag )
        {
        }

    //-------------------------------------------------------------------------------

    void
    SDF_Generator::raycast(
            mtk::Mesh*          aMesh,
            Matrix< IndexMat >& aElementsAtSurface )
    {
        // create mesh wrapper
        Mesh tMesh( aMesh );

        // create core
        Core tCore( tMesh, mObject );

        // perform raycast
        tCore.calculate_raycast( aElementsAtSurface );
    }
    //-------------------------------------------------------------------------------

    void
    SDF_Generator::raycast(
            const std::shared_ptr< mtk::Mesh >& aMesh,
            Matrix< IndexMat >&                 aElementsAtSurface )
    {
        // create mesh wrapper
        Mesh tMesh( aMesh );

        // create core
        Core tCore( tMesh, mObject );

        // perform raycast
        tCore.calculate_raycast( aElementsAtSurface );
    }

    //-------------------------------------------------------------------------------

    void
    SDF_Generator::raycast(
            mtk::Mesh*          aMesh,
            Matrix< IndexMat >& aElementsAtSurface,
            Matrix< IndexMat >& aElementsInVolume )
    {
        // create mesh wrapper
        Mesh tMesh( aMesh );

        // create core
        Core tCore( tMesh, mObject );

        // perform raycast
        tCore.calculate_raycast( aElementsAtSurface, aElementsInVolume );
    }

    //-------------------------------------------------------------------------------

    void
    SDF_Generator::raycast(
            const std::shared_ptr< mtk::Mesh >& aMesh,
            Matrix< IndexMat >&                 aElementsAtSurface,
            Matrix< IndexMat >&                 aElementsInVolume )
    {
        // create mesh wrapper
        Mesh tMesh( aMesh, mVerboseFlag );

        // create core
        Core tCore( tMesh, mObject, mVerboseFlag );

        // perform raycast
        tCore.calculate_raycast( aElementsAtSurface, aElementsInVolume );
    }

    //-------------------------------------------------------------------------------

    void
    SDF_Generator::calculate_sdf(
            mtk::Mesh*        aMesh,
            Matrix< DDRMat >& aSDF )
    {
        // trace this function
        Tracer tTracer( "GEN", "SDF-Generator", "Compute SDF" );

        // create mesh wrapper
        Mesh tMesh( aMesh, mVerboseFlag );

        // create core
        Core tCore( tMesh, mObject, mVerboseFlag );

        // calculate SDF
        tCore.calculate_raycast_and_sdf( aSDF );

        // tCore.save_to_vtk( "sdf_mesh.vtk");
    }
    //-------------------------------------------------------------------------------

    void
    SDF_Generator::calculate_sdf(
            const std::shared_ptr< mtk::Mesh >& aMesh,
            Matrix< DDRMat >&                   aSDF )
    {
        // create mesh wrapper
        Mesh tMesh( aMesh, mVerboseFlag );

        // create core
        Core tCore( tMesh, mObject, mVerboseFlag );

        // calculate SDF
        tCore.calculate_raycast_and_sdf( aSDF );
    }

    //-------------------------------------------------------------------------------
}    // namespace moris::sdf
