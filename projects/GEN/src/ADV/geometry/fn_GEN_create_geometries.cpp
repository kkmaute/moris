/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_geometries.cpp
 *
 */

#include "fn_GEN_create_geometries.hpp"
#include "fn_Parsing_Tools.hpp"

#include "fn_GEN_create_field.hpp"
#include "cl_GEN_Combined_Field.hpp"
#include "cl_GEN_Voxel_Input.hpp"

namespace moris::ge
{
    template< typename Vector_Type >
    inline Cell< std::shared_ptr< Level_Set_Geometry > >
    create_geometries(
            Cell< ParameterList >         aGeometryParameterLists,
            Vector_Type&                  aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMTKMesh )
    {
        // Create geometry cell
        Cell< std::shared_ptr< Level_Set_Geometry > >      tGeometries( 0 );
        Cell< std::shared_ptr< Combined_Field > > tMultigeometries( 0 );

        // Create individual geometries
        for ( uint tGeometryIndex = 0; tGeometryIndex < aGeometryParameterLists.size(); tGeometryIndex++ )
        {
            // Create geometry
            std::shared_ptr< Level_Set_Geometry > tGeometry = create_geometry(
                    aGeometryParameterLists( tGeometryIndex ),
                    aADVs,
                    aLibrary,
                    nullptr,
                    aMTKMesh );

            // Determine if to add to multigeometry
            bool        tMultigeometryFound = false;
            std::string tGeometryName       = tGeometry->get_name();
            if ( tGeometryName != "" )
            {
                // Loop to see if this multigeometry ID exists already FIXME
//                for ( uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++ )
//                {
//                    if ( tMultigeometries( tMultigeometryIndex )->get_name() == tGeometryName )
//                    {
//                        tMultigeometryFound = true;
//                        tMultigeometries( tMultigeometryIndex )->add_geometry( tGeometry );
//                        break;
//                    }
//                }

                // Check for creating new multigeometry FIXME
//                if ( not tMultigeometryFound )
//                {
//                    for ( uint tCreatedGeometryIndex = 0; tCreatedGeometryIndex < tGeometries.size(); tCreatedGeometryIndex++ )
//                    {
//                        if ( tGeometries( tCreatedGeometryIndex )->get_name() == tGeometryName )
//                        {
//                            tMultigeometryFound = true;
//                            tMultigeometries.push_back( std::make_shared< Multigeometry >(
//                                    Cell< std::shared_ptr< Level_Set_Geometry > >( { tGeometries( tCreatedGeometryIndex ), tGeometry } ) ) );
//                            tGeometries.erase( tCreatedGeometryIndex );
//                            break;
//                        }
//                    }
//                }
            }

            // If no multigeometry, add as regular geometry
            if ( not tMultigeometryFound )
            {
                tGeometries.push_back( tGeometry );
            }

            // TODO generalize this
//            if ( aGeometryParameterLists( tGeometryIndex ).get< std::string >( "type" ) == "voxel" )
//            {
//                uint tNumVoxelIDs = reinterpret_cast< Voxel_Input* >( tGeometry.get() )->get_num_voxel_Ids();
//
//                for ( uint tVoxelID = 1; tVoxelID <= tNumVoxelIDs; tVoxelID++ )
//                {
//                    std::shared_ptr< Level_Set_Geometry > tGeometrySingleGrain =
//                            create_geometry(
//                                    aGeometryParameterLists( tGeometryIndex ),
//                                    aADVs,
//                                    aLibrary,
//                                    tGeometry,
//                                    aMTKMesh,
//                                    tVoxelID );
//
//                    tGeometries.push_back( tGeometrySingleGrain );
//                }
//            }
        }

        // Add multigeometries at the end FIXME
//        for ( uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++ )
//        {
//            tGeometries.push_back( tMultigeometries( tMultigeometryIndex ) );
//        }

        return tGeometries;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Vector_Type >
    std::shared_ptr< Level_Set_Geometry >
    create_geometry(
            ParameterList                 aGeometryParameterList,
            Vector_Type&                  aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            std::shared_ptr< Level_Set_Geometry >   aGeometry,
            mtk::Mesh*                    aMTKMesh,
            uint                          aIndex )
    {
        // Create field
        std::shared_ptr< Field > tField = create_field( aGeometryParameterList, aADVs, {}, aLibrary, aMTKMesh, aIndex );

        // Geometry parameters
        Level_Set_Parameters tParameters;
        tParameters.mName                     = aGeometryParameterList.get< std::string >( "name" );
        tParameters.mNumRefinements           = aGeometryParameterList.get< std::string >( "number_of_refinements" );
        tParameters.mRefinementMeshIndices    = aGeometryParameterList.get< std::string >( "refinement_mesh_index" );
        tParameters.mRefinementFunctionIndex  = aGeometryParameterList.get< sint >( "refinement_function_index" );
        tParameters.mDiscretizationIndex      = aGeometryParameterList.get< sint >( "discretization_mesh_index" );
        tParameters.mDiscretizationLowerBound = aGeometryParameterList.get< real >( "discretization_lower_bound" );
        tParameters.mDiscretizationUpperBound = aGeometryParameterList.get< real >( "discretization_upper_bound" );
        tParameters.mIsocontourThreshold      = aGeometryParameterList.get< real >( "isocontour_threshold" );
        tParameters.mIsocontourTolerance      = aGeometryParameterList.get< real >( "isocontour_tolerance" );
        tParameters.mIntersectionTolerance    = aGeometryParameterList.get< real >( "intersection_tolerance" );
        if ( aGeometryParameterList.get< bool >( "multilinear_intersections" ) )
        {
            tParameters.mIntersectionInterpolation = Int_Interpolation::MULTILINEAR;
        }

        // Create geometry
        return std::make_shared< Level_Set_Geometry >( tField, tParameters );
    }

    //--------------------------------------------------------------------------------------------------------------
    // Explicit template instantiation
    //--------------------------------------------------------------------------------------------------------------

    template Cell< std::shared_ptr< Level_Set_Geometry > > create_geometries(
            Cell< ParameterList >         aGeometryParameterLists,
            Matrix< DDRMat >&             aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMTKMesh );

    template Cell< std::shared_ptr< Level_Set_Geometry > > create_geometries(
            Cell< ParameterList >         aGeometryParameterLists,
            sol::Dist_Vector*&            aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMTKMesh );

    //--------------------------------------------------------------------------------------------------------------

}
