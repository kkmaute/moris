/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * fn_MTK_QuadraturePointMapper_Ray_ArborX_Details.hpp
 *
 */
#pragma once
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_MappingResult.hpp"

#include <ArborX.hpp>
#include <ArborX_Box.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_View.hpp>
#include <decl/Kokkos_Declare_SERIAL.hpp>
#include <functional>
#include <iostream>
#include <sys/types.h>
#include <tuple>
#include <unordered_map>

using moris::moris_index;

namespace moris::mtk::arborx
{
    template< typename T >
    using index_map = std::unordered_map< moris_index, T >;

    /**
     * @brief Datastructure to organize and store the mapping of multiple rays onto cells (boxes) in multiple meshes.
     * The first index is the mesh index, the second index is the cell index in that mesh (assuming that it has a continuous numbering as in the mtk::SurfaceMesh).
     * The values are the indices of the rays that intersect with the respective cell.
     */
    using cell_locator_map = index_map< index_map< moris::Vector< moris_index > > >;

    template< typename MemorySpace >
    struct QueryBoxes
    {
        explicit QueryBoxes(
                Kokkos::View< ArborX::Box *, MemorySpace > aBoxes,
                Kokkos::View< moris_index *, MemorySpace > aMeshIndices,
                Kokkos::View< moris_index *, MemorySpace > aCellIndices )
                : mBoxes( aBoxes )
                , mMeshIndices( aMeshIndices )
                , mCellIndices( aCellIndices )
        {
        }
        KOKKOS_FUNCTION
        [[nodiscard]] std::size_t size() const
        {
            return mBoxes.extent( 0 );
        }
        KOKKOS_FUNCTION
        ArborX::Box const &operator()( std::size_t i ) const
        {
            return mBoxes( i );
        }

        Kokkos::View< ArborX::Box *, MemorySpace > mBoxes;
        Kokkos::View< moris_index *, MemorySpace > mMeshIndices;
        Kokkos::View< moris_index *, MemorySpace > mCellIndices;
    };

    template< typename MemorySpace >
    struct QueryRays
    {
        explicit QueryRays( Kokkos::View< ArborX::Experimental::Ray *, MemorySpace > aRays,
                Kokkos::View< moris_index *, MemorySpace >                           aCellIndices = {} )
                : mRays( aRays )
                , mCellIndices( aCellIndices )
        {
        }

        KOKKOS_FUNCTION
        [[nodiscard]] std::size_t size() const
        {
            return mRays.extent( 0 );
        }

        KOKKOS_FUNCTION
        ArborX::Experimental::Ray const &operator()( std::size_t i ) const
        {
            return mRays( i );
        }
        Kokkos::View< ArborX::Experimental::Ray *, MemorySpace > mRays;
        Kokkos::View< moris_index *, MemorySpace >               mCellIndices;
    };

    struct QueryResult
    {
        moris_index mPointIndex;
        //        moris_index mRayIndex;
        moris_index mBoxIndex;
    };


    template< typename MemorySpace >
    struct RayIntersectionCallback
    {
        template< typename Predicate, typename OutputFunctor >
        KOKKOS_FUNCTION void operator()( Predicate const &predicate, int const primitive_index, OutputFunctor const &out ) const
        {
            int const predicate_index = ArborX::getData( predicate );
            //            uint const tBocCellIndex   = mQueryBoxes.mCellIndices( primitive_index );
            //            uint const tBoxMeshIndex   = mQueryBoxes.mMeshIndices( primitive_index );
            //            uint const tRayCellIndex   = mQueryRays.mCellIndices( predicate_index );

            moris_index const tPointIndex = predicate_index;
            //            std::cout << "Intersection found between ray " << predicate_index << " (" << tRayCellIndex << ") and box " << primitive_index << " (" << tBocCellIndex << " on mesh " << tBoxMeshIndex << ")" << std::endl;
            out( QueryResult{ tPointIndex, primitive_index } );
        }
        QueryRays< MemorySpace > mQueryRays;
    };

    template< typename MemorySpace, typename ExecutionSpace >
    QueryBoxes< MemorySpace > construct_query_boxes( ExecutionSpace const &aExecutionSpace, moris::Vector< std::pair< moris_index, moris::mtk::Surface_Mesh > > const &aTargetSurfaceMeshes );

    template< typename MemorySpace, typename ExecutionSpace >
    QueryRays< MemorySpace > construct_query_rays( ExecutionSpace const &iRay, moris::mtk::MappingResult const &aMappingResult );

    //    /**
    //     * @brief To be able to use the tuple as a key in the unordered_map, we need to define a hash function for it.
    //     * The logic is based on https://stackoverflow.com/questions/20834838/using-tuple-in-unordered-map and the cantor pairing function
    //     * from https://stackoverflow.com/questions/38965931/hash-function-for-3-integers.
    //     */
    //    struct cell_locator_hash : public std::unary_function< cell_locator_tuple, std::size_t >
    //    {
    //        static std::size_t cantor_pairing( std::size_t a, std::size_t b )
    //        {
    //            return ( a + b ) * ( a + b + 1 ) / 2 + b;
    //        }
    //        std::size_t operator()( const cell_locator_tuple &k ) const
    //        {
    //            return cantor_pairing( std::get< 0 >( k ), cantor_pairing( std::get< 1 >( k ), std::get< 2 >( k ) ) );
    //        }
    //    };

    //    std::unordered_map< cell_locator_tuple, moris::Vector< moris_index >, cell_locator_hash >
    cell_locator_map
    map_rays_to_boxes(
            moris::mtk::MappingResult const                                           &aMappingResult,
            moris::Vector< std::pair< moris_index, moris::mtk::Surface_Mesh > > const &aTargetSurfaceMeshes );

}    // namespace moris::mtk::arborx

namespace ArborX
{
    using moris::mtk::arborx::QueryBoxes;
    using moris::mtk::arborx::QueryRays;

    template< typename MemorySpace >
    struct AccessTraits< QueryRays< MemorySpace >, PredicatesTag >
    {
        using memory_space = MemorySpace;

        // Function to return the number of queries
        static KOKKOS_FUNCTION std::size_t size( QueryRays< MemorySpace > const &rays )
        {
            return rays.size();
        }

        // Function to construct a predicate from a query at a given index
        static KOKKOS_FUNCTION auto get( QueryRays< MemorySpace > const &rays, std::size_t i )
        {
            return attach( intersects( rays( i ) ), i );    // returns this predicate value and attaches the predicates index as extra data to it.
        }
    };

    template< typename MemorySpace >
    struct AccessTraits< QueryBoxes< MemorySpace >, PrimitivesTag >
    {
        using memory_space = MemorySpace;

        // Function to return the number of primitives
        static KOKKOS_FUNCTION std::size_t size( QueryBoxes< MemorySpace > const &primitives )
        {
            return primitives.size();
        }

        // Function to return a primitive at a given index
        static KOKKOS_FUNCTION auto get( QueryBoxes< MemorySpace > const &primitives, std::size_t i )
        {
            return primitives( i );    // Assuming primitives(i) returns a Box
        }
    };
}    // namespace ArborX
