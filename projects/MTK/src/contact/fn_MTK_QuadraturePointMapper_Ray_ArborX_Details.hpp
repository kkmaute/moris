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
#include "cl_MTK_Surface_Mesh.hpp"

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
#include <utility>

using moris::moris_index;

// Forward declare surface mesh class
namespace moris::mtk
{
    class Surface_Mesh;
}

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

    /**
     * @brief Converts a moris::Matrix< moris::DDRMat > to an ArborX::Point or ArborX::Vector.
     * @tparam T The type of the ArborX object to be returned (Point or Vector)
     * @param aMatrix The matrix to be converted (either 3x1 or 2x1)
     * @return The converted ArborX object
     */
    template< typename T >
    T coordinate_to_arborx_point( moris::Matrix< moris::DDRMat > const &aMatrix );

    template< typename MemorySpace >
    struct QueryBoxes
    {
        explicit QueryBoxes(
                Kokkos::View< ArborX::Box *, MemorySpace > aBoxes,
                Kokkos::View< moris_index *, MemorySpace > aMeshIndices,
                Kokkos::View< moris_index *, MemorySpace > aCellIndices )
                : mBoxes( std::move( aBoxes ) )
                , mMeshIndices( std::move( aMeshIndices ) )
                , mCellIndices( std::move( aCellIndices ) )
        {
        }
        [[nodiscard]] KOKKOS_FUNCTION
                std::size_t
                size() const
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
                Kokkos::View< moris_index *, MemorySpace >                           aCellIndices )
                : mRays( std::move( aRays ) )
                , mCellIndices( std::move( aCellIndices ) )
        {
        }

        [[nodiscard]] KOKKOS_FUNCTION
                std::size_t
                size() const
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
    struct IntersectionCallback
    {
        template< typename Predicate, typename OutputFunctor >
        KOKKOS_FUNCTION void operator()( Predicate const &predicate, int const primitive_index, OutputFunctor const &out ) const
        {
            int const predicate_index = ArborX::getData( predicate );
            //            uint const tBocCellIndex   = mQueryBoxes.mCellIndices( primitive_index );
            //            uint const tBoxMeshIndex   = mQueryBoxes.mMeshIndices( primitive_index );
            //            uint const tRayCellIndex   = mQueryRays.mCellIndices( predicate_index );

            /**
             * since ArborX treats Rays as directional objects (i.e. they can only hit on the positive side of the normal), the rays are doubled. The first ray is the original ray and the second ray is the negative ray.
             * Therefore, the indices 0 and 1 both belong to the first point, 2 and 3 to the second point and so on...
             * Since the detailed check if and where the ray actually intersects the box does not depend on the direction, we are only interested in the index of the point, not the direction (i.e. which of the two rays hit).
             *
             * The index of the ray is the index of the predicate with the direction in either positive (p) or negative (n) direction.
             *   index:    0   1   2   3   4   5   6   7   8   9
             *   rays:     0p  0n  1p  1n  2p  2n  3p  3n  4p  4n
             *   points:   0   0   1   1   2   2   3   3   4   4
             *
             *  The point index can therefore be easily calculated from the ray index by dividing by 2 and rounding down (which happens automatically in C++!).
             */
            moris_index const tPointIndex = predicate_index / 2;
            //            std::cout << "Intersection found between ray " << predicate_index << " (" << tRayCellIndex << ") and box " << primitive_index << " (" << tBocCellIndex << " on mesh " << tBoxMeshIndex << ")" << std::endl;
            out( QueryResult{ tPointIndex, primitive_index } );
        }
        QueryRays< MemorySpace > mQueryRays;
    };

    template< typename MemorySpace, typename ExecutionSpace >
    QueryBoxes< MemorySpace > construct_query_boxes( ExecutionSpace const &aExecutionSpace, moris::Vector< std::pair< moris_index, moris::mtk::Surface_Mesh > > const &aTargetSurfaceMeshes );

    template< typename MemorySpace, typename ExecutionSpace >
    QueryRays< MemorySpace > construct_query_rays( ExecutionSpace const &aEcecutionSpace, moris::mtk::MappingResult const &aMappingResult );

    /**
     * @brief Constructs a QueryRays object from a set of origin and direction matrices.
     *
     * @param aOrigins Matrix of origin points for the rays. Each column is a point, size <dim> x <number of origins>
     * @param aDirections Matrix of directions for the rays. Each column is a direction, size <dim> x <number of directions>
     */
    template< typename MemorySpace, typename ExecutionSpace >
    QueryRays< MemorySpace > construct_query_rays_from_primitives(
            ExecutionSpace const &aExecutionSpace,
            Matrix< DDRMat >     &aOrigins,
            Matrix< DDRMat >     &aDirections )
    {
        // Get the number of origins, directions, and total number of rays
        uint const tNumPoints     = aOrigins.n_cols();
        uint const tNumDirections = aDirections.n_cols();
        uint const tNumRays       = tNumPoints * tNumDirections;

        // Initialize Kokkos views for the rays and cell indices to construct the output
        Kokkos::View< ArborX::Experimental::Ray *, MemorySpace > tRays( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), tNumRays );
        Kokkos::View< moris_index *, MemorySpace >               tCellIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), tNumRays );

        // Initialize the rays from the input primitives
        Kokkos::parallel_for(
                "initialize_rays",
                Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumRays ),
                KOKKOS_LAMBDA( size_t const iRay ) {
                    // Calculate the origin and direction indices for this ray
                    size_t const iOrigin    = iRay / tNumDirections;
                    size_t const iDirection = iRay % tNumDirections;

                    // Get the origin point for this ray
                    ArborX::Point const tOrigin = coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

                    // Get the direction for this ray
                    auto const &tNormal = aDirections.get_column( iDirection );

                    // Initialize the ray and set the cell index
                    tRays( iRay )        = ArborX::Experimental::Ray{ tOrigin, coordinate_to_arborx_point< ArborX::Experimental::Vector >( tNormal ) };
                    tCellIndices( iRay ) = 0;    // FIXME: update QueryRays struct to make this optional
                } );

        return QueryRays< MemorySpace >{ tRays, tCellIndices };
    }

    /**
     * @brief Constructs a QueryRays object from a set of origin and direction matrices.
     *
     * @param aOrigins Matrix of origin points for the rays. Each column is a point, size <dim> x <number of origins>
     * @param aDirections Matrices of directions for the rays. The Vector size must be equal to the number of origin points, or aOrigins.n_cols().
     * For the inner matrix, each column is a direction, size <dim> x <number of directions>
     */
    template< typename MemorySpace, typename ExecutionSpace >
    QueryRays< MemorySpace > construct_query_rays_from_primitives(
            ExecutionSpace const       &aExecutionSpace,
            Matrix< DDRMat >           &aOrigins,
            Vector< Matrix< DDRMat > > &aDirections )
    {
        // Get the number of origins
        uint const tNumPoints = aOrigins.n_cols();

        // Calculate the total number of rays and starting indices for each origin
        uint             tNumRays = 0;
        Vector< size_t > tStartIndices( tNumPoints );
        for ( uint i = 0; i < tNumPoints; ++i )
        {
            tStartIndices( i ) = tNumRays;
            tNumRays += aDirections( i ).n_cols();
        }

        // Initialize Kokkos views for the rays and cell indices to construct the output
        Kokkos::View< ArborX::Experimental::Ray *, MemorySpace > tRays( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:rays" ), tNumRays );
        Kokkos::View< moris_index *, MemorySpace >               tCellIndices( Kokkos::view_alloc( aExecutionSpace, Kokkos::WithoutInitializing, "view:cell_indices" ), tNumRays );

        {    // Initialize the rays from the input primitives
            // Loop over all the origins and directions
            Kokkos::parallel_for(
                    "initialize_rays",
                    Kokkos::RangePolicy< ExecutionSpace >( aExecutionSpace, 0, tNumPoints ),
                    KOKKOS_LAMBDA( size_t const iOrigin ) {
                        // Get the origin point for this set of rays
                        ArborX::Point const tOrigin = coordinate_to_arborx_point< ArborX::Point >( aOrigins.get_column( iOrigin ) );

                        // Get the number of directions for this origin
                        uint const tNumDirections = aDirections( iOrigin ).n_cols();

                        // Get the starting index for this origin
                        size_t const tStartIndex = tStartIndices( iOrigin );

                        // Loop over all the directions for this origin
                        for ( size_t iDirection = 0; iDirection < tNumDirections; ++iDirection )
                        {
                            // Get the index for this ray
                            size_t const tIndex = tStartIndex + iDirection;

                            // Get the direction for this ray
                            auto const &tNormal = aDirections( iOrigin ).get_column( iDirection );

                            tRays( tIndex )        = ArborX::Experimental::Ray{ tOrigin, coordinate_to_arborx_point< ArborX::Experimental::Vector >( tNormal ) };
                            tCellIndices( tIndex ) = 0;    // FIXME: update QueryRays struct to make this optional
                        }
                    } );
        }

        return QueryRays< MemorySpace >{ tRays, tCellIndices };
    }

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
