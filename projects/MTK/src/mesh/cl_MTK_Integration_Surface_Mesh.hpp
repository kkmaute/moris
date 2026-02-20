/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Surface_Mesh.hpp
 *
 */

#pragma once

#include "fn_MTK_Integration_Surface_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_XTK_Enums.hpp"
#include "moris_typedefs.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_Json_Object.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include <ostream>

namespace moris::mtk
{
    /**
     * @brief This class is used to extract (or better: to provide a view of) a surface mesh from a given integration mesh.
     * This class does not store any vertex or cell data but only provides the necessary information to access the data in the integration mesh.
     * The surface mesh will use local indices that only refer to the vertices/facets in the surface mesh.
     * A corresponding mapping between the global and local indices is provided.
     */
    class Integration_Surface_Mesh : public Surface_Mesh
    {

      public:    // constructors
        Integration_Surface_Mesh(
                Integration_Surface_Mesh_Data const &aData );

        // methods

        [[nodiscard]] Matrix< DDRMat > initialize_vertex_coordinates( Integration_Mesh const *aIGMesh );

        [[nodiscard]] Matrix< DDRMat > initialize_vertex_coordinates_from_side_sets(
                Integration_Mesh const           *aIGMesh,
                const Vector< Side_Set const * > &aSideSets );

        [[nodiscard]] Vector< Vector< moris_index > > get_cell_to_vertex_indices(
                Integration_Mesh const           *aIGMesh,
                const Vector< Side_Set const * > &aSideSets ) const;

        Vector< Side_Set const * >
        obtain_sidesets_from_names( Integration_Mesh_DataBase_IG const *aIGMesh, const Vector< std::string > &aSideSetNames );

        void set_all_displacements( const Matrix< DDRMat > &aDisplacements ) override;

        /**
         * @brief Returns the indices of all neighboring vertices for each vertex in the surface mesh.
         * @return A list of lists. The outer list contains the neighbor-lists for each vertex. The inner list contains the indices of the neighbors.
         */
        [[nodiscard]] Vector< Vector< moris_index > > get_vertex_neighbors() const;

        /**
         * @brief Returns the indices of all neighboring vertices for the vertex with the given local index.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return A list of indices of the neighbors of the vertex with the given index.
         */
        [[nodiscard]] Vector< moris_index > get_vertex_neighbors( moris_index aLocalVertexIndex ) const;


        /**
         * @brief Returns the facet measure (length/area) for each facet in the surface mesh.
         * @return A (n x 1) matrix where n is the number of facets in the surface mesh.
         */
        [[nodiscard]] Vector< real > compute_facet_measure() const final;

        /**
         * @brief Returns the averaged vertex normals for each vertex in the surface mesh.
         * @details The vertex normals are averaged over all facets that are connected to the vertex, weighted by the respective facet measure.
         * @return A (d x n) matrix where d is the dimension of the mesh (holding the normal components) and n is the number of vertices in the surface mesh.
         */
        [[nodiscard]] const Matrix< DDRMat > &get_vertex_normals() const;

        /**
         * @brief Returns the global index of a vertex with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalVertexIndex The local index of the vertex in the surface mesh.
         * @return The global index of the vertex in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_vertex_index( moris_index aLocalVertexIndex ) const;

        /**
         * @brief Returns the global index of a cell with the given local index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return The global index of the cell in the whole mesh.
         */
        [[nodiscard]] moris_index get_global_cell_index( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns the local index of a vertex with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalVertexIndex The global index of the vertex in the whole mesh.
         * @return The local index of the vertex in the surface mesh.
         */
        [[nodiscard]] moris_index get_local_vertex_index( moris_index aGlobalVertexIndex ) const;

        [[nodiscard]] moris_index get_cluster_of_cell( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns the (averaged) vertex normals of all vertices that are part of the cell with the given local index.
         * @param aLocalCellIndex The local index of the cell in the surface mesh.
         * @return A (d x n) matrix where d is the dimension of the mesh and n is the number of vertices in the cell.
         */
        [[nodiscard]] Matrix< DDRMat > get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const;

        /**
         * @brief Returns the local index of a cell with the given global index. Global refers to the whole mesh while local is only valid for the surface mesh.
         * @param aGlobalCellIndex The global index of the cell in the whole mesh.
         * @return The local index of the cell in the surface mesh.
         */
        // [[nodiscard]] moris_index get_local_cell_index( moris_index aGlobalCellIndex ) const;

        [[nodiscard]] uint get_spatial_dimension() const override;

        Json to_json() const;

        //--------------------------------------------------------------------------------
        // XQI Related functions
        //--------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------

        /**
         * Computes the sensitivities of the requested XQI type wrt to PDVs.
         *
         * @param aType The type of XQI for which the sensitivities are requested.
         * @param aVertexPDVIDs The PDV IDs associated with each vertex in the surface mesh. -1 if no PDV is associated with the vertex.
         * @param aSensitivities The distributed vector where the sensitivities will be stored. It is assumed that this vector is already initialized and has the correct map.
         * @param aRequestIndex The vector index in the Dist_Vector to store the sensitivities.
         * @param aExtra Extra arguments that may be required for specific XQI types (e.g., agglomeration functions).
         */
        template< typename... ExtraArgs >
        void compute_XQI_sensitivities(
                const xtk::XQI_Type                    aType,
                const Vector< Vector< moris_index > > &aVertexPDVIDs,
                sol::Dist_Vector                      *aSensitivities,
                const uint                             aRequestIndex,
                ExtraArgs &&...aExtra ) const
        {
            // Need a unified function that computes the sensitivity of the requested XQI wrt to a given vertex
            // takes only the local vertex index as input and returns a (d x 1) matrix with the sensitivity components
            std::function< Matrix< DDRMat >( uint ) > get_dXQI_dvertex = nullptr;

            switch ( aType )
            {
                case xtk::XQI_Type::VOLUME:
                    // no extra args required
                    get_dXQI_dvertex = [ this ]( uint aV ) -> Matrix< DDRMat > {
                        return this->compute_dvolume_dvertex( aV );
                    };
                    break;

                case xtk::XQI_Type::SHAPE_DIAMETER:
                {
                    // Capture extra args into a tuple
                    auto tExtras     = std::make_tuple( std::forward< ExtraArgs >( aExtra )... );
                    get_dXQI_dvertex = [ this, tExtras ]( uint aV ) -> Matrix< DDRMat > {
                        // apply the tuple to a helper that calls the member function with the extra args
                        return std::apply(
                                [ this, aV ]( auto &&...args ) -> Matrix< DDRMat > {
                                    return this->compute_ddiameter_dvertex( aV, std::forward< decltype( args ) >( args )... );
                                },
                                tExtras );
                    };
                }
                break;

                default:
                    MORIS_ERROR( false, "XQI type not implemented for surface mesh geometry." );
                    break;
            }

            MORIS_ASSERT( get_dXQI_dvertex, "Internal error: no callable assigned for XQI sensitivity" );

            // Loop over surface mesh vertices
            for ( uint iV = 0; iV < this->get_number_of_vertices(); iV++ )
            {
                // Check that this vertex has at least one PDV associated with it
                bool tHasPDV = false;
                for ( uint iDim = 0; iDim < aVertexPDVIDs.size(); iDim++ )
                {
                    if ( aVertexPDVIDs( iDim )( iV ) != -1 )
                    {
                        tHasPDV = true;
                        break;
                    }
                }

                if ( !tHasPDV )
                {
                    continue;
                }

                // Compute the sensitivity wrt to the vertex via the unified callable
                Matrix< DDRMat > tdXQI_dvertex = get_dXQI_dvertex( iV );

                // Sum into the distributed sensitivity vector
                for ( uint iDim = 0; iDim < aVertexPDVIDs.size(); iDim++ )
                {
                    moris_index tPDVID = aVertexPDVIDs( iDim )( iV );
                    if ( tPDVID != -1 )
                    {
                        real &tValue = ( *aSensitivities )( tPDVID, aRequestIndex );
                        tValue += tdXQI_dvertex( iDim );
                    }
                }
            }
        }

      private:    // methods
        void initialize_facet_measure();

        void initialize_vertex_normals();

        /**
         * @brief Contains information about the surface mesh including mapping to the original IG mesh and other useful maps
         */
        Integration_Surface_Mesh_Data mData;

        /**
         * @brief Stores the averaged vertex normals for each vertex in the surface mesh. The indices are the indices of the vertices in the surface mesh, not the global indices!
         */
        Matrix< DDRMat > mVertexNormals = Matrix< DDRMat >( 0, 0 );

        /**
         * @brief Is used to store the measure of each facet.
         */
        Vector< real > mFacetMeasure = Vector< real >( 0 );
    };
}    // namespace moris::mtk
