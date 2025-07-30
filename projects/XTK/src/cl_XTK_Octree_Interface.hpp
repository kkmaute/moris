/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Octree_Interface.hpp
 *
 */

#ifndef SRC_cl_XTK_Octree_Interface
#define SRC_cl_XTK_Octree_Interface

#include "cl_XTK_Decomposition_Algorithm.hpp"
#include <functional>
#include <cmath>

namespace moris::xtk
{
    inline Vector< moris::real >
    linspace(
            moris::real tLowerVal,
            moris::real tEndIn,
            moris::lint tNumIn )
    {

        Vector< moris::real > linspaced;

        const moris::real tStart = tLowerVal;
        const moris::real tEnd   = tEndIn;
        const moris::lint tNum   = tNumIn;

        if ( tNum == 0 )
        {
            return linspaced;
        }

        if ( tNum == 1 )
        {
            linspaced.push_back( tStart );
            return linspaced;
        }

        moris::real tDelta = ( tEnd - tStart ) / ( tNum - 1 );

        for ( moris::lint i = 0; i < tNum - 1; ++i )
        {
            linspaced.push_back( tStart + tDelta * i );
        }

        linspaced.push_back( tEnd );
        return linspaced;
    }

    class IJK
    {
      public:
        IJK(
                const lint aI,
                const lint aJ,
                const lint aK )
                : mI( aI )
                , mJ( aJ )
                , mK( aK )
        {
        }
        lint mI;
        lint mJ;
        lint mK;

        lint
        get_I() const
        {
            return mI;
        }

        lint
        get_J() const
        {
            return mJ;
        }

        lint
        get_K() const
        {
            return mK;
        }
    };
    // https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1ijkMesh.html

    class IJK_Mesh
    {
      public:
        const moris::uint mNumVertX;
        const moris::uint mNumVertY;
        const moris::uint mNumVertZ;

      public:
        IJK_Mesh()
                : mNumVertX( 0 )
                , mNumVertY( 0 )
                , mNumVertZ( 0 ) {};
        IJK_Mesh(
                const moris::uint aNumX,
                const moris::uint aNumY,
                const moris::uint aNumZ )
                : mNumVertX( aNumX )
                , mNumVertY( aNumY )
                , mNumVertZ( aNumZ )
        {
        }

        //
        Vertex_Ancestry
        get_vertex_parent_entities() const;

        uint
        num_vert_x() const
        {
            return mNumVertX;
        }

        uint
        num_vert_y() const
        {
            return mNumVertY;
        }

        uint
        num_vert_z() const
        {
            return mNumVertZ;
        }

        uint
        num_verts() const
        {
            return mNumVertX * mNumVertY * mNumVertZ;
        }

        moris_index
        get_vertex_index(
                moris_index aI,
                moris_index aJ,
                moris_index aK ) const
        {
            return ( aI + ( ( this->num_vert_x() ) * ( aJ + ( this->num_vert_y() ) * aK ) ) );
        }

        IJK
        get_vertex_ijk( moris_index aVertexIndex ) const
        {
            moris::lint tI = aVertexIndex % this->num_vert_x();
            moris::lint tK = aVertexIndex / ( this->num_vert_x() * this->num_vert_y() );
            moris::lint tJ = ( aVertexIndex / this->num_vert_x() ) - tK * this->num_vert_y();

            return IJK( tI, tJ, tK );
        }

        lint
        min_vert_i() const
        {
            return 0;
        }

        lint
        max_vert_i() const
        {
            return this->num_vert_x() - 1;
        }

        lint
        min_vert_j() const
        {
            return 0;
        }

        lint
        max_vert_j() const
        {
            return this->num_vert_y() - 1;
        }

        lint
        min_vert_k() const
        {
            return 0;
        }

        lint
        max_vert_k() const
        {
            return this->num_vert_z() - 1;
        }

        void
        get_vertex_parent(
                IJK const   &aVertIJK,
                moris_index &aParentOrdinal,
                moris_index &aParentRank ) const
        {
            // mark internal vertices
            if ( this->is_internal_vertex( aVertIJK ) )
            {
                aParentOrdinal = 0;
                aParentRank    = 3;
                return;
            }
        }

        bool
        is_internal_vertex( IJK const &aVertIJK ) const
        {
            if ( ( aVertIJK.mI != this->min_vert_i() && aVertIJK.mI != this->max_vert_i() ) && ( aVertIJK.mJ != this->min_vert_j() && aVertIJK.mJ != this->max_vert_j() ) && ( aVertIJK.mK != this->min_vert_k() && aVertIJK.mK != this->max_vert_k() ) )
            {
                return true;
            }

            return false;
        }

        uint
        num_cell_x() const
        {
            return mNumVertX - 1;
        }

        uint
        num_cell_y() const
        {
            return mNumVertY - 1;
        }

        uint
        num_cell_z() const
        {
            return mNumVertZ - 1;
        }

        uint
        num_cells() const
        {
            return num_cell_x() * num_cell_y() * num_cell_z();
        }

        uint
        num_verts_per_cell() const
        {
            return 8;
        }

        moris_index
        get_cell_index(
                moris_index aI,
                moris_index aJ,
                moris_index aK ) const
        {
            return ( aI + ( ( this->num_cell_x() ) * ( aJ + ( this->num_cell_y() ) * aK ) ) );
        }

        IJK
        get_cell_ijk( moris_index aCellIndex ) const
        {
            moris::lint tI = aCellIndex % this->num_cell_x();
            moris::lint tK = aCellIndex / ( this->num_cell_x() * this->num_cell_y() );
            moris::lint tJ = ( aCellIndex / this->num_cell_x() ) - tK * this->num_cell_y();

            return IJK( tI, tJ, tK );
        }

        Matrix< IndexMat >
        get_cell_to_vertex() const
        {
            // connectivity of cell i,j,k
            Matrix< IndexMat > tCellToNode( this->num_cells(), this->num_verts_per_cell() );

            for ( moris::uint k = 0; k < this->num_cell_z(); k++ )
            {
                for ( moris::uint j = 0; j < this->num_cell_y(); j++ )
                {
                    for ( moris::uint i = 0; i < this->num_cell_x(); i++ )
                    {
                        const moris::uint tIndex = this->get_cell_index( i, j, k );

                        tCellToNode( tIndex, 0 ) = this->get_vertex_index( i, j, k );
                        tCellToNode( tIndex, 1 ) = this->get_vertex_index( i + 1, j, k );
                        tCellToNode( tIndex, 2 ) = this->get_vertex_index( i + 1, j + 1, k );
                        tCellToNode( tIndex, 3 ) = this->get_vertex_index( i, j + 1, k );
                        tCellToNode( tIndex, 4 ) = this->get_vertex_index( i, j, k + 1 );
                        tCellToNode( tIndex, 5 ) = this->get_vertex_index( i + 1, j, k + 1 );
                        tCellToNode( tIndex, 6 ) = this->get_vertex_index( i + 1, j + 1, k + 1 );
                        tCellToNode( tIndex, 7 ) = this->get_vertex_index( i, j + 1, k + 1 );
                    }
                }
            }

            return tCellToNode;
        }

        uint
        num_faces() const
        {
            if ( this->num_cells() == 0 )
            {
                return 0;
            }
            const moris_index tNumCellsX = this->num_cell_x();
            const moris_index tNumCellsY = this->num_cell_y();
            const moris_index tNumCellsZ = this->num_cell_z();

            return ( ( tNumCellsX + 1 ) * tNumCellsY * tNumCellsZ )
                 + ( ( tNumCellsY + 1 ) * tNumCellsZ * tNumCellsX )
                 + ( ( tNumCellsZ + 1 ) * tNumCellsX * tNumCellsY );
        }

        uint
        num_faces_internal() const
        {
            if ( this->num_cells() == 0 )
            {
                return 0;
            }

            const moris_index tNumCellsX = this->num_cell_x();
            const moris_index tNumCellsY = this->num_cell_y();
            const moris_index tNumCellsZ = this->num_cell_z();

            return ( ( tNumCellsX - 1 ) * tNumCellsY * tNumCellsZ )
                 + ( ( tNumCellsY - 1 ) * tNumCellsZ * tNumCellsX )
                 + ( ( tNumCellsZ - 1 ) * tNumCellsX * tNumCellsY );
        }

        bool
        check_ijk_to_vertex_index() const
        {
            for ( moris::uint k = 0; k < this->num_vert_z(); k++ )
            {
                for ( moris::uint j = 0; j < this->num_vert_y(); j++ )
                {
                    for ( moris::uint i = 0; i < this->num_vert_x(); i++ )
                    {
                        moris_index tIndex = this->get_vertex_index( i, j, k );
                        IJK         tIJK   = this->get_vertex_ijk( tIndex );

                        if ( (moris_index)i != tIJK.mI )
                        {
                            return false;
                        }
                        if ( (moris_index)j != tIJK.mJ )
                        {
                            return false;
                        }
                        if ( (moris_index)k != tIJK.mK )
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }

        bool
        check_ijk_to_cell_index() const
        {
            for ( moris::uint k = 0; k < this->num_cell_z(); k++ )
            {
                for ( moris::uint j = 0; j < this->num_cell_y(); j++ )
                {
                    for ( moris::uint i = 0; i < this->num_cell_x(); i++ )
                    {
                        moris_index tIndex = this->get_cell_index( i, j, k );
                        IJK         tIJK   = this->get_cell_ijk( tIndex );

                        if ( (moris_index)i != tIJK.mI )
                        {
                            return false;
                        }
                        if ( (moris_index)j != tIJK.mJ )
                        {
                            return false;
                        }
                        if ( (moris_index)k != tIJK.mK )
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
    };

    class Octree_Template
    {
        // backend mesh for template (very simple ijk mesh)
        const moris_index           mLevel = -1;
        const IJK_Mesh              mOctreeMeshGrid;
        const Matrix< DDRMat >      mOctreeParamCoords;
        const Matrix< IndexMat >    mOctreeCells;
        const Vertex_Ancestry       mVertexAncestry;
        const Vector< moris_index > mVertexHash;

      public:
        Octree_Template() {};
        Octree_Template( const moris_index aTemplateLevel )
                : mLevel( aTemplateLevel )
                , mOctreeMeshGrid( this->create_mesh_grid( aTemplateLevel ) )
                , mOctreeParamCoords( this->generate_octree_param_coords( aTemplateLevel ) )
                , mOctreeCells( mOctreeMeshGrid.get_cell_to_vertex() )
                , mVertexAncestry( mOctreeMeshGrid.get_vertex_parent_entities() )
                , mVertexHash( this->generate_vertex_hashes() )
        {
            MORIS_ASSERT( mOctreeMeshGrid.check_ijk_to_vertex_index(), "Issue with octree numbering" );
            MORIS_ASSERT( mOctreeMeshGrid.check_ijk_to_cell_index(), "Issue with octree numbering" );
        }

        IJK_Mesh
        create_mesh_grid( const moris_index aTemplateLevel )
        {
            moris::uint tNum1DVerts = 1 + std::pow( 2, aTemplateLevel );

            return IJK_Mesh( tNum1DVerts, tNum1DVerts, tNum1DVerts );
        }

        const IJK_Mesh *
        get_mesh_grid() const
        {
            return &mOctreeMeshGrid;
        }

        const Vector< moris_index > *
        get_vertex_hashes() const
        {
            return &mVertexHash;
        }

        const Matrix< DDRMat > &
        get_vertex_param_coords() const
        {
            return mOctreeParamCoords;
        }

        const Vertex_Ancestry *
        get_vertex_ancestry() const
        {
            return &mVertexAncestry;
        }
        const Matrix< IndexMat > *
        get_cells() const
        {
            return &mOctreeCells;
        }

        moris_index
        get_level() const
        {
            return mLevel;
        }

        Matrix< DDRMat >
        generate_octree_param_coords( const moris_index aTemplateLevel )
        {
            // generate a linspace from -1.0 to 1.0 - this corresponds to the local coords of a hex or quad
            Vector< moris::real > const tLinSpace = xtk::linspace( -1.0, 1.0, mOctreeMeshGrid.num_vert_x() );

            // allocate the output matrix
            Matrix< DDRMat > tParamCoords( mOctreeMeshGrid.num_verts(), 3 );

            // fill the matrix
            for ( moris::uint k = 0; k < mOctreeMeshGrid.num_vert_z(); k++ )
            {
                for ( moris::uint j = 0; j < mOctreeMeshGrid.num_vert_y(); j++ )
                {
                    for ( moris::uint i = 0; i < mOctreeMeshGrid.num_vert_x(); i++ )
                    {
                        moris_index tIndex        = mOctreeMeshGrid.get_vertex_index( i, j, k );
                        tParamCoords( tIndex, 0 ) = tLinSpace( i );
                        tParamCoords( tIndex, 1 ) = tLinSpace( j );
                        tParamCoords( tIndex, 2 ) = tLinSpace( k );
                    }
                }
            }
            return tParamCoords;
        }

        moris_index
        get_num_corner_points() const
        {
            std::cout << "WARNING NEEDs ABSTRACTION" << '\n';
            return 8;
        }

        Vector< moris_index >
        generate_vertex_hashes()
        {
            Vector< moris_index > tVertexHashes( mOctreeMeshGrid.num_verts() );

            std::cout << "WARNING TODO MAKE PROPER VERTEX HASHES" << '\n';
            // I need to make these hashes in a way that makes consistent requests across facets

            MORIS_ERROR( mOctreeMeshGrid.max_vert_i() < 10000 && mOctreeMeshGrid.max_vert_j() < 10000 && mOctreeMeshGrid.max_vert_k() < 10000, "Current limit imposed by hashing in octree" );

            const Vertex_Ancestry *tVertexAncestry = this->get_vertex_ancestry();

            // iterate through vertices
            for ( moris::uint iVertex = 0; iVertex < mOctreeMeshGrid.num_verts(); iVertex++ )
            {
                const IJK tIJK = mOctreeMeshGrid.get_vertex_ijk( iVertex );

                moris_index tIInt = (moris_index)std::trunc( 10000 * mOctreeParamCoords( iVertex, 0 ) );
                moris_index tJInt = (moris_index)std::trunc( 10000 * mOctreeParamCoords( iVertex, 1 ) );
                moris_index tKInt = (moris_index)std::trunc( 10000 * mOctreeParamCoords( iVertex, 2 ) );

                // std::cout << " mOctreeParamCoords( iVertex, 0 ) =  " << mOctreeParamCoords( iVertex, 0 ) << " | tIInt = " << tIInt << std::endl;

                // iternal vertex hash
                if ( tVertexAncestry->get_vertex_parent_rank( iVertex ) == mtk::EntityRank::ELEMENT )
                {
                    // tVertexHashes( iVertex ) = 100 * tIJK.mI + 10000 * tIJK.mJ + 100000000 * tIJK.mK;
                    tVertexHashes( iVertex ) = std::hash< int >{}( tIInt + 1000 * tJInt + 100000 * tKInt );
                }
                else if ( tVertexAncestry->get_vertex_parent_rank( iVertex ) == mtk::EntityRank::NODE )
                {
                    // intentionally left out, so we don't hash vertices that are already in the mesh
                }

                else if ( tVertexAncestry->get_vertex_parent_rank( iVertex ) == mtk::EntityRank::FACE )
                {

                    if ( tIJK.mI == mOctreeMeshGrid.min_vert_i() || tIJK.mI == mOctreeMeshGrid.max_vert_i() )
                    {
                        // create a hash using only the k and j components  and a fixed i signature - 1
                        // tVertexHashes( iVertex ) = 1 + 100 * tIJK.mJ + 1000000 * tIJK.mK;

                        tVertexHashes( iVertex ) = std::hash< int >{}( 1 + 10 * tJInt + 10000 * tKInt );
                    }

                    if ( tIJK.mJ == mOctreeMeshGrid.min_vert_j() || tIJK.mJ == mOctreeMeshGrid.max_vert_j() )
                    {
                        // create a hash using only the k and j components  and a fixed i signature - 1
                        // tVertexHashes( iVertex ) = 2 + 100 * tIJK.mI + 1000000 * tIJK.mK;
                        tVertexHashes( iVertex ) = std::hash< int >{}( 2 + 10 * tIInt + 10000 * tKInt );
                    }

                    if ( tIJK.mK == mOctreeMeshGrid.min_vert_k() || tIJK.mK == mOctreeMeshGrid.max_vert_k() )
                    {
                        // create a hash using only the k and j components  and a fixed i signature - 1
                        // tVertexHashes( iVertex ) = 3 + 100 * tIJK.mI + 1000000 * tIJK.mJ;
                        tVertexHashes( iVertex ) = std::hash< int >{}( 3 + 10 * tIInt + 10000 * tJInt );
                    }
                }
                else if ( tVertexAncestry->get_vertex_parent_rank( iVertex ) == mtk::EntityRank::EDGE )
                {
                    // hash edges at the are invariant to j k - edges 0, 2, 4 6
                    if ( ( tIJK.mJ == mOctreeMeshGrid.min_vert_j() || tIJK.mJ == mOctreeMeshGrid.max_vert_j() ) && ( tIJK.mK == mOctreeMeshGrid.min_vert_k() || tIJK.mK == mOctreeMeshGrid.max_vert_k() ) )
                    {
                        // tVertexHashes( iVertex ) = 4 + 100 * tIJK.mI;
                        tVertexHashes( iVertex ) = std::hash< int >{}( 4 + 100 * tIInt );
                    }

                    // hash edges at the are invariant to i k - edges 1,3,5,7
                    if ( ( tIJK.mI == mOctreeMeshGrid.min_vert_i() || tIJK.mI == mOctreeMeshGrid.max_vert_i() ) && ( tIJK.mK == mOctreeMeshGrid.min_vert_k() || tIJK.mK == mOctreeMeshGrid.max_vert_k() ) )
                    {
                        // tVertexHashes( iVertex ) = 5 + 100 * tIJK.mJ;
                        tVertexHashes( iVertex ) = std::hash< int >{}( 5 + 100 * tJInt );
                    }
                    // hash edges at the are invariant to i j - edges 8 9 10 11
                    if ( ( tIJK.mI == mOctreeMeshGrid.min_vert_i() || tIJK.mI == mOctreeMeshGrid.max_vert_i() ) && ( tIJK.mJ == mOctreeMeshGrid.min_vert_j() || tIJK.mJ == mOctreeMeshGrid.max_vert_j() ) )
                    {
                        // tVertexHashes( iVertex ) = 6 + 100 * tIJK.mK;
                        tVertexHashes( iVertex ) = std::hash< int >{}( 6 + 100 * tKInt );
                    }
                }
            }

            return tVertexHashes;
        }
    };

    class Octree_Comparator
    {
      public:
        const Octree_Template *mOctreeTemplateCoarse;
        const Octree_Template *mOctreeTemplateFine;
        const moris_index      mCachedLevelDiff;
        const moris_index      mCachedPower;

        Octree_Comparator(
                const Octree_Template *aOctreeTemplateCoarse,
                const Octree_Template *aOctreeTemplateFine )
                : mOctreeTemplateCoarse( aOctreeTemplateCoarse )
                , mOctreeTemplateFine( aOctreeTemplateFine )
                , mCachedLevelDiff( mOctreeTemplateFine->get_level() - mOctreeTemplateCoarse->get_level() )
                , mCachedPower( std::pow( 2, mCachedLevelDiff ) )
        {
            MORIS_ERROR( mOctreeTemplateCoarse->get_level() <= mOctreeTemplateFine->get_level(), "Coarse to fine mismatch" );
        }
        /**
         * @brief Translates IJK from template 0 to template 1
         *
         * @param aIJK Coarse IJK
         * @return IJK Fine IJK
         */

        IJK
        translate( const IJK &aCoarseIJK ) const
        {

            return IJK(
                    aCoarseIJK.mI * mCachedPower,
                    aCoarseIJK.mJ * mCachedPower,
                    aCoarseIJK.mK * mCachedPower );
        }

        bool
        vertex_exists_in_coarse( const IJK &aFineIJK ) const
        {
            if ( aFineIJK.mI % mCachedPower == 0 && aFineIJK.mJ % mCachedPower == 0 && aFineIJK.mK % mCachedPower == 0 )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        /**
         * @brief List of IJK positions unique to the fine grid
         *
         * @return Cell<IJK>  container of IJK positions in fine grid
         */
        Vector< IJK >
        unique_to_fine_grid() const
        {
            Vector< IJK > tCell;
            tCell.reserve( mOctreeTemplateFine->get_mesh_grid()->num_verts() - mOctreeTemplateCoarse->get_mesh_grid()->num_verts() );

            // iterate through fine grid points
            for ( moris::uint k = 0; k < mOctreeTemplateFine->get_mesh_grid()->num_vert_z(); k++ )
            {
                for ( moris::uint j = 0; j < mOctreeTemplateFine->get_mesh_grid()->num_vert_y(); j++ )
                {
                    for ( moris::uint i = 0; i < mOctreeTemplateFine->get_mesh_grid()->num_vert_x(); i++ )
                    {
                        IJK tIJK( (moris_index)i, (moris_index)j, (moris_index)k );
                        if ( !this->vertex_exists_in_coarse( tIJK ) )
                        {
                            tCell.push_back( tIJK );
                        }
                    }
                }
            }
            return tCell;
        }
    };

    class Octree_Interface : public Decomposition_Algorithm
    {
      private:
        moris::lint                                              mOctreeRefinementLevel;
        Vector< std::shared_ptr< Octree_Template const > >       mOctreeTemplates;
        Vector< moris_index >                                    mDifference;
        Vector< std::unordered_map< moris_index, moris_index > > mIgVertexGroupIndexToIjkIndex;

        // useful data that is passed in (avoid passing around to every function)
        Integration_Mesh_Generation_Data *mIgMeshGenData;
        Decomposition_Data               *mDecompData;
        Cut_Integration_Mesh             *mCutIgMesh;
        moris::mtk::Mesh                 *mBackgroundMesh;
        Integration_Mesh_Generator       *mIgMeshGen;

      public:
        Octree_Interface( Parameter_List &aParameterList );

        // set of
        void
        perform(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator ) override;

        enum Decomposition_Algorithm_Type
        get_algorithm_type() const override;

        moris_index
        get_signature() const override;

        bool
        has_geometric_dependent_vertices() const override;

        Vector< moris_index > get_decomposed_cell_indices() override;

        void
        perform_impl_vertex_requests(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator ) override;

        void
        perform_impl_generate_mesh(
                Integration_Mesh_Generation_Data *aMeshGenerationData,
                Decomposition_Data               *aDecompositionData,
                Cut_Integration_Mesh             *aCutIntegrationMesh,
                moris::mtk::Mesh                 *aBackgroundMesh,
                Integration_Mesh_Generator       *aMeshGenerator ) override;

      private:
        Vector< std::shared_ptr< Octree_Template const > >
        generate_octree_templates();

        /**
         * @brief Determine the lowest and highest element level that we need to refine
         *
         * @return Vector<moris_index> const - {LowerBound,UpperBound}
         */
        Vector< moris_index >
        determine_octree_bounds();

        mtk::CellTopology
        get_ig_cell_topology() const;

        Vector< std::unordered_map< moris_index, moris_index > >
        generate_octree_template_vertex_group_to_ijk_map();
    };

}    // namespace moris::xtk

#endif /* cl_XTK_Octree_Interface.hpp */
