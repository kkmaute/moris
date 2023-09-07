/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Integration_Mesh.hpp
 *
 */

#ifndef MORIS_cl_XTK_Cut_Integration_Mesh_HPP_
#define MORIS_cl_XTK_Cut_Integration_Mesh_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Vertex_XTK_Impl.hpp"

#include "cl_XTK_Subphase_Group.hpp"

#include "cl_Tracer.hpp"

#include "cl_Communication_Tools.hpp"
#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace moris;

namespace xtk
{
    // ----------------------------------------------------------------------------------

    struct IG_Cell_Group
    {
        IG_Cell_Group( moris_index aNumCellsInGroup );

        IG_Cell_Group();

        void
        add_Cell( moris::mtk::Cell* aCell );

        moris_index
        get_cell_group_ordinal( moris_index aCell );

        void
        remove_cell( moris_index aCell );

        bool
        cell_is_in_group( moris_index aCell );

        void
        shift_indices( moris_index aCell );

        moris::Cell< moris::mtk::Cell* >               mIgCellGroup;
        std::unordered_map< moris_index, moris_index > mIgCellIndexToCellOrdinal;

    };    // struct IG_Cell_Group

    // ----------------------------------------------------------------------------------

    struct IG_Cell_Side_Group
    {
        IG_Cell_Side_Group( moris_index aEstimatedNumCells );

        moris::Cell< moris::mtk::Cell* > mIgCells;               // over allocated
        moris::Cell< moris_index >       mIgCellSideOrdinals;    // over allocated

    };    // struct IG_Cell_Side_Group

    // ----------------------------------------------------------------------------------

    struct IG_Cell_Double_Side_Group
    {
        IG_Cell_Double_Side_Group( moris_index aEstimatedNumCells );

        moris::Cell< moris::mtk::Cell* > mLeaderIgCells;               // over allocated
        moris::Cell< moris_index >       mLeaderIgCellSideOrdinals;    // over allocated

        moris::Cell< moris::mtk::Cell* > mFollowerIgCells;               // over allocated
        moris::Cell< moris_index >       mFollowerIgCellSideOrdinals;    // over allocated

        void
        print()
        {
            std::cout << "Number of Leaders:   " << mLeaderIgCells.size() << std::endl;
            std::cout << "Number of Followers: " << mFollowerIgCells.size() << std::endl;

            int tStrLen = std::string( "Lead Cell Id   | " ).size();

            std::cout << "Lead Cell Id   | ";
            std::cout << "Side Ord       | ";
            std::cout << "Follow Cell Id | ";
            std::cout << "Side Ord       " << std::endl;
            // iterate through pairs
            for ( uint i = 0; i < mLeaderIgCells.size(); i++ )
            {
                std::cout << std::setw( tStrLen ) << mLeaderIgCells( i )->get_id();
                std::cout << std::setw( tStrLen ) << mLeaderIgCellSideOrdinals( i );
                std::cout << std::setw( tStrLen ) << mFollowerIgCells( i )->get_id();
                std::cout << std::setw( tStrLen ) << mFollowerIgCellSideOrdinals( i ) << std::endl;
            }
        }

    };    // struct IG_Cell_Double_Side_Group

    // ----------------------------------------------------------------------------------

    struct IG_Vertex_Group
    {
      private:
        moris::Cell< moris::mtk::Vertex const * >          mIgVertexGroup;
        std::unordered_map< moris_index, moris_index >     mIgVertexIndexToVertexOrdinal;
        moris::Cell< std::shared_ptr< Matrix< DDRMat > > > mIgVertexLocalCoords;

      public:
        IG_Vertex_Group( moris_index aNumVerticesInGroup );

        // destructor explicitly
        ~IG_Vertex_Group() {}

        std::size_t
        size();

        void
        reserve( std::size_t aReserveSize );

        void
        add_vertex(
                moris::mtk::Vertex const *          aVertex,
                std::shared_ptr< Matrix< DDRMat > > aVertexLocalCoord );

        void
        add_vertex_local_coord_pointers();

        moris::mtk::Vertex const *
        get_vertex( moris_index aGroupVertexOrdinal );

        moris_index
        get_vertex_group_ordinal( moris_index aVertex );

        std::shared_ptr< Matrix< DDRMat > >
        get_vertex_local_coords( moris_index aVertex );

        uint
        get_vertex_local_coords_dim() const;

        bool
        vertex_is_in_group( moris_index aVertex );

        void
        remove_vertex( moris_index aVertex );

        void
        shift_indices( moris_index aVertex );

        void
        print();

    };    // struct IG_Vertex_Group

    // ----------------------------------------------------------------------------------

    struct Edge_Based_Connectivity
    {
        moris::Cell< moris::Cell< moris::mtk::Vertex* > > mEdgeVertices;             // input: edge || output: list of vertices on edge
        moris::Cell< moris::Cell< moris::mtk::Cell* > >   mEdgeToCell;               // input: edge || output: list of cells attached to edge
        moris::Cell< moris::Cell< moris_index > >         mEdgeToCellEdgeOrdinal;    // input: edge || output: ?
        moris::Cell< moris::Cell< moris_index > >         mCellToEdge;               // input: cell || output: list of edge indices on cell
    };

    // ----------------------------------------------------------------------------------

    struct Edge_Based_Ancestry
    {
        moris::Cell< moris_index > mEdgeParentEntityIndex;
        moris::Cell< moris_index > mEdgeParentEntityRank;
        moris::Cell< moris_index > mEdgeParentEntityOrdinalWrtBackgroundCell;
    };

    // ----------------------------------------------------------------------------------

    struct Vertex_Ancestry
    {
        Vertex_Ancestry() = default;
        Vertex_Ancestry(
                moris::Cell< moris_index > const &     aVertexParentEntityIndices,
                moris::Cell< mtk::EntityRank > const & aVertexParentEntityRank )
                : mVertexParentEntityIndex( aVertexParentEntityIndices )
                , mVertexParentEntityRank( aVertexParentEntityRank )
        {
        }

        moris::Cell< moris_index >     mVertexParentEntityIndex;
        moris::Cell< mtk::EntityRank > mVertexParentEntityRank;

        moris_index
        get_vertex_parent_index( moris_index aVertexIndex ) const
        {
            return mVertexParentEntityIndex( aVertexIndex );
        }

        mtk::EntityRank
        get_vertex_parent_rank( moris_index aVertexIndex ) const
        {
            return mVertexParentEntityRank( aVertexIndex );
        }
    };

    // ----------------------------------------------------------------------------------

    struct Facet_Based_Connectivity
    {
        // in: index of facet || out: list of vertices (pointers) living on facet with the inputted index
        moris::Cell< moris::Cell< moris::mtk::Vertex* > > mFacetVertices;    // over allocated

        // in: index of facet || out: list of mtk::Cells (pointers) attached to facet with the inputted index
        moris::Cell< moris::Cell< moris::mtk::Cell* > > mFacetToCell;    // over allocated

        // in(1): index of facet; in(2): how many-eth mtk::Cell attached to facet || out: side ordinal (index) of this facet relative to this mtk::Cell
        moris::Cell< moris::Cell< moris_index > > mFacetToCellEdgeOrdinal;    // over allocated

        // in: index of mtk::Cell in List || out: list of facet-indices attached to it (facet indices as defined within this Facet_Based_Connectivity-object)
        moris::Cell< moris::Cell< moris_index > > mCellToFacet;    // over allocated

        // map relating Cell index in List to Cell index in Mesh
        std::unordered_map< moris_index, moris_index > mCellIndexToCellOrdinal;    // over allocated

        // in: index of vertex || out: list of facets (indices) connected to vertex with the inputted index
        moris::Cell< moris::Cell< moris_index > > mVertexFacets;

        moris_index
        get_cell_ordinal( const moris_index& aCellIndex )
        {
            auto tIter = mCellIndexToCellOrdinal.find( aCellIndex );
            MORIS_ASSERT( tIter != mCellIndexToCellOrdinal.end(), "Cell not in facet connectivity" );
            return tIter->second;
        }

        void
        merge_facets( moris_index aCellInd, moris_index aVInd1, moris::Cell< moris_index >& aFacetMergeInds, Facet_Based_Connectivity* aOldFacetConnectivity )
        {
            moris_index tCellOrdinal = this->get_cell_ordinal( aCellInd );

            moris::Cell< moris_index > tFacetInds = mCellToFacet( tCellOrdinal );
            tFacetInds.remove( MORIS_INDEX_MAX );


            moris_index FacetInd1 = tFacetInds( 1 );    // facet to delete
            moris_index FacetInd2 = tFacetInds( 0 );    // facet to keep
            for ( moris_index iV = 0; (uint)iV < mFacetVertices( tFacetInds( 0 ) ).size(); iV++ )
            {
                if ( mFacetVertices( tFacetInds( 0 ) )( iV )->get_index() == aVInd1 )
                {
                    FacetInd1 = tFacetInds( 0 );
                    FacetInd2 = tFacetInds( 1 );
                    break;
                }
            }


            if ( tFacetInds.size() > 2 )
            {
                std::cout << tFacetInds.size() << std::endl;
                MORIS_ERROR( false, "Error: cell being merged has more than 2 facets." );
            }

            aFacetMergeInds.push_back( FacetInd1 );


            // mVertexFacets
            // remove FacetInd1 from all mVertexFacets lists
            for ( uint iVert = 0; iVert < aOldFacetConnectivity->mFacetVertices( FacetInd1 ).size(); iVert++ )
            {
                moris_index tVertInd = aOldFacetConnectivity->mFacetVertices( FacetInd1 )( iVert )->get_index();
                if ( tVertInd != MORIS_INDEX_MAX )
                {
                    for ( uint iFacet = 0; iFacet < mVertexFacets( tVertInd ).size(); iFacet++ )
                    {
                        if ( mVertexFacets( tVertInd )( iFacet ) == FacetInd1 )
                        {
                            // mVertexFacets(tVertInd)(iFacet) = -1;
                            mVertexFacets( tVertInd ).erase( iFacet );
                            iFacet--;
                        }
                    }
                }
            }


            // mFacetVertices.erase(FacetInd1);
            if ( mFacetVertices( FacetInd1 ).size() == 2 )    // 2D
            {
                mFacetVertices( FacetInd1 ) = { NULL, NULL };
            }
            else if ( mFacetVertices( FacetInd1 ).size() == 3 )    // 3D
            {
                mFacetVertices( FacetInd1 ) = { NULL, NULL, NULL };
            }
            else
            {
                MORIS_ASSERT( false, "Error in mFacetVertices" );
            }


            // copy references to cells from facet being deleted to the one being merged
            for ( moris_index iC = 0; (uint)iC < mFacetToCell( FacetInd1 ).size(); iC++ )
            {
                if ( mFacetToCell( FacetInd1 )( iC ) != NULL )
                {
                    if ( mFacetToCell( FacetInd1 )( iC )->get_index() != aCellInd )
                    {

                        mFacetToCell( FacetInd2 ).push_back( mFacetToCell( FacetInd1 )( iC ) );    // todo: duplicates in 3d?
                        mFacetToCellEdgeOrdinal( FacetInd2 ).push_back( mFacetToCellEdgeOrdinal( FacetInd1 )( iC ) );
                    }
                }
            }


            // remove reference to deleted cell from mFacetToCell
            for ( moris_index iC = 0; (uint)iC < mFacetToCell( FacetInd2 ).size(); iC++ )
            {
                if ( mFacetToCell( FacetInd2 )( iC ) != NULL )
                {
                    if ( mFacetToCell( FacetInd2 )( iC )->get_index() == aCellInd )
                    {

                        mFacetToCell( FacetInd2 )( iC ) = NULL;

                        mFacetToCellEdgeOrdinal( FacetInd2 )( iC ) = MORIS_INDEX_MAX;
                    }
                }
            }

            // will eventually be deleted, but for merges necessary as reference
            mFacetToCell( FacetInd1 ) = mFacetToCell( FacetInd2 );


            // mFacetToCellEdgeOrdinal.erase(tFacetInds(0));
            mFacetToCellEdgeOrdinal( FacetInd1 ) = { MORIS_INDEX_MAX };


            // replace reference to facet
            for ( moris_index iC = 0; (uint)iC < mCellToFacet.size(); iC++ )
            {
                for ( moris_index iF = 0; (uint)iF < mCellToFacet( iC ).size(); iF++ )
                {
                    if ( mCellToFacet( iC )( iF ) == FacetInd1 )
                    {
                        mCellToFacet( iC )( iF ) = FacetInd2;
                    }
                }
            }
            // delete cell being merged
            // mCellToFacet.erase(this->get_cell_ordinal(aCellInd));
            mCellToFacet( tCellOrdinal ) = { MORIS_INDEX_MAX };

            mCellIndexToCellOrdinal.at( aCellInd ) = MORIS_INDEX_MAX;
        }

        moris::Cell< moris_index >
        verts_to_facets( moris_index aVertInd1, moris_index aVertInd2 )
        {

            moris::Cell< moris_index > tFacetList;
            tFacetList.reserve( 20 );    // estimate of how many facets connect two vertices. Always 2 for 2D, more for 3D

            // moris::Cell<moris_index> tFacetList2;
            moris::Cell< moris_index > tCommonFacets;
            tCommonFacets.resize( mVertexFacets( aVertInd1 ).size() + mVertexFacets( aVertInd2 ).size() );

            std::sort( mVertexFacets( aVertInd1 ).begin(), mVertexFacets( aVertInd1 ).end() );
            std::sort( mVertexFacets( aVertInd2 ).begin(), mVertexFacets( aVertInd2 ).end() );

            moris::Cell< moris::mtk::Cell* > tNullCell = { NULL };

            auto it = set_intersection( mVertexFacets( aVertInd1 ).begin(),
                    mVertexFacets( aVertInd1 ).end(),
                    mVertexFacets( aVertInd2 ).begin(),
                    mVertexFacets( aVertInd2 ).end(),
                    tCommonFacets.begin() );

            for ( auto st = tCommonFacets.begin(); st != it; ++st )
            {
                if ( mFacetVertices( *st )( 0 ) != NULL )
                {
                    tFacetList.push_back( *st );
                }
            }


            if ( tFacetList.size() == 0 )
            {
                MORIS_ERROR( false, "Error: No facet found for two vertices." );
            }

            return tFacetList;
        }

        void
        merge_vertices( moris_index aVInd1, moris_index aVInd2, moris::Cell< moris_index >& aFacetMergeInds, Facet_Based_Connectivity* aOldFacetConnectivity, mtk::Vertex* aVertex2 )
        {
            // facet indices being merged
            moris::Cell< moris_index > aFacetIndices = this->verts_to_facets( aVInd1, aVInd2 );


            // fix mVertexFacets ----------------

            // add facets from aVert1 to aVert2
            for ( uint iFacet = 0; iFacet < mVertexFacets( aVInd1 ).size(); iFacet++ )
            {
                if ( std::find( aFacetIndices.begin(), aFacetIndices.end(), mVertexFacets( aVInd1 )( iFacet ) ) == aFacetIndices.end() )
                {
                    mVertexFacets( aVInd2 ).push_back( mVertexFacets( aVInd1 )( iFacet ) );
                }
            }

            // set merged facets to -1
            for ( uint iDelFacets = 0; iDelFacets < aFacetIndices.size(); iDelFacets++ )
            {
                for ( uint iVert = 0; iVert < mFacetVertices( iDelFacets ).size(); iVert++ )
                {
                    if ( mFacetVertices( iDelFacets )( iVert ) != NULL )
                    {
                        moris_index tVertInd = mFacetVertices( iDelFacets )( iVert )->get_index();
                        if ( tVertInd != MORIS_INDEX_MAX )
                        {
                            for ( uint iFacet = 0; iFacet < mVertexFacets( tVertInd ).size(); iFacet++ )
                            {
                                if ( mVertexFacets( tVertInd )( iFacet ) == aFacetIndices( iDelFacets ) )
                                {
                                    // mVertexFacets(tVertInd)(iFacet) = -1;
                                    mVertexFacets( tVertInd ).erase( iFacet );
                                    iFacet--;
                                }
                            }
                        }
                    }
                }
            }

            mVertexFacets( aVInd1 ) = {};


            // replace pointer to aVInd1 with pointer to aVInd2 in mFacetVertices
            for ( uint iFacet = 0; iFacet < aOldFacetConnectivity->mVertexFacets( aVInd1 ).size(); iFacet++ )
            {
                moris_index tFacetInd = aOldFacetConnectivity->mVertexFacets( aVInd1 )( iFacet );
                for ( uint iVert = 0; iVert < mFacetVertices( tFacetInd ).size(); iVert++ )
                {
                    if ( mFacetVertices( tFacetInd )( 0 ) != NULL )    // todo: check mFacetVertices(iF) != {NULL, NULL}
                    {
                        moris_index tVertInd = mFacetVertices( tFacetInd )( iVert )->get_index();
                        if ( tVertInd == aVInd1 )
                        {
                            mFacetVertices( tFacetInd )( iVert ) = aVertex2;
                        }
                    }
                }
            }


            // remove facets ------------------------------------------------------------

            for ( uint iFacet = 0; iFacet < aFacetIndices.size(); iFacet++ )    // go through facets to delete
            {
                aFacetMergeInds.push_back( aFacetIndices( iFacet ) );

                // go through cells attached to facet
                for ( uint iCell = 0; iCell < aOldFacetConnectivity->mFacetToCell( aFacetIndices( iFacet ) ).size(); iCell++ )
                {
                    if ( aOldFacetConnectivity->mFacetToCell( aFacetIndices( iFacet ) )( iCell ) != NULL )
                    {
                        moris_index tCellInd = mFacetToCell( aFacetIndices( iFacet ) )( iCell )->get_index();

                        if ( tCellInd != MORIS_INDEX_MAX )
                        {
                            moris_index tCellOrd = this->get_cell_ordinal( tCellInd );
                            for ( uint iCellFacet = 0; iCellFacet < mCellToFacet( tCellOrd ).size(); iCellFacet++ )
                            {
                                if ( mCellToFacet( tCellOrd )( iCellFacet ) == aFacetIndices( iFacet ) )
                                {
                                    // mCellToFacet(get_cell_ordinal(tCellInd)).erase(iFf);
                                    mCellToFacet( tCellOrd )( iCellFacet ) = MORIS_INDEX_MAX;
                                }
                            }
                        }
                    }
                }

                if ( mFacetVertices( aFacetIndices( iFacet ) ).size() == 2 )    // 2D
                {
                    mFacetVertices( aFacetIndices( iFacet ) ) = { NULL, NULL };
                }
                else if ( mFacetVertices( aFacetIndices( iFacet ) ).size() == 3 )    // 3D
                {
                    mFacetVertices( aFacetIndices( iFacet ) ) = { NULL, NULL, NULL };
                }

                // mFacetToCell.erase(aFacetIndices(iFacet));
                mFacetToCell( aFacetIndices( iFacet ) ) = { NULL };

                // mFacetToCellEdgeOrdinal.erase(aFacetIndices(iFacet));
                mFacetToCellEdgeOrdinal( aFacetIndices( iFacet ) ) = { MORIS_INDEX_MAX };
            }
        }

    };    // struct Facet_Based_Connectivity

    // ----------------------------------------------------------------------------------

    struct Facet_Based_Ancestry
    {
        moris::Cell< moris_index > mFacetParentEntityIndex;
        moris::Cell< moris_index > mFacetParentEntityRank;
        moris::Cell< moris_index > mFacetParentEntityOrdinalWrtBackgroundCell;
    };

    // ----------------------------------------------------------------------------------

    struct Cell_Neighborhood_Connectivity
    {
        // Cells of Cells (i.e. matrix-lists) for Cell-Connectivity information
        // first index is mtk::Cell for which the neighborhood is to be defined
        // second index is List of mtk::Cells connected to mtk::Cell with first index (connection through a facet)

        // pointers to connected mtk::Cells
        moris::Cell< std::shared_ptr< moris::Cell< moris::mtk::Cell* > > > mNeighborCells;

        // indices of side ordinals through which the mtk::Cell of first index connects to the mtk::Cells of second indices
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mMySideOrdinal;

        // fixme: this can be deleted, as it is not used ?!
        // indices of side ordinals through which the mtk::Cells of second indices connects to the mtk::Cells of the first index
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mNeighborSideOrdinal;
    };

    // ----------------------------------------------------------------------------------

    struct Cell_Connectivity
    {
        Cell_Connectivity(){};
        Cell_Connectivity(
                Matrix< IndexMat > const & aCellVertexInds,
                Matrix< IndexMat > const & aCellEdgesInds,
                Matrix< IndexMat > const & aCellFacesInds )
                : mCellVertexInds( aCellVertexInds )
                , mCellEdgesInds( aCellEdgesInds )
                , mCellFacesInds( aCellFacesInds )
        {
        }

        moris_index
        get_entity_index(
                moris_index     aEntityOrdinal,
                mtk::EntityRank aEntityRank ) const
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::FACE:
                    return mCellFacesInds( aEntityOrdinal );
                    break;
                case mtk::EntityRank::EDGE:
                    return mCellEdgesInds( aEntityOrdinal );
                    break;
                case mtk::EntityRank::NODE:
                    return mCellVertexInds( aEntityOrdinal );
                    break;
                default:
                    MORIS_ERROR( 0, "UNDEFINED RANK" );
                    return 0;
                    break;
            }
        }

        moris_index
        get_entity_ordinal(
                moris_index const & aEntityIndex,
                mtk::EntityRank     aEntityRank ) const
        {
            switch ( aEntityRank )
            {
                case mtk::EntityRank::FACE:
                    for ( uint iEnt = 0; iEnt < mCellFacesInds.numel(); iEnt++ )
                    {
                        if ( mCellFacesInds( iEnt ) == aEntityIndex )
                        {
                            return iEnt;
                        }
                    }
                    MORIS_ERROR( 0, "Face not found" );
                    return MORIS_INDEX_MAX;
                    break;
                case mtk::EntityRank::EDGE:
                    for ( uint iEnt = 0; iEnt < mCellEdgesInds.numel(); iEnt++ )
                    {
                        if ( mCellEdgesInds( iEnt ) == aEntityIndex )
                        {
                            return iEnt;
                        }
                    }
                    MORIS_ERROR( 0, "Edge not found" );
                    return MORIS_INDEX_MAX;
                    break;
                case mtk::EntityRank::NODE:
                    for ( uint iEnt = 0; iEnt < mCellVertexInds.numel(); iEnt++ )
                    {
                        if ( mCellVertexInds( iEnt ) == aEntityIndex )
                        {
                            return iEnt;
                        }
                    }
                    MORIS_ERROR( 0, "Node not found" );
                    return MORIS_INDEX_MAX;
                    break;
                default:
                    MORIS_ERROR( 0, "UNDEFINED RANK" );
                    return 0;
                    break;
            }
        }

        const Matrix< IndexMat > mCellVertexInds;
        const Matrix< IndexMat > mCellEdgesInds;
        const Matrix< IndexMat > mCellFacesInds;
    };

    // ----------------------------------------------------------------------------------

    struct Subphase_Neighborhood_Connectivity
    {
        // input: sub-phase index || output: list of sub-phases connected to it
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhase;

        // input: sub-phase index || output: list of facet ordinals belonging to parent cell through which the sub-phases are connected
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhaseMySideOrds;

        // input: sub-phase index || output: list of facet ordinals belonging to parent cell of the neighboring sub-phase through which the sub-phases are connected
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mSubphaseToSubPhaseNeighborSideOrds;

        // input: sub-phase index || output: // TODO: some info needed when having a refinement boundary
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mTransitionNeighborCellLocation;

        // auxiliary array to perform depth first search,  it is empty by default to save memory
        // input: sub-phase index || output: flag indicating if the vertex is visited
        moris::Cell< bool > mVisitedFlag;

        void
        print_subphase_neighborhood()
        {

            std::cout << "Subphases" << std::endl;
            for ( uint iC = 0; iC < mSubphaseToSubPhase.size(); iC++ )
            {
                std::cout << std::setw( 6 ) << iC << " | ";

                for ( uint iN = 0; iN < mSubphaseToSubPhase( iC )->size(); iN++ )
                {
                    std::cout << std::setw( 6 ) << ( *mSubphaseToSubPhase( iC ) )( iN );
                }
                std::cout << std::endl;
            }

            std::cout << "Subphases My Side Ordinals" << std::endl;
            for ( uint iC = 0; iC < mSubphaseToSubPhaseMySideOrds.size(); iC++ )
            {
                std::cout << std::setw( 6 ) << iC << " | ";

                for ( uint iN = 0; iN < mSubphaseToSubPhaseMySideOrds( iC )->size(); iN++ )
                {
                    std::cout << std::setw( 6 ) << ( *mSubphaseToSubPhaseMySideOrds( iC ) )( iN );
                }
                std::cout << std::endl;
            }

            std::cout << "Subphases Neighbor Side Ordinals" << std::endl;
            for ( uint iC = 0; iC < mSubphaseToSubPhaseNeighborSideOrds.size(); iC++ )
            {
                std::cout << std::setw( 6 ) << iC << " | ";

                for ( uint iN = 0; iN < mSubphaseToSubPhaseNeighborSideOrds( iC )->size(); iN++ )
                {
                    std::cout << std::setw( 6 ) << ( *mSubphaseToSubPhaseNeighborSideOrds( iC ) )( iN );
                }
                std::cout << std::endl;
            }

            std::cout << "Transition Neighbor Locations" << std::endl;
            for ( uint iC = 0; iC < mTransitionNeighborCellLocation.size(); iC++ )
            {
                std::cout << std::setw( 6 ) << iC << " | ";

                for ( uint iN = 0; iN < mTransitionNeighborCellLocation( iC )->size(); iN++ )
                {
                    std::cout << std::setw( 12 ) << ( *mTransitionNeighborCellLocation( iC ) )( iN );
                }
                std::cout << std::endl;
            }
        }

        // ----------------------------------------------------------------------------------

        void
        depth_first_search( moris_index const & aSubphaseIndex, moris_index const & aDegree, moris::Cell< moris_index >& aNeighbors )
        {
            // set the starting vertex as visited
            mVisitedFlag( aSubphaseIndex ) = true;

            // break out of recursive function if it reached to vertex itself
            if ( 0 == aDegree )
            {
                aNeighbors.push_back( aSubphaseIndex );
                return;
            }

            // loop ove the neighbour and if not visited then perform depth first search
            for ( const auto& iNeighbour : *( mSubphaseToSubPhase( aSubphaseIndex ) ) )
            {
                if ( !mVisitedFlag( iNeighbour ) )
                {
                    depth_first_search( iNeighbour, aDegree - 1, aNeighbors );
                }
            }
        }

        // ----------------------------------------------------------------------------------

        void
        get_kth_degree_neighbours( moris_index const & aSubphaseIndex, moris_index const & aDegree, moris::Cell< moris_index >& aNeighbors )
        {
            // mark all the vertices as not visited
            mVisitedFlag.resize( mSubphaseToSubPhase.size());

            // reset all vertices to false for a new traversal
            mVisitedFlag.assign(mSubphaseToSubPhase.size(), false );

            // perform depth first search
            depth_first_search( aSubphaseIndex, aDegree, aNeighbors );
        }
    };

    // ----------------------------------------------------------------------------------

    class Child_Mesh_Experimental;
    class Model;
    class Cell_XTK_No_CM;
    class Cut_Integration_Mesh : public moris::mtk::Mesh
    {
        friend class Integration_Mesh_Generator;
        friend class Integration_Mesh_Cleanup;
        friend class Model;

      protected:
        uint mSpatialDimension;

        bool mSameLevelChildMeshes = true;

        // integration cells
        moris_index                                           mFirstControlledCellIndex;
        moris::Cell< moris::mtk::Cell* >                      mIntegrationCells;
        moris::Cell< std::shared_ptr< xtk::Cell_XTK_No_CM > > mControlledIgCells;

        // quantities related to integration cells
        moris::Cell< moris::Cell< moris_index > > mIntegrationCellToCellGroupIndex;
        moris::Cell< moris_index >                mIntegrationCellToSubphaseIndex;
        moris::Cell< moris_index >                mIntegrationCellBulkPhase;

        // integration vertices
        moris_index                                              mFirstControlledVertexIndex;
        moris::Cell< moris::mtk::Vertex* >                       mIntegrationVertices;
        moris::Cell< std::shared_ptr< moris::mtk::Vertex_XTK > > mControlledIgVerts;

        // vertex ancestry
        moris::Cell< moris_index > mIgVertexParentEntityIndex;
        moris::Cell< moris_index > mIgVertexParentEntityRank;
        moris::Cell< moris_index > mIgVertexConnectedCell;

        // vertex quantities
        moris::Cell< std::shared_ptr< Matrix< DDRMat > > > mVertexCoordinates;

        // all data is stored in the current mesh. pointers are in the child mesh
        // as well as accessor functions are provided there
        moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mChildMeshes;
        moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mOwnedChildMeshes;
        moris::Cell< std::shared_ptr< Child_Mesh_Experimental > > mNotOwnedChildMeshes;

        // communication map
        moris::Matrix< IdMat >            mCommunicationTable;
        std::map< moris_id, moris_index > mCommunicationMap;
        bool                              mCommMapHasBeenConstructed = false;

        // Integration - Lagrange Mesh relation
        // group of all integration cells in a single parent cell
        moris::Cell< std::shared_ptr< IG_Cell_Group > >   mIntegrationCellGroups;
        moris::Cell< std::shared_ptr< IG_Vertex_Group > > mIntegrationVertexGroups;
        moris::Cell< moris::mtk::Cell* >                  mIntegrationCellGroupsParentCell;
        moris::Cell< moris_index >                        mParentCellCellGroupIndex;

        moris::Cell< moris_index > mOwnedIntegrationCellGroupsInds;
        moris::Cell< moris_index > mNotOwnedIntegrationCellGroups;

        // Lagrange Mesh B-Spline Mesh relation
        moris::Cell< Bspline_Mesh_Info* > mBsplineMeshInfos;

        // B-spline mesh index with the coarsest B-spline element containing a given base IP cell
        moris::Cell< moris_index > mCoarsestBsplineMesh;    // input: base IP cell index || output: coarsest B-spline mesh

        // union multisets of void MSD indices for each Lagrange element
        // input: Lagrange/ base IP cell index || output: Union multiset of MSD indices for that base IP cell
        moris::Cell< moris::Cell< moris_index > > mUnionVoidMsdIndices;

        // bulk-phases the void material sub-domains belong to
        // input: Lagrange/ base IP cell index || output: bulk phase indices corresponding to union void MSD indices
        moris::Cell< moris::Cell< moris_index > > mUnionVoidMsdIndexBulkPhases;

        // subphase groupings
        moris::Cell< moris_index >                      mSubPhaseIds;              // input: sub-phase index || output: global sub-phase ID
        moris::Cell< std::shared_ptr< IG_Cell_Group > > mSubPhaseCellGroups;       // input: sub-phase index || output: pointer to IG-Cell group on which subphase lives
        moris::Cell< moris_index >                      mSubPhaseBulkPhase;        // input: sub-phase index || output: index of bulk-phase (i.e. material phase)
        moris::Cell< moris::mtk::Cell* >                mSubPhaseParentCell;       // input: sub-phase index || output: index of Bg-cell sub-phase lives on
        moris::Cell< moris::Cell< moris_index > >       mParentCellToSubphase;     // input: Bg-cell index   || output: list of sub-phase indices present in Bg-cell
        moris::Cell< moris_index >                      mParentCellHasChildren;    // input: Bg-cell index   || output: bool, whether Bg-cell has Ig-cells living on it

        moris::Cell< moris_index >                  mOwnedSubphaseGroupsInds;
        moris::Cell< moris_index >                  mNotOwnedSubphaseGroupsInds;
        std::unordered_map< moris_id, moris_index > mGlobalToLocalSubphaseMap;

        // subphase connectivity
        std::shared_ptr< Subphase_Neighborhood_Connectivity > mSubphaseNeighborhood;

        // subphase-group connectivity
        moris::Cell< std::shared_ptr< Subphase_Neighborhood_Connectivity > > mSubphaseGroupNeighborhood;

        // face connectivity
        std::shared_ptr< Facet_Based_Connectivity > mIgCellFaceConnectivity;

        // face ancestry
        std::shared_ptr< Facet_Based_Ancestry > mIgCellFaceAncestry;

        // interface facets - indexed based on mIgCellFaceConnectivity facet indices
        moris::Cell< moris_index > mInterfaceFacets;

        // double side interface groups
        // outer cell - bulk phase 0
        // inner cell - bulk phase 1
        // IG_Cell_Double_Side_Group pairings between integration cells
        moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > mBpToBpDblSideInterfaces;

        // background facet to child facet connectivity
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > mBGFacetToChildFacet;

        // block set data
        std::unordered_map< std::string, moris_index >  mBlockSetLabelToOrd;
        moris::Cell< std::string >                      mBlockSetNames;
        moris::Cell< std::shared_ptr< IG_Cell_Group > > mBlockSetCellGroup;
        moris::Cell< mtk::CellTopology >                mBlockCellTopo;

        // Side Set Data
        std::unordered_map< std::string, moris_index >       mSideSideSetLabelToOrd;
        Cell< std::string >                                  mSideSetLabels;
        moris::Cell< std::shared_ptr< IG_Cell_Side_Group > > mSideSetCellSides;

        // connectivity from vertex to child mesh
        // outer cell vertex
        // inner cell child meshes associated with the vertex
        moris::Cell< moris::Cell< moris_index > > mVertexToChildMeshIndex;

        // connectivity from vertex to child mesh
        // outer cell integration cell index
        // inner cell child meshes associated with the integration cell
        moris::Cell< moris::Cell< moris_index > > mCellToChildMeshIndex;

        // outer cell geometry index
        // if vertex index is in the map then it is a member of the geometric interface
        moris::Cell< std::unordered_map< moris_index, moris_index > > mGeometryInterfaceVertexIndices;

        moris_index mGlobalMaxVertexId;
        moris_index mGlobalMaxCellId;

        std::unordered_map< moris_id, moris_index > mIntegrationCellIdToIndexMap;
        std::unordered_map< moris_id, moris_index > mIntegrationVertexIdToIndexMap;

        moris::Cell< moris_index > mIntegrationCellIndexToId;
        moris::Cell< moris_index > mIntegrationVertexIndexToId;

        moris::mtk::Mesh* mBackgroundMesh;
        Model*            mXTKModel;

      public:
        // ----------------------------------------------------------------------------------

        Cut_Integration_Mesh(
                moris::mtk::Mesh* aBackgroundMesh,
                Model*            aXTKModel );

        // ----------------------------------------------------------------------------------

        ~Cut_Integration_Mesh();

        // ----------------------------------------------------------------------------------

        void delete_Bspline_mesh_info();

        // ----------------------------------------------------------------------------------

        // Core Mesh Functions
        uint
        get_spatial_dim() const;

        // ----------------------------------------------------------------------------------

        mtk::MeshType
        get_mesh_type() const;

        // ----------------------------------------------------------------------------------

        uint
        get_num_entities(
                mtk::EntityRank   aEntityRank,
                const moris_index aIndex ) const;

        // ----------------------------------------------------------------------------------

        uint
        get_num_base_ip_cells() const;

        // ----------------------------------------------------------------------------------

        uint get_num_sets() const;

        // ----------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_node_coordinate( moris_index aNodeIndex ) const;

        // ----------------------------------------------------------------------------------

        uint
        get_node_owner( moris_index aNodeIndex ) const;

        // ----------------------------------------------------------------------------------

        uint
        get_element_owner( moris_index aElementIndex ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IdMat >
        get_communication_table() const;

        // ----------------------------------------------------------------------------------

        std::map< moris_id, moris_index >
        get_communication_map();

        // ----------------------------------------------------------------------------------

        void
        add_proc_to_comm_table( moris_index aProcRank );

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex );

        // ----------------------------------------------------------------------------------

        mtk::CellTopology
        get_blockset_topology( const std::string& aSetName );

        // ----------------------------------------------------------------------------------

        mtk::CellShape
        get_IG_blockset_shape( const std::string& aSetName );

        // ----------------------------------------------------------------------------------

        mtk::CellShape
        get_IP_blockset_shape( const std::string& aSetName );

        // ----------------------------------------------------------------------------------

        moris_id
        get_glb_entity_id_from_entity_loc_index(
                moris_index     aEntityIndex,
                mtk::EntityRank aEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;

        // ----------------------------------------------------------------------------------

        moris_index
        get_loc_entity_ind_from_entity_glb_id(
                moris_id        aEntityId,
                mtk::EntityRank aEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_entity_connected_to_entity_loc_inds(
                moris_index       aEntityIndex,
                mtk::EntityRank   aInputEntityRank,
                mtk::EntityRank   aOutputEntityRank,
                const moris_index aDiscretizationIndex = 0 ) const;

        // ----------------------------------------------------------------------------------

        moris::Cell< std::string >
        get_set_names( mtk::EntityRank aSetEntityRank ) const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_block_set_index( std::string aBlockSetLabel ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_block_entity_loc_inds( std::string aSetName ) const;

        // ----------------------------------------------------------------------------------

        moris_index
        get_side_set_index( std::string aSideSetLabel ) const;

        // ----------------------------------------------------------------------------------

        void
        get_sideset_elems_loc_inds_and_ords(
                const std::string&  aSetName,
                Matrix< IndexMat >& aElemIndices,
                Matrix< IndexMat >& aSidesetOrdinals ) const;

        // ----------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_set_entity_loc_inds(
                mtk::EntityRank aSetEntityRank,
                std::string     aSetName ) const;

        // ----------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Matrix< DDRMat > > >*
        get_all_vertex_coordinates_loc_inds();

        // ----------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< IG_Cell_Group > >&
        get_all_cell_groups();

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Child_Mesh_Experimental >
        get_child_mesh( moris_index aChildMeshIndex );

        // ----------------------------------------------------------------------------------

        uint
        get_num_child_meshes() const;

        // ----------------------------------------------------------------------------------

        std::unordered_map< moris_id, moris_index >
        get_vertex_glb_id_to_loc_vertex_ind_map() const;

        // ----------------------------------------------------------------------------------

        bool
        vertex_exists( moris_index tId ) const;

        // ----------------------------------------------------------------------------------

        mtk::Cell const &
        get_mtk_cell( moris_index aElementIndex ) const;

        // ----------------------------------------------------------------------------------

        mtk::Cell&
        get_mtk_cell( moris_index aElementIndex );

        // ----------------------------------------------------------------------------------

        mtk::Vertex&
        get_mtk_vertex( moris_index aVertexIndex );

        // ----------------------------------------------------------------------------------

        mtk::Vertex const &
        get_mtk_vertex( moris_index aVertexIndex ) const;

        // ----------------------------------------------------------------------------------

        moris::mtk::Vertex*
        get_mtk_vertex_pointer( moris_index aVertexIndex );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< IG_Vertex_Group >
        get_vertex_group( moris_index aVertexGroupIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_parent_cell_group_index( moris_index aParentCellIndex );

        // ----------------------------------------------------------------------------------

        void
        replace_controlled_ig_cell(
                moris_index                              aCellIndex,
                moris_id                                 aCellId,
                std::shared_ptr< moris::mtk::Cell_Info > aCellInfo,
                moris::Cell< moris::mtk::Vertex* >&      aVertexPointers );

        // ----------------------------------------------------------------------------------

        void
        set_integration_cell(
                moris_index                            aCellIndex,
                std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell );

        // ----------------------------------------------------------------------------------

        void
        add_integration_cell(
                moris_index                            aCellIndex,
                std::shared_ptr< xtk::Cell_XTK_No_CM > aNewCell );

        // ----------------------------------------------------------------------------------

        moris_index
        get_integration_cell_controlled_index(
                moris_index aCellIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_integration_vertex_controlled_index(
                moris_index aVertIndex );

        // ----------------------------------------------------------------------------------

        void
        add_cell_to_cell_group(
                moris_index aCellIndex,
                moris_index aCellGroupIndex );

        // ----------------------------------------------------------------------------------

        moris_id
        allocate_entity_ids( moris::size_t aNumIdsToAllocate,
                mtk::EntityRank            aEntityRank );

        // ----------------------------------------------------------------------------------

        /**
         * @brief allocate new IDs for Subphases globally across all procs
         *
         * @param aNumIdsToAllocate number of Subphase IDs the current processor would like to allocate
         * @return moris_id first ID of the range of IDs that gets assigned to the current processor's Subphases
         */
        moris_id
        allocate_subphase_ids( moris::size_t aNumIdsToAllocate );

        // ----------------------------------------------------------------------------------

        moris_index
        get_first_available_index( mtk::EntityRank aEntityRank ) const;

        // ----------------------------------------------------------------------------------

        uint
        get_num_ig_cell_groups();

        // ----------------------------------------------------------------------------------

        std::shared_ptr< IG_Cell_Group >
        get_ig_cell_group( moris_index aGroupIndex );

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index > const &
        get_ig_cell_group_memberships( moris_index aIgCellIndex );

        // ----------------------------------------------------------------------------------

        moris::mtk::Cell*
        get_ig_cell_group_parent_cell( moris_index aGroupIndex );

        // ----------------------------------------------------------------------------------

        mtk::CellTopology
        get_child_element_topology();
        void

        // ----------------------------------------------------------------------------------

        set_child_mesh_subphase(
                moris_index                 aCMIndex,
                moris::Cell< moris_index >& aSubphasesGroups );

        // ----------------------------------------------------------------------------------

        uint
        get_num_subphases();

        // ----------------------------------------------------------------------------------

        uint
        get_num_subphase_groups( moris_index aMeshListIndex );

        // ----------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Child_Mesh_Experimental > >&
        get_owned_child_meshes();

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index >&
        get_owned_subphase_indices();

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index >&
        get_not_owned_subphase_indices();

        // ----------------------------------------------------------------------------------

        moris::mtk::Cell*
        get_subphase_parent_cell( moris_index aSubPhaseIndex );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< IG_Cell_Group >
        get_subphase_ig_cells( moris_index aSubPhaseIndex );

        // ----------------------------------------------------------------------------------

        const moris::Cell< moris_index >&
        get_ig_cells_in_SPG(
                moris_index aMeshIndexInList,
                moris_index aSubphaseGroupIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_id( moris_index aSubphaseIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_index( moris_id aSubphaseId );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_bulk_phase( moris_index aSubphaseIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_group_id(
                moris_index aSubphaseGroupIndex,
                moris_index aBsplineMeshListIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_group_index(
                moris_id    aSubphaseGroupId,
                moris_index aBsplineMeshListIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_subphase_group_bulk_phase(
                moris_index aSubphaseGroupIndex,
                moris_index aBsplineMeshListIndex );

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index > const &
        get_parent_cell_subphases( moris_index aParentCellIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_ig_cell_subphase_index( moris_index aIgCellIndex );

        // ----------------------------------------------------------------------------------

        bool
        parent_cell_has_children( moris_index aParentCellIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_vertex_parent_index( moris_index const & aVertexIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_vertex_parent_rank( moris_index const & aVertexIndex );

        // ----------------------------------------------------------------------------------

        void
        finalize_cut_mesh_construction();

        // ----------------------------------------------------------------------------------

        void
        deduce_ig_cell_group_ownership();

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief Assigns IDs to all IG cells (owned and not-owned). The IDs from not owned entities are obtained from other processors.
         */
        void
        assign_controlled_ig_cell_ids();

        // ----------------------------------------------------------------------------------

        /**
         * @brief assigns IDs to the IG cells owned by the executing processor
         *
         * @param aFirstFreeId first free IG cell ID to be used by the executing processor
         */
        void
        assign_IDs_to_owned_IG_cells( moris_id aFirstFreeId );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Prepares lists of IG cells the executing proc has access to but doesn't own for each processor that actually owns these subphases.
         *
         * @param aNotOwnedIgCellGroups outer cell: owner proc index in XTK comm-table || inner cell: list of non-owned IG cell groups whose IG cell's IDs need to be requested from that owning proc
         * @param aParentCellIds outer cell: owner proc index in XTK comm-table || inner cell: parent Cell IDs corresponding to the IG cell groups in the lists above
         * @param aNumIgCellsInParentCell outer cell: owner proc index in XTK comm-table || inner cell: number of IG cells in each of the IG cell groups in the lists above
         */
        void
        prepare_requests_for_not_owned_IG_cell_IDs(
                Cell< Cell< moris_index > >& aNotOwnedIgCellGroups,
                Cell< Matrix< IdMat > >&     aParentCellIds,
                Cell< Matrix< IndexMat > >&  aNumIgCellsInParentCell );

        // ----------------------------------------------------------------------------------

        /**
         * @brief find and answer with the IDs of the IG cells owned by the executing processor and requested from other processors
         *
         * @param aFirstIgCellIdsInCellGroups output: outer cell: owner proc index in XTK comm-table || inner cell: IDs of the first IG cell in the respective IG cell groups specified below
         * @param aReceivedParentCellIds input: outer cell: owner proc index in XTK comm-table || inner cell: parent Cell IDs corresponding to the IG cell groups whose IDs are requested
         * @param aReceivedNumIgCellsInParentCell input: outer cell: owner proc index in XTK comm-table || inner cell: number of IG cells in each of the IG cell groups whose IG cell's IDs are requested
         */
        void
        prepare_answers_for_owned_IG_cell_IDs(
                Cell< Matrix< IdMat > >&           aFirstIgCellIdsInCellGroups,
                Cell< Matrix< IdMat > > const &    aReceivedParentCellIds,
                Cell< Matrix< IndexMat > > const & aReceivedNumIgCellsInParentCell );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Assign the received IG cell IDs to on the executing processor
         *
         * @param aNotOwnedIgCellGroups // outer cell: owner proc index in XTK comm-table || inner cell: list of not owned IG cell groups whose IG cell's IDs was requested from that owning proc
         * @param aReceivedFirstIgCellIdsInCellGroups // outer cell: owner proc index in XTK comm-table || inner cell: IDs of the not owned first IG cell within the group specified
         */
        void
        handle_requested_IG_cell_ID_answers(
                Cell< Cell< moris_index > > const & aNotOwnedIgCellGroups,
                Cell< Matrix< IdMat > > const &     aReceivedFirstIgCellIdsInCellGroups );

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        void
        set_face_connectivity( std::shared_ptr< Facet_Based_Connectivity > aFaceConnectivity );

        // ----------------------------------------------------------------------------------

        void
        set_face_ancestry( std::shared_ptr< Facet_Based_Ancestry > aFaceAncestry );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Facet_Based_Connectivity >
        get_face_connectivity();

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Facet_Based_Ancestry >
        get_face_ancestry();

        // ----------------------------------------------------------------------------------

        void
        set_interface_facets( moris::Cell< moris_index >& aInterfaces );

        // ----------------------------------------------------------------------------------

        moris::Cell< moris_index > const &
        get_interface_facets();

        // ----------------------------------------------------------------------------------

        void
        set_bulk_phase_to_bulk_phase_dbl_side_interface( moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > >& aBpToBpDblSideInterfaces );

        // ----------------------------------------------------------------------------------

        moris::Cell< moris::Cell< std::shared_ptr< IG_Cell_Double_Side_Group > > > const &
        get_bulk_phase_to_bulk_phase_dbl_side_interface();

        // ----------------------------------------------------------------------------------

        void
        set_background_facet_to_child_facet_connectivity( moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & aBgToChildFacet );

        // ----------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const &
        get_background_facet_to_child_facet_connectivity();

        // ----------------------------------------------------------------------------------

        void
        set_subphase_neighborhood( std::shared_ptr< Subphase_Neighborhood_Connectivity > aSubphaseNeighborhood );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Subphase_Neighborhood_Connectivity >
        get_subphase_neighborhood();

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Subphase_Neighborhood_Connectivity >
        get_subphase_group_neighborhood( moris_index aMeshIndex );

        // ----------------------------------------------------------------------------------
        moris::Cell< Bspline_Mesh_Info* >&
        get_bspline_mesh_info();

        // ----------------------------------------------------------------------------------

        void
        setup_glob_to_loc_subphase_map();

        // ----------------------------------------------------------------------------------

        void
        construct_spg_id_to_index_map( moris_index aBsplineMeshListIndex );

        // ----------------------------------------------------------------------------------

        moris_index
        get_cell_bulk_phase( moris_index aCellIndex );

        // ----------------------------------------------------------------------------------

        Cell< moris_index >
        register_side_set_names( moris::Cell< std::string > const & aSideSetNames );

        // ----------------------------------------------------------------------------------

        Cell< moris_index >
        register_block_set_names( moris::Cell< std::string > const & aBlockSetNames,
                mtk::CellTopology                                    aCellTopo );

        // ----------------------------------------------------------------------------------

        void
        write_mesh( std::string aOutputPath,
                std::string     aOutputFile );

        // ----------------------------------------------------------------------------------

        Cell_Connectivity
        get_background_cell_connectivity( moris_index aBGCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the B-spline mesh index with the coarsest B-spline element
         * which contains the requested base IP cell
         *
         * @param aBaseIpCellIndex
         * @return moris_index
         */
        moris_index
        get_coarsest_bspline_mesh_index_on_base_ip_cell( moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the union MSD indices for base IP cell
         *
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& list of Union MSD indices on this IP cell
         */
        moris::Cell< moris_index > const &
        get_union_MSD_indices_for_base_IP_cell( const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the bulk-phases associated with the union MSDs for a given base IP cell
         *
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& list of bulk phases corresponding to the union MSDs on this IP cell
         */
        moris::Cell< moris_index > const &
        get_bulk_phases_for_union_MSD_indices_for_base_IP_cell( const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the material SPG indices for base IP cell
         *
         * @param aBsplineMeshListIndex index of the B-spline mesh in the list of B-spline meshes to be enriched
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& SPGs wrt. which material clusters will be constructed
         */
        moris::Cell< moris_index > const &
        get_material_SPG_indices_for_base_IP_cell(
                const moris_index aBsplineMeshListIndex,
                const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the material MSD indices for base IP cell
         *
         * @param aBsplineMeshListIndex index of the B-spline mesh in the list of B-spline meshes to be enriched
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& MSD indices corresponding to SPGs wrt. which material clusters will be constructed
         */
        moris::Cell< moris_index > const &
        get_material_MSD_indices_for_base_IP_cell(
                const moris_index aBsplineMeshListIndex,
                const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the void SPG indices for base IP cell
         *
         * @param aBsplineMeshListIndex index of the B-spline mesh in the list of B-spline meshes to be enriched
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& SPGs wrt. which void clusters need to be constructed
         */
        moris::Cell< moris_index > const &
        get_void_SPG_indices_for_base_IP_cell(
                const moris_index aBsplineMeshListIndex,
                const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the void MSD indices for base IP cell
         *
         * @param aBsplineMeshListIndex index of the B-spline mesh in the list of B-spline meshes to be enriched
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& MSD indices corresponding to SPGs wrt. which void clusters need to be constructed
         */
        moris::Cell< moris_index > const &
        get_void_MSD_indices_for_base_IP_cell(
                const moris_index aBsplineMeshListIndex,
                const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the free void MSD indices for base IP cell
         *
         * @param aBsplineMeshListIndex index of the B-spline mesh in the list of B-spline meshes to be enriched
         * @param aBaseIpCellIndex index of the base IP cell / Lagrange element
         * @return moris::Cell< moris_index > const& MSD indices without associated SPGs for which void clusters need to be constructed
         */
        moris::Cell< moris_index > const &
        get_free_void_MSD_indices_for_base_IP_cell(
                const moris_index aBsplineMeshListIndex,
                const moris_index aBaseIpCellIndex ) const;

        // ----------------------------------------------------------------------------------

        void
        print()
        {
            this->print_cells();
            this->print_vertices();
            this->print_block_sets();
            this->print_groupings();
            this->print_vertex_ancestry();
        }

        // ----------------------------------------------------------------------------------

        void
        print_cells(
                bool        aOmitIndex = false,
                std::string aFile      = "" );

        // ----------------------------------------------------------------------------------

        void
        print_vertices(
                bool        aOmitIndex = false,
                std::string aFile      = "" );

        // ----------------------------------------------------------------------------------

        void
        print_block_sets()
        {
            const char separator = ' ';
            const int  nameWidth = 24;
            const int  numWidth  = 12;

            std::cout << "Num Block Sets: " << this->mBlockSetCellGroup.size() << std::endl;

            std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << "Name";
            std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Ordinal";
            std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Num Cells";
            std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << "Cells";
            std::cout << std::endl;

            for ( uint iBS = 0; iBS < mBlockSetCellGroup.size(); iBS++ )
            {
                std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetNames( iBS );
                std::cout << std::left << std::setw( numWidth ) << std::setfill( separator ) << mBlockSetLabelToOrd.find( mBlockSetNames( iBS ) )->second;
                std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetCellGroup( iBS )->mIgCellGroup.size();

                //
                for ( uint iIgCell = 0; iIgCell < mBlockSetCellGroup( iBS )->mIgCellGroup.size(); iIgCell++ )
                {
                    std::cout << std::left << std::setw( nameWidth ) << std::setfill( separator ) << mBlockSetCellGroup( iBS )->mIgCellGroup( iIgCell )->get_id();
                }

                std::cout << std::endl;
            }
        }

        // ----------------------------------------------------------------------------------

        void
        print_groupings( std::string aFile = "" );

        // ----------------------------------------------------------------------------------

        void
        print_vertex_ancestry() const
        {
            moris::print( mIgVertexParentEntityIndex, "Vertex Parent Entity Index" );
            moris::print( mIgVertexParentEntityRank, "Vertex Parent Entity Rank" );

            // moris::Cell< moris_index > mIgVertexParentEntityIndex;
            // moris::Cell< moris_index > mIgVertexParentEntityRank;
        }

        // ----------------------------------------------------------------------------------

        void
        trim_data();

        // ----------------------------------------------------------------------------------

        /**
         * @brief updates the memeber data mCommunicationTable, it is accessed by the basis processor object
         * 
         * @param aNewCommunicationTable 
         */
        void
        update_communication_table( moris::Cell< moris_id > const & aNewCommunicationTable );

      private:
        // ----------------------------------------------------------------------------------

        void
        create_base_cell_blocks();

        // ----------------------------------------------------------------------------------

        void
        setup_comm_map();

        // ----------------------------------------------------------------------------------

    };    // class Cut_Integration_Mesh

    // ----------------------------------------------------------------------------------

    class Child_Mesh_Experimental
    {
        friend class Cut_Integration_Mesh;
        friend class Integration_Mesh_Generator;

      public:
        moris::mtk::Cell*                  mParentCell;
        moris_index                        mChildMeshIndex;
        std::shared_ptr< IG_Cell_Group >   mIgCells;
        std::shared_ptr< IG_Vertex_Group > mIgVerts;

        // subphases
        moris::Cell< std::shared_ptr< IG_Cell_Group > > mSubphaseCellGroups;

      public:
        Child_Mesh_Experimental()
        {
        }

        moris_index
        get_parent_element_index()
        {
            return this->get_parent_cell()->get_index();
        }

        moris::mtk::Cell*
        get_parent_cell()
        {
            return mParentCell;
        }

        moris_index
        get_child_mesh_index()
        {
            return mChildMeshIndex;
        }

        void
        set_subphase_groups( moris::Cell< std::shared_ptr< IG_Cell_Group > >& aSubphasesGroups )
        {
            mSubphaseCellGroups = aSubphasesGroups;
        }

        uint
        get_num_subphase_cell_groups() const
        {
            return mSubphaseCellGroups.size();
        }

        std::shared_ptr< IG_Cell_Group >
        get_subphase_cell_group( moris_index aLocalSpIndex )
        {
            return mSubphaseCellGroups( aLocalSpIndex );
        }

    };    // class Child_Mesh_Experimental

    // ----------------------------------------------------------------------------------

}    // namespace xtk

#endif
