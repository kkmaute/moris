/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Interpolation_Mesh_Editor.cpp
 *
 */

#include "cl_MTK_Interpolation_Mesh_Editor.hpp"
#include "cl_MTK_Mesh_DataBase_IP.hpp"
#include "cl_Tracer.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Vertex_Interpolation_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_Tracer.hpp"

namespace moris::mtk
{
    // ----------------------------------------------------------------------------
    Interpolation_Mesh_Editor::Interpolation_Mesh_Editor( moris::mtk::Interpolation_Mesh& aMTKMesh )
            : mInputMesh( aMTKMesh )
    {
        mIPMeshInfo = new Interpolation_Mesh_Info;

        // copy the necessary data for the old mesh
        mIPMeshInfo->mVertices               = mInputMesh.get_all_vertices();
        mIPMeshInfo->mNumInterpolations      = mInputMesh.get_num_interpolations();
        mIPMeshInfo->mSpatialDim             = mInputMesh.get_spatial_dim();
        mIPMeshInfo->mNumCells               = mInputMesh.get_num_elems();
        mIPMeshInfo->mNumLocalInterpolations = mInputMesh.get_enriched_mesh_indices().numel();
    }

    // ----------------------------------------------------------------------------

    Interpolation_Mesh_Editor::~Interpolation_Mesh_Editor()
    {
        delete mIPMeshInfo;
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::generate_vertex_data()
    {
        // initialize the sizes
        this->initialize_vertex_data();

        // loop over the lagrange vertices
        for ( uint iVertex = 0; iVertex < mIPMeshInfo->mVertices.size(); iVertex++ )
        {
            // get the matrix coordinates
            Matrix< DDRMat > tVertexCoords = mIPMeshInfo->mVertices( iVertex )->get_coords();

            // copy the matrix data
            std::copy( tVertexCoords.data(), tVertexCoords.data() + mIPMeshInfo->mSpatialDim, mOutputMesh->mVertexCoordinates.colptr( iVertex ) );

            // loop over the ip vertices corresponding a vertex
            for ( uint iOrder = 0; iOrder < mIPMeshInfo->mNumLocalInterpolations; iOrder++ )
            {
                // this converts the local consecutive mesh order to the bspline mesh index observed in global scheme
                uint tGlobalOrder = mOutputMesh->mMeshIndices( iOrder );

                // put the interpolation vertex in the list
                mIPMeshInfo->mVertexInterpolations.push_back( mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder ) );

                // obtain t-matrix info
                Matrix< DDRMat > const * tWeights = mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder )->get_weights();
                Matrix< IdMat >          tIds     = mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder )->get_ids();
                Matrix< IndexMat >       tIndices = mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder )->get_indices();
                Matrix< IdMat >          tOwners  = mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder )->get_owners();

                // insert t-matrix info into a big cell
                mOutputMesh->mWeights( iOrder ).insert( mOutputMesh->mOffSetTMatrix( iOrder )( iVertex ), tWeights->cbegin(), tWeights->cend() );

                mOutputMesh->mBasisIds( iOrder ).insert( mOutputMesh->mOffSetTMatrix( iOrder )( iVertex ), tIds.begin(), tIds.end() );

                mOutputMesh->mBasisIndices( iOrder ).insert( mOutputMesh->mOffSetTMatrix( iOrder )( iVertex ), tIndices.begin(), tIndices.end() );

                mOutputMesh->mBasisOwners( iOrder ).insert( mOutputMesh->mOffSetTMatrix( iOrder )( iVertex ), tOwners.begin(), tOwners.end() );
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::initialize_vertex_data()
    {
        uint tNumInterpolations = mIPMeshInfo->mNumLocalInterpolations;

        // set the size of mOffSetTMatrix
        mOutputMesh->mOffSetTMatrix.resize( tNumInterpolations );
        mOutputMesh->mWeights.resize( tNumInterpolations );
        mOutputMesh->mBasisIds.resize( tNumInterpolations );
        mOutputMesh->mBasisOwners.resize( tNumInterpolations );
        mOutputMesh->mBasisIndices.resize( tNumInterpolations );

        // This loop is strictly to determine the size
        //  loop over the interpolation orders ( bspline meshes )
        for ( uint iLocalOrder = 0; iLocalOrder < tNumInterpolations; iLocalOrder++ )
        {
            // set the size of the offset cell for each bspline mesh
            mOutputMesh->mOffSetTMatrix( iLocalOrder ).resize( mIPMeshInfo->mVertices.size() + 1 );

            // initialize the offset value
            size_t tOffSet = 0;

            // first entry in the offset vector is always zero
            mOutputMesh->mOffSetTMatrix( iLocalOrder )( 0 ) = tOffSet;

            // loop ove the Lagrange nodes
            for ( uint iVertex = 0; iVertex < mIPMeshInfo->mVertices.size(); iVertex++ )
            {
                // this converts the local consecutive mesh order to the bspline mesh index observed in global scheme
                uint tGlobalOrder = mOutputMesh->mMeshIndices( iLocalOrder );

                tOffSet += mIPMeshInfo->mVertices( iVertex )->get_interpolation( tGlobalOrder )->get_weights()->n_rows();
                mOutputMesh->mOffSetTMatrix( iLocalOrder )( iVertex + 1 ) = tOffSet;
            }

            // reserve the vector of the weight, ids, indices, owners
            mOutputMesh->mWeights( iLocalOrder ).reserve( tOffSet );
            mOutputMesh->mBasisIds( iLocalOrder ).reserve( tOffSet );
            mOutputMesh->mBasisOwners( iLocalOrder ).reserve( tOffSet );
            mOutputMesh->mBasisIndices( iLocalOrder ).reserve( tOffSet );
        }

        // set the coordinate matrix size
        mOutputMesh->mVertexCoordinates.set_size( mIPMeshInfo->mSpatialDim, mIPMeshInfo->mVertices.size() );

        // reserve space for ip vertices
        mIPMeshInfo->mVertexInterpolations.reserve( mIPMeshInfo->mNumLocalInterpolations * mIPMeshInfo->mVertices.size() );
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::generate_cell_data()
    {
        this->initialize_cell_data();

        // loop over the cells
        for ( size_t iCell = 0; iCell < mIPMeshInfo->mNumCells; iCell++ )
        {
            // get the mtk cell
            const mtk::Cell& tCell = mInputMesh.get_mtk_cell( iCell );

            // get vertex pointers
            Vector< Vertex* > tVertexPointers = tCell.get_vertex_pointers();

            // transfer pointers to the indices
            std::transform( tVertexPointers.begin(),
                    tVertexPointers.end(),
                    mIPMeshInfo->mCellToVertexIndicies.begin() + mOutputMesh->mCellToVertexOffSet( iCell ),
                    []( Vertex*& aVertex ) -> moris_index { return aVertex->get_index(); } );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::initialize_cell_data()
    {
        // get number of vertices and number of vertices for each cell
        uint tNumCellVertex = mInputMesh.get_mtk_cell( 0 ).get_number_of_vertices();

        // right now we assume that all the IP cells have the same toplogy and order
        mIPMeshInfo->mCellToVertexIndicies.resize( mIPMeshInfo->mNumCells * tNumCellVertex );

        // initialize the offset
        mOutputMesh->mCellToVertexOffSet.resize( mIPMeshInfo->mNumCells + 1 );

        size_t tOffSet                        = 0;
        mOutputMesh->mCellToVertexOffSet( 0 ) = tOffSet;

        // loop over the offsets to populate the data
        for ( uint iOffSet = 0; iOffSet < mIPMeshInfo->mNumCells; iOffSet++ )
        {
            tOffSet += tNumCellVertex;
            mOutputMesh->mCellToVertexOffSet( iOffSet + 1 ) = tOffSet;
        }
    }

    // ----------------------------------------------------------------------------

    Interpolation_Mesh_DataBase_IP*
    Interpolation_Mesh_Editor::perform()
    {
        Tracer tTracer( "MTK", "IP mesh building", "Build" );

        mOutputMesh = new Interpolation_Mesh_DataBase_IP( mInputMesh.get_background_mesh() );

        this->create_enriched_mesh_indices();

        this->generate_vertex_data();

        this->generate_cell_data();

        // create the vertex interpolaions
        this->create_ip_vertices();

        // create the vertecies
        this->create_vertices();

        // create the vertecies
        this->create_cells();

        // create the communication table
        this->create_communication_table();

        // create the vertex global id to index map for GEN
        this->create_vertex_glb_id_to_loc_vertex_ind_map();

        // create the adof map, needed for FEM analysis
        this->create_adof_map();

        // check that input and output mesh are compatible
        this->check_input_output_mesh();

        return mOutputMesh;
    }

    // ----------------------------------------------------------------------------

    moris::Memory_Map
    Interpolation_Mesh_Editor::get_memory_usage()
    {
        moris::Memory_Map tMemoryMap;

        return tMemoryMap;
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::free_memory()
    {
        delete mIPMeshInfo;
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_enriched_mesh_indices()
    {
        mOutputMesh->mMeshIndices = mInputMesh.get_enriched_mesh_indices();

        // generate the inverse map that corresponds to global mesh index to local mesh index in the database
        for ( uint iLocalOrder = 0; iLocalOrder < mOutputMesh->mMeshIndices.numel(); iLocalOrder++ )
        {
            mOutputMesh->mGlobalMeshIndexToLocalMeshIndex[ mOutputMesh->mMeshIndices( iLocalOrder ) ] = iLocalOrder;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_ip_vertices()
    {
        // create vertex interpolations
        mOutputMesh->mVertexInterpoltions.resize( mIPMeshInfo->mVertexInterpolations.size() );

        // loop over the order of interpolations
        for ( uint iOrder = 0; iOrder < mIPMeshInfo->mNumLocalInterpolations; iOrder++ )
        {
            // loop over the lagrange vertices
            for ( size_t iVertex = 0; iVertex < mIPMeshInfo->mVertices.size(); iVertex++ )
            {
                // construct the vertices
                mOutputMesh->mVertexInterpoltions( mIPMeshInfo->mNumLocalInterpolations * iVertex + iOrder ) =
                        Vertex_Interpolation_DataBase( iVertex,
                                iOrder,
                                mOutputMesh );
            }
        }

        // convert vertex interpolations to pointers
        mOutputMesh->mVertexInterpoltionsPtrs.reserve( mIPMeshInfo->mVertexInterpolations.size() );

        // fill in the vertex interpolation values
        for ( auto& iVertex : mOutputMesh->mVertexInterpoltions )
        {
            mOutputMesh->mVertexInterpoltionsPtrs.push_back( &iVertex );
        }
    }

    //-----------------------------------------------------------------------------
    void
    Interpolation_Mesh_Editor::create_vertices()
    {
        // allocate size for the vertices
        mOutputMesh->mVertices.reserve( mIPMeshInfo->mVertices.size() );
        mOutputMesh->mVertexIdList.reserve( mIPMeshInfo->mVertices.size() );
        mOutputMesh->mVertexOwnerList.reserve( mIPMeshInfo->mVertices.size() );

        // counter for the vertices
        uint iCounter = 0;

        // loop over old vertices to transfer data
        for ( const auto& iVertex : mIPMeshInfo->mVertices )
        {
            // Assert to ensure consecutive vertices
            MORIS_ASSERT( iVertex->get_index() == (moris_index)iCounter, "Index alignment issue in vertices" );

            // constrcut the vertex and put it in the list
            mOutputMesh->mVertices.emplace_back( Vertex_DataBase( iCounter,
                    mOutputMesh ) );

            // increase the count
            iCounter++;

            mOutputMesh->mVertexIdList.push_back( iVertex->get_id() );
            mOutputMesh->mVertexOwnerList.push_back( iVertex->get_owner() );
        }
    }

    //-----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_cells()
    {
        // create pointers with the parent class for the connectivity
        mOutputMesh->mCellToVertices.reserve( mIPMeshInfo->mCellToVertexIndicies.size() );

        // populate the cell to vertex connectivity list
        for ( const int& iVertex : mIPMeshInfo->mCellToVertexIndicies )
        {
            mOutputMesh->mCellToVertices.push_back( &mOutputMesh->mVertices( iVertex ) );
        }

        // FIXME : we assume all IP cells have the same topology
        // get a refernce cell to the get the IP cell geometry
        mtk::Cell const & tCell = mInputMesh.get_mtk_cell( 0 );

        // FIXME: will be changed in the xtk_refactor branch
        //  enum CellTopology tCellTopology = tCell.get_cell_info()->get_cell_topology();
        //  moris::mtk::Cell_Info_Factory tFactory;
        //  mCellInfo = tFactory.create_cell_info_sp( tCellTopology );

        // set the interpolation mesh cell info and create a cell info
        mtk::Cell_Info_Factory tFactory;
        mOutputMesh->mCellInfo = tFactory.create_cell_info_sp( tCell.get_geometry_type(), tCell.get_interpolation_order() );

        // reserve enough space for cells
        mOutputMesh->mCells.reserve( mIPMeshInfo->mNumCells );
        mOutputMesh->mCellIdList.reserve( mIPMeshInfo->mNumCells );
        mOutputMesh->mCellOwnerList.reserve( mIPMeshInfo->mNumCells );

        // loop over the old cells to create new cells
        for ( size_t iCell = 0; iCell < mIPMeshInfo->mNumCells; iCell++ )
        {
            mOutputMesh->mCells.push_back( Cell_DataBase( mInputMesh.get_mtk_cell( iCell ),
                    mOutputMesh->mCellInfo,
                    iCell,
                    mOutputMesh ) );

            // get cell id and index
            mOutputMesh->mCellIdList.push_back( mInputMesh.get_mtk_cell( iCell ).get_id() );
            mOutputMesh->mCellOwnerList.push_back( mInputMesh.get_mtk_cell( iCell ).get_owner() );
        }
    }

    //-----------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_communication_table()
    {
        // copy the communication table
        mOutputMesh->mCommunicationTable = mInputMesh.get_communication_table();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_vertex_glb_id_to_loc_vertex_ind_map()
    {
        // reserve the space for unordered map
        mOutputMesh->mVertexGlobalIdToLocalIndex.reserve( mIPMeshInfo->mVertices.size() );

        // create the vertex map used in gen based on the new vertex
        for ( uint iCounter = 0; iCounter < mOutputMesh->mVertices.size(); ++iCounter )
        {
            MORIS_ASSERT( mOutputMesh->mVertices( iCounter ).get_index() == iCounter, "Index alignment issue in vertices" );

            mOutputMesh->mVertexGlobalIdToLocalIndex[ mOutputMesh->mVertices( iCounter ).get_id() ] = iCounter;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::create_adof_map()
    {
        // resize enough space fo the maps on different mesh
        mOutputMesh->mAdofMap.resize( mIPMeshInfo->mNumInterpolations );

        // loop over the bspline meshes to get the adof map
        // this is copied from the old mesh directly because we don't assign id and index
        // for the vertex interpolations
        for ( uint iBSpline = 0; iBSpline < mIPMeshInfo->mNumInterpolations; iBSpline++ )
        {
            mInputMesh.get_adof_map( iBSpline, mOutputMesh->mAdofMap( iBSpline ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_Mesh_Editor::check_maps()
    {
        // obtain vertex map from the mesh
        std::unordered_map< moris_id, moris_index > tVertexGlobalIdToLocalIndex = mInputMesh.get_vertex_glb_id_to_loc_vertex_ind_map();

        // compare the vertex map
        bool tVertexMapEqual = mOutputMesh->mVertexGlobalIdToLocalIndex == tVertexGlobalIdToLocalIndex;

        if ( !tVertexMapEqual )
        {
            MORIS_LOG_INFO( "Vertex Global To Local Map does not match!" );
            return tVertexMapEqual;
        }

        // initialize the adof map equality to be true
        bool tAdofMapEqual = true;

        // initialize the adof map for the old mesh
        map< moris_id, moris_index > tAdofMapOldMesh;

        // loop over bspline meshes to
        for ( uint iOrder = 0; iOrder < mOutputMesh->mAdofMap.size(); iOrder++ )
        {
            // populate the adof map for the old mesh
            mInputMesh.get_adof_map( iOrder, tAdofMapOldMesh );

            // check to see if they are equal
            tAdofMapEqual = tAdofMapOldMesh.data() == mOutputMesh->mAdofMap( iOrder ).data();

            if ( !tAdofMapEqual )
            {
                MORIS_LOG_INFO( "Adof map is mismacthing with the old mesh for bspine mesh number %u!", iOrder );
                return tAdofMapEqual;
            }
        }

        return tVertexMapEqual && tAdofMapEqual;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_Mesh_Editor::check_vertices()
    {

        // check to see if the vertices are stored in consecutive manner and have the same id and indices
        bool tVertexIdAndIndexEqual = std::equal( mOutputMesh->mVertices.begin(),
                mOutputMesh->mVertices.end(),
                mIPMeshInfo->mVertices.begin(),
                []( Vertex_DataBase a, mtk::Vertex const * b ) -> bool { return a.get_id() == b->get_id() and a.get_index() == b->get_index(); } );

        // check if old vertices and new vertices have the same coords
        bool tEqualCoords = std::equal( mOutputMesh->mVertices.begin(),
                mOutputMesh->mVertices.end(),
                mIPMeshInfo->mVertices.begin(),
                []( Vertex_DataBase a, mtk::Vertex const * b ) -> bool {
                    return std::equal( a.get_coords().begin(), a.get_coords().end(), b->get_coords().begin() );
                } );

        return ( tEqualCoords and tVertexIdAndIndexEqual );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_Mesh_Editor::check_t_matrices()
    {
        bool tOutput = true;
        for ( const auto& iGloablOrder : mOutputMesh->mMeshIndices )
        {
            for ( uint i = 0; i < mOutputMesh->mVertices.size(); i++ )
            {
                // get member data for the old and new mesh , t-matrices
                Matrix< IndexMat > tBasisIndicesOldMesh = mInputMesh.get_mtk_vertex( i ).get_interpolation( iGloablOrder )->get_indices();
                Matrix< IndexMat > tBaisIndicesNewMesh  = mOutputMesh->mVertices( i ).get_interpolation( iGloablOrder )->get_indices();

                Matrix< IdMat > tBasisIdsOldMesh = mInputMesh.get_mtk_vertex( i ).get_interpolation( iGloablOrder )->get_ids();
                Matrix< IdMat > tBasisIdsNewMesh = mOutputMesh->mVertices( i ).get_interpolation( iGloablOrder )->get_ids();

                Matrix< IdMat > tBasisOwnerOldMesh = mInputMesh.get_mtk_vertex( i ).get_interpolation( iGloablOrder )->get_owners();
                Matrix< IdMat > tBasisOwnerNewMesh = mOutputMesh->mVertices( i ).get_interpolation( iGloablOrder )->get_owners();

                const Matrix< DDRMat >* tBasisWeightsOldMesh = mInputMesh.get_mtk_vertex( i ).get_interpolation( iGloablOrder )->get_weights();
                const Matrix< DDRMat >* tBasisWightsNewMesh  = mOutputMesh->mVertices( i ).get_interpolation( iGloablOrder )->get_weights();

                // check if t-matrices are equal
                bool tEqualIndex  = std::equal( tBasisIndicesOldMesh.begin(), tBasisIndicesOldMesh.end(), tBaisIndicesNewMesh.begin() );
                bool tEqualId     = std::equal( tBasisIdsOldMesh.begin(), tBasisIdsOldMesh.end(), tBasisIdsNewMesh.begin() );
                bool tEqualOwner  = std::equal( tBasisOwnerOldMesh.begin(), tBasisOwnerOldMesh.end(), tBasisOwnerNewMesh.begin() );
                bool tEqualWeight = std::equal( tBasisWeightsOldMesh->cbegin(), tBasisWeightsOldMesh->cend(), tBasisWightsNewMesh->cbegin() );

                tOutput = tEqualIndex && tEqualId && tEqualOwner && tEqualWeight;

                if ( !tOutput )
                {
                    MORIS_LOG_ERROR( "T-matrices are not matching for vertex %u of the mesh, the bspline local mesh index is %u", i, (uint)0 );
                    return tOutput;
                }
            }
        }
        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Interpolation_Mesh_Editor::check_cells()
    {
        bool tOutput = true;

        uint tNumCells = mInputMesh.get_num_elems();

        // loop over the old cells to create new cells
        for ( size_t iCell = 0; iCell < tNumCells; iCell++ )
        {
            // get the cells of the old and new mesh
            mtk::Cell const & tCellOldMesh = mInputMesh.get_mtk_cell( iCell );
            mtk::Cell const & tCellNewMesh = mOutputMesh->get_mtk_cell( iCell );

            // check to cells have the same id
            bool tSameCellId = tCellOldMesh.get_id() == tCellNewMesh.get_id();

            // get the vertices of the old and new mesh
            Vector< Vertex* > tOldVertices = tCellOldMesh.get_vertex_pointers();
            Vector< Vertex* > tNewVertices = tCellNewMesh.get_vertex_pointers();

            // coordinate of the cells for the old and new mesh
            Matrix< DDRMat > tVertexCoordsOld = tCellOldMesh.get_vertex_coords();
            Matrix< DDRMat > tVertexCoordsNew = tCellNewMesh.get_vertex_coords();

            // vertices have the same id and index
            bool tVertexIdAndIndexEqual = std::equal( tOldVertices.begin(),
                    tOldVertices.end(),
                    tNewVertices.begin(),
                    []( Vertex* a, Vertex* b ) -> bool {
                        return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                    } );

            bool tCoordsEqual = std::equal( tVertexCoordsOld.begin(), tVertexCoordsOld.end(), tVertexCoordsNew.begin() );

            // combine the outputs
            tOutput = tCoordsEqual && tSameCellId && tVertexIdAndIndexEqual && tOldVertices.size() == tNewVertices.size();

            if ( !tOutput )
            {
                MORIS_LOG_ERROR( "Cell number %zu of the mesh, is inconsistent, tCoordsEqual: %d ,tSameCellId: %d,tVertexIdAndIndexEqual: %d  ",
                        iCell,
                        tCoordsEqual,
                        tSameCellId,
                        tVertexIdAndIndexEqual );
                return tOutput;
            }
        }
        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_Mesh_Editor::check_input_output_mesh()
    {
        MORIS_ASSERT( this->check_vertices(), "Replaced vertices have different id and indices" );

        MORIS_ASSERT( this->check_maps(), "Maps of adof and vertex global to local are the same" );

        MORIS_ASSERT( this->check_t_matrices(), " T matrices are not the same for the newly created ip vertices" );

        MORIS_ASSERT( this->check_cells(), " Cells are not matching" );
    }
}    // namespace moris::mtk
