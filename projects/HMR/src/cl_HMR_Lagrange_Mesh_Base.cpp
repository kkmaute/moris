#include <cstdio>
#include <fstream>

#include "cl_Stopwatch.hpp" //CHR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src


namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------
//   public:
//------------------------------------------------------------------------------

        Lagrange_Mesh_Base::Lagrange_Mesh_Base (
                const Parameters     * aParameters,
                Background_Mesh_Base * aBackgroundMesh,
                BSpline_Mesh_Base    * aBSplineMesh,
                const uint           & aOrder ) :
                        Mesh_Base(
                                aParameters,
                                aBackgroundMesh,
                                aOrder ),
                         mBSplineMesh( aBSplineMesh )

        {
            // sanity check
            if ( aBSplineMesh != NULL )
            {
                MORIS_ERROR( aBSplineMesh->get_order() >= aOrder,
                        "Error while creating Lagrange mesh. Linked B-Spline mesh must have same or higher order.");
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::update_mesh()
        {
            // start timer
            tic tTimer;

            // activate pattern on background mesh
            this->select_activation_pattern();

            // tidy up memory
            this->delete_pointers();

            // create Lagrange Elements from Background Elements
            this->create_elements();

            // create nodes
            this->create_nodes();

            // update list of used nodes
            this->update_node_list();

            // link elements to B-Spline meshes
            if ( mBSplineMesh != NULL )
            {
                this->link_twins();
            }

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created Lagrange mesh.\n               Mesh has %lu active and refined elements and %lu nodes.\n               Creation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) this->get_number_of_elements(),
                        ( long unsigned int ) this->get_number_of_nodes_on_proc(),
                        ( double ) tElapsedTime / 1000 );
            }
        }

// ----------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::save_to_file( const std::string & aFilePath )
        {
            // get the file extension
            auto tFileExt = aFilePath.substr(aFilePath.find_last_of(".")+1,
                                             aFilePath.length() );

            // guess routine from extension
            if (tFileExt == "vtk")
            {
                this->save_to_vtk( aFilePath );
            }
            else if(tFileExt == "msh")
            {
                this->save_to_gmsh( aFilePath );
            }
            else
            {
                MORIS_ERROR( false,
                        "Wrong file type passed to Lagrange mesh.\n If you want to save an exodus file, create an MTK object first." );
            }
        }

//------------------------------------------------------------------------------
//   protected:
//------------------------------------------------------------------------------

        Mat< real > &
        Lagrange_Mesh_Base::create_field_data( const std::string & aLabel )
        {
            // first field is always element mesh
            mFieldLabels.push_back( aLabel );

            uint tIndex = mFieldData.size();

            // initialize empty matrix. It is populated later
            Mat< real > tEmpty;
            mFieldData.push_back( tEmpty );

            return mFieldData( tIndex );
        }

// -----------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::create_nodes_on_higher_levels()
        {
            // get max level of mest
            uint tMaxLevel = mBackgroundMesh->get_max_level();

            // loop over all levels
            for ( uint l=0; l<tMaxLevel; ++l )
            {
                // get all elements from this level
                Cell< Background_Element_Base* > tElements;

                mBackgroundMesh->collect_elements_on_level_including_aura(
                        l, tElements );

                // loop over all elements on this level
                for( auto tElement : tElements )
                {

                    // test if this element has children and is not padding
                    // and is refined
                    if ( tElement->has_children() && ! tElement->is_padding() &&
                            tElement->is_refined( mActivePattern ) )
                    {
                        // calculate nodes of children
                        mAllElementsOnProc( tElement->get_memory_index() )
                                ->create_basis_for_children(
                                        mAllElementsOnProc,
                                        mNumberOfAllBasis );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
//   private:
//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::create_nodes()
        {
            // nodes on first level are created separately
            this->create_nodes_on_level_zero();

            // it is easier to create all zero level nodes first
            // and remove obsolete ones. Don't worry, there are not that many.
            this->delete_nodes_on_padding_elements();

            // now we do the nodes on higher levels
            this->create_nodes_on_higher_levels();

            // fill cell with all remaining nodes
            this->collect_nodes();

            // calculate system wide IDs (helpful for debugging)
            this->calculate_node_ids();

            // find out node to element connectivity
            this->determine_elements_connected_to_basis();

            // and determine the node ownership to be the smallest proc ID
            // of all elements connected to this node
            this->guess_basis_ownership();

            // Make sure that node ownership is correct. Correct otherwise.
            this->confirm_basis_ownership();

            // calculate node coordinates with respect to user defined offset
            this->calculate_node_coordinates();

            // create node numbers
            this->calculate_node_indices();


        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::delete_nodes_on_padding_elements()
        {
            // get number of elements on coarsest level
            auto tNumberOfElements = mAllCoarsestElementsOnProc.size();

            // counter for nodes
            luint tCount = 0;

            // loop over all elements on coarsest level
            for( luint e=0; e<tNumberOfElements; ++e)
            {
                // get pointer to Lagrange element
                Element* tElement
                    = mAllCoarsestElementsOnProc( e );

                // test if element is not padding
                if ( ! tElement->is_padding() )
                {
                    // loop over all nodes of element
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to node
                        Basis* tNode = tElement->get_basis( k );

                        if ( ! tNode->is_flagged() )
                        {
                            // set basis as active
                            tNode->flag();

                            // increment counter
                            ++tCount;
                        }
                    }
                }
            }

            // ask background mesh for number of elements per direction
            Mat< luint > tNumberOfElementsPerDirection =
                mBackgroundMesh->get_number_of_elements_per_direction_on_proc();

            // assign Cell for nodes to be deleted
            Cell< Basis* > tNodes( mNumberOfAllBasis-tCount, nullptr );

            // reset counter
            tCount = 0;

            // loop over all elements on coarsest level
            for( luint e=0; e<tNumberOfElements; ++e)
            {
                // get pointer to Lagrange element
                Element* tElement
                    = mAllCoarsestElementsOnProc( e );

                // test if element is padding
                if ( tElement->is_padding() )
                {
                    // loop over all nodes of element
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to node
                        Basis* tNode = tElement->get_basis( k );

                        // test if node exists
                        if ( tNode != NULL )
                        {
                            // test if node is not flagged
                            if ( ! tNode->is_flagged() )
                            {
                                // flag node
                                tNode->flag();

                                // copy node to delete list
                                tNodes( tCount++ ) = tNode;
                            }
                        }
                    }
                }

            }

            // delete nodes
            for ( uint k=0; k<tCount; ++k )
            {
                delete tNodes( k );

                // decrement number of nodes
                --mNumberOfAllBasis;
            }

            // tidy up: unflag all remaining nodes
            // loop over all elements on coarsest level
            for( luint e=0; e<tNumberOfElements; ++e)
            {
                // get pointer to Lagrange element
                Element* tElement
                = mAllCoarsestElementsOnProc( e );

                // test if element is not padding
                if ( ! tElement->is_padding() )
                {
                    // loop over all nodes of element
                    for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to node
                        Basis* tNode = tElement->get_basis( k );

                        // unset flag
                        tNode->unflag();
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::collect_nodes()
        {
            // clear node list
            mAllBasisOnProc.clear();

            // reserve size
            mAllBasisOnProc.resize( mNumberOfAllBasis, nullptr );

            // initialize counter
            luint tCount = 0;

            // get number of active elements on proc
            luint tNumberOfElements = mBackgroundMesh
                ->get_number_of_active_elements_on_proc_including_aura();

            // get rank
            uint tMyRank = par_rank();

            // reset element counter
            mNumberOfElements = 0;

            // loop over all active elements on proc
            for ( luint e=0; e<tNumberOfElements; ++e )
            {
                // get pointer to background element
                Background_Element_Base* tBackElement = mBackgroundMesh
                        ->get_element_from_proc_domain_including_aura( e );

                // get pointer to Lagrange element
                Element* tElement
                    = mAllElementsOnProc( tBackElement->get_memory_index() );

                if ( ! tBackElement->is_deactive( mActivePattern )  )
                {
                    // flag nodes that are used by this proc
                    if ( tBackElement->get_owner() == tMyRank )

                    {
                        for ( uint k=0; k<mNumberOfBasisPerElement; ++k )
                        {
                            tElement->get_basis( k )->use();
                        }

                        // increment element counter
                        ++mNumberOfElements;
                    }

                    // loop over all nodes of this element
                    for ( uint k=0; k<mNumberOfBasisPerElement; ++k )
                    {
                        // get pointer to node
                        Basis* tNode = tElement->get_basis( k );

                        // test if node is flagged
                        if ( ! tNode->is_flagged() )
                        {
                            // set index in memory
                            tNode->set_memory_index( tCount );

                            // add node to list
                            mAllBasisOnProc( tCount ++ ) = tNode;

                            // flag node
                            tNode->flag();
                        }
                    }
                }
            }
            // make sure that number of nodes is correct
            MORIS_ERROR( tCount == mNumberOfAllBasis, "Number of Nodes does not match." );

           /* this is moved to calculate_node_indices
             // reset node counter
            mNumberOfUsedNodes = 0;

            // count number of nodes
            for( auto tNode : mAllBasisOnProc )
            {
                if ( tNode->is_used() )
                {
                    ++mNumberOfUsedNodes;
                }
            } */
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::calculate_node_ids()
        {

            switch ( mParameters->get_number_of_dimensions() )
            {
                case( 1 ):
                {
                    for( auto tNode  : mAllBasisOnProc )
                    {
                        // get ij position of node
                        const luint * tI = tNode->get_ijk();

                        // calculate ID and write to node
                        tNode->set_domain_id(
                                this->calculate_node_id(
                                        tNode->get_level(),
                                        tI[0] ) );
                    }

                    break;
                }
                case( 2 ):
                {
                    for( auto tNode  : mAllBasisOnProc )
                    {
                        // get ij position of node
                        const luint * tIJ = tNode->get_ijk();

                        // calculate ID and write to node
                        tNode->set_domain_id(
                                this->calculate_node_id(
                                        tNode->get_level(),
                                        tIJ[0],
                                        tIJ[1]) );
                    }

                    break;
                }
                case( 3 ):
                {
                    // 3D case
                    for( auto tNode  : mAllBasisOnProc )
                    {
                        // get ij position of node
                        const luint * tIJK = tNode->get_ijk();

                        // calculate ID and write to node
                        tNode->set_domain_id(
                                this->calculate_node_id(
                                        tNode->get_level(),
                                        tIJK[0],
                                        tIJK[1],
                                        tIJK[2]) );
                    }

                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Lagrange_Mesh: Invalid number of dimensions");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * calculates system wide unique node indices for MTK
         */
        void
        Lagrange_Mesh_Base::calculate_node_indices()
        {
            // reset node counters
            mNumberOfUsedAndOwnedNodes = 0;
            mNumberOfUsedNodes = 0;

            // get number of ranks
            uint tNumberOfProcs = par_size();

            // initialize local index of node
            // reset counter


            if( tNumberOfProcs == 1 ) // serial mode
            {
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if node is used by current setup
                    if ( tNode->is_used() )
                    {
                        // in serial local index and domain index are identical
                        tNode->set_domain_index( mNumberOfUsedAndOwnedNodes++ );
                        tNode->set_local_index( mNumberOfUsedNodes++ );
                    }
                }
            }
            else // parallel mode
            {

                // get my rank
                uint tMyRank = par_rank();

                for( auto tNode : mAllBasisOnProc )
                {
                    // test if node is used by current setup
                    if ( tNode->is_used() )
                    {
                        // test if node is owned
                        if ( tNode->get_owner() == tMyRank )
                        {
                            tNode->set_domain_index( mNumberOfUsedAndOwnedNodes++ );
                        }

                        // set local index of node
                        tNode->set_local_index( mNumberOfUsedNodes++ );
                    }

                    // make sure that this basis is not flagged
                    tNode->unflag();
                }

                // communicate number of owned nodes with other procs
                Mat< luint > tNodesOwnedPerProc
                    = comm_gather_and_broadcast( mNumberOfUsedAndOwnedNodes );

                // get proc neighbors from background mesh
                auto tProcNeighbors = mBackgroundMesh->get_proc_neigbors();

                // calculate node offset table
                Mat< luint > tNodeOffset( tNumberOfProcs, 1, 0 );
                for( uint p=1; p<tNumberOfProcs; ++p )
                {
                    tNodeOffset( p ) =   tNodeOffset( p-1 )
                                       + tNodesOwnedPerProc( p-1 );
                }

                // get my offset

                luint tMyOffset = tNodeOffset( tMyRank );

                // loop over all nodes on proc
                for( auto tNode : mAllBasisOnProc )
                {
                    // test if the is used and node belongs to me
                    if ( tNode->is_used() )
                    {
                        if ( tNode->get_owner() == tMyRank )
                        {
                            // set global node index
                            tNode->set_domain_index(
                                    tNode->get_domain_index()
                                    + tMyOffset );
                        }
                    }
                }

                // now the global node indices of used and owned nodes
                // must be communicated to the other procs

                // get number of proc neighbors
                uint tNumberOfProcNeighbors
                    = mBackgroundMesh->get_number_of_proc_neighbors();

                // create cell of matrices to send
                Mat< luint > tEmpty;
                Cell< Mat< luint > > tSendIndex( tNumberOfProcNeighbors, tEmpty );

                // loop over all proc neighbors
                for ( uint p = 0; p<tNumberOfProcNeighbors; ++p )
                {
                    if (    tProcNeighbors( p ) < tNumberOfProcs
                         && tProcNeighbors( p ) != tMyRank )
                    {
                        // cell containing node pointers
                        Cell< Basis* > tNodes;

                        // collect nodes within inverse aura
                        this->collect_basis_from_aura( p, 1, tNodes );

                        // initialize node counter
                        luint tCount = 0;

                        // loop over all nodes
                        for( auto tNode : tNodes )
                        {
                            // test if node belongs to me
                            if (  tNode->get_owner() == tMyRank )
                            {
                                // increment counter
                                ++tCount;
                            }
                        }

                        // assign memory for send matrix
                        tSendIndex( p ).set_size( tCount, 1 );

                        // reset counter
                        tCount = 0;

                        // loop over all nodes
                        for( auto tNode : tNodes )
                        {
                            // test if node belongs to me
                            if ( tNode->get_owner() == tMyRank )
                            {
                                // write index of node into array
                                tSendIndex( p )( tCount++ ) = tNode->get_domain_index();
                            }
                        }
                    } // end proc exists and is not me
                } // end loop over all procs

                // matrices to receive
                Cell< Mat< luint > > tReceiveIndex;

                // communicate ownership to neighbors
                communicate_mats(
                        tProcNeighbors,
                        tSendIndex,
                        tReceiveIndex );

                // loop over all proc neighbors
                for ( uint p = 0; p<tNumberOfProcNeighbors; ++p )
                {
                    // get rank of neighbor
                    auto tNeighborRank = tProcNeighbors( p );

                    if (    tNeighborRank < tNumberOfProcs
                            && tNeighborRank != tMyRank )
                    {
                        // cell containing node pointers
                        Cell< Basis* > tNodes;

                        // collect nodes within aura
                        this->collect_basis_from_aura( p, 0, tNodes );

                        // initialize node counter
                        luint tCount = 0;

                        // loop over all nodes
                        for( auto tNode : tNodes )
                        {
                            // test if this node belongs to neighbor
                            if ( tNode->get_owner() == tNeighborRank )
                            {
                                // assign index to node
                                tNode->set_domain_index(
                                        tReceiveIndex( p )( tCount++ ) );
                            }
                        }

                    }
                } // end loop over all procs
            } // end parallel
        }

//------------------------------------------------------------------------------

        Element *
        Lagrange_Mesh_Base::get_child( Element * aElement,
                const uint            & aChildIndex )
        {
            // get pointer to background element
            Background_Element_Base* tBackElement = aElement->get_background_element();

            if ( tBackElement->has_children() )
            {
                // get child of background element
                Background_Element_Base* tBackChild = tBackElement
                        ->get_child( aChildIndex );

                // grab child from element list
                return mAllElementsOnProc( tBackChild->get_memory_index() );
            }
            else
            {
                // return nothing
                return nullptr;
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::reset_fields()
        {
            mFieldLabels.clear();
            mFieldData.clear();

            // first field is always element mesh
            mFieldLabels.push_back("Element_Level");

            // initialize empty matrix. It is populated later
            Mat< real > tEmpty;
            mFieldData.push_back( tEmpty );

            // second field is always vertex IDs
            mFieldLabels.push_back("Vertex_IDs");
            mFieldData.push_back( tEmpty );

        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::add_field( const std::string & aLabel,
                                       const Mat< real > & aData )
        {
            mFieldLabels.push_back( aLabel );
            mFieldData.push_back( aData );
        }

//------------------------------------------------------------------------------

        MTK *
        Lagrange_Mesh_Base::create_mtk_object()
        {
            MORIS_ERROR( mOrder <= 2 , "Tried to create an MTK object for third or higher order. \n This is not supported by Exodus II.");



            // create new MTK object
            MTK* aMTK = new MTK( this );

            // create data
            aMTK->create_mesh_data();

            // return MTK object
            return aMTK;
        }

//------------------------------------------------------------------------------

        bool
        Lagrange_Mesh_Base::test_for_double_nodes()
        {
            // strategy: fill a matrix with node IDs. Make them unique.
            // each node must appear only once

            this->determine_elements_connected_to_basis();

            // get numnber of nodes
            luint tNumberOfNodes = mAllBasisOnProc.size();
            // matrix which will contain node IDs
            Mat< luint > tNodeIDs( tNumberOfNodes, 1 );

            // loop over all nodes
            for( luint k=0; k<tNumberOfNodes; ++k )
            {
                // get node
                Basis* tNode = mAllBasisOnProc( k );

                // get level of node
                luint tLevel = tNode->get_level();

                if ( tLevel == 0 )
                {
                    tNodeIDs( k ) = tNode->get_domain_id();
                }
                else
                {
                    // get ijk of node
                    const luint* tNodeIJK = tNode->get_ijk();

                    // copy array into writable array
                    luint tIJK[ 3 ];
                    for( uint i = 0; i<mNumberOfDimensions; ++i )
                    {
                        tIJK[ i ] = tNodeIJK[ i ];
                    }

                    bool tCheck = true;

                    // now see if there is any node above
                    while ( tLevel > 0 && tCheck )
                    {
                        for( uint i = 0; i<mNumberOfDimensions; ++i )
                        {
                            tCheck = tCheck && ( tIJK[ i ] % 2 == 0 );
                        }
                        if ( tCheck )
                        {
                            // go up
                            for( uint i = 0; i<mNumberOfDimensions; ++i )
                            {
                                tIJK[ i ] /= 2;
                            }

                            // decrement level
                            --tLevel;
                        }
                    }

                    // calculate new id of node
                    if ( mNumberOfDimensions == 1 )
                    {
                        tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ] );
                    }
                    else if ( mNumberOfDimensions == 2 )
                    {
                        tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ], tIJK[ 1 ] );
                    }
                    else if ( mNumberOfDimensions == 3 )
                    {
                        tNodeIDs( k ) =  this->calculate_node_id( tLevel, tIJK[ 0 ], tIJK[ 1 ], tIJK[ 2 ] );
                    }
                }
            }

            /*        // init a boost bitset
                    BoostBitset tBitset( tNodeIDs.max()+1 );

                    // loop over all nodes
                    for( luint k=0; k<tNumberOfNodes; ++k )
                    {
                        // write id of node into array

                        // test if ID was already flagged
                        if ( ! tBitset.test(  tNodeIDs( k ) ) )
                        {
                            tBitset.set( tNodeIDs( k ) );
                        }
                        else
                        {
                            // other coordinate
                            luint j=0;
                            // search for corresponding node
                            for( luint i=0; i<tNumberOfNodes; ++i )
                            {
                                if ( tNodeIDs( i ) == tNodeIDs( k ) )
                                {
                                    j = i;
                                    break;
                                }

                            }

                            std::fprintf( stdout,
                                    "Error: Node %lu has the ID %lu, which is already used bu node %lu.\n\n",
                                    ( long unsigned int ) k,
                                    ( long unsigned int ) tNodeIDs( k ),
                                    ( long unsigned int ) j
                            );

                            // print elements of node 1
                            std::fprintf( stdout, "Elements connected to Node %lu : \n", ( long unsigned int ) k );

                            Mat<luint> tElements = mElementsPerNode( mAllBasisOnProc( k )->get_domain_index() );
                             for( luint i=0; i<tElements.length(); ++i )
                            {
                                // get element
                                Element* tElement = mAllElementsOnProc ( tElements( i ) );

                                std::fprintf( stdout, "    Element %4lu    ID %4lu      Parent:  %4lu\n",
                                        ( long unsigned int ) tElement->get_background_element()->get_subdomain_index(),
                                        ( long unsigned int ) tElement->get_background_element()->get_domain_id(),
                                        ( long unsigned int ) tElement->get_background_element()->get_parent()->get_domain_id() );
                            }
                            std::fprintf( stdout, "\n" );

                            // print elements of node2
                            std::fprintf( stdout, "Elements connected to Node %lu : \n", ( long unsigned int ) j );

                            tElements = mElementsPerNode( mAllBasisOnProc( j )->get_domain_index() );

                            for( luint i=0; i<tElements.length(); ++i )
                            {
                                // get element
                                Element* tElement = mAllElementsOnProc ( tElements( i ) );

                                std::fprintf( stdout, "    Element %4lu    ID %4lu     Parent:  %4lu\n",
                                        ( long unsigned int ) tElement->get_background_element()->get_subdomain_index(),
                                        ( long unsigned int ) tElement->get_background_element()->get_domain_id(),
                                        ( long unsigned int ) tElement->get_background_element()->get_parent()->get_domain_id()  );

                            }
                            std::fprintf( stdout, "\n" );

                            exit(-1);
                        }

                    } */

            // make matrix unique
            tNodeIDs = unique( tNodeIDs );

            // make sure that number of nodes is the same
            return tNodeIDs.length() == tNumberOfNodes;
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::save_to_gmsh( const std::string & aFilePath )
        {
            // start timer
            tic tTimer;

            // get my rank
            uint tMyRank = par_rank();

            // modify filename
            std::string tFilePath;
            if ( moris::par_size() > 1 )
            {
                tFilePath = aFilePath + "." +  std::to_string( par_size() ) + "." +  std::to_string( tMyRank );
            }
            else
            {
                tFilePath = aFilePath;
            }

            // create output file
            std::FILE * tFile = std::fopen( tFilePath.c_str(), "w+");

            // write header
            std::fprintf( tFile, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

            // write node coordinates
            std::fprintf( tFile, "$Nodes\n%lu\n", ( long unsigned int ) mNumberOfUsedNodes );

            // get mesh scale factor
            real tScale = mParameters->get_gmsh_scale();

            switch ( mParameters->get_number_of_dimensions() )
            {
                case( 2 ) :
                {
                    // loop over all nodes
                    for( auto tNode : mAllBasisOnProc )
                    {
                        // test if this node is relevant
                        if ( tNode->is_used() )
                        {
                            // get coordinates of node
                            const real* tXY = tNode->get_xyz();

                            // write coordinates to ASCII file
                            std::fprintf( tFile,
                                    "%lu %.17f %.17f 0\n",
                                    ( long unsigned int ) tNode->get_domain_index()+1,
                                    ( double ) tXY[ 0 ]*tScale,
                                    ( double ) tXY[ 1 ]*tScale );
                        }
                    }
                    break;
                }
                case( 3 ) :
                {
                    // loop over all nodes
                    for( auto tNode : mAllBasisOnProc )
                    {
                        // test if this node is relevant
                        if ( tNode->is_used() )
                        {
                            // get coordinates of node
                            const real* tXYZ = tNode->get_xyz();

                            // write coordinates to ASCII file
                            std::fprintf( tFile,
                                    "%lu %.17f %.17f %.17f\n",
                                    ( long unsigned int ) tNode->get_domain_index()+1,
                                    ( double ) tXYZ[ 0 ]*tScale,
                                    ( double ) tXYZ[ 1 ]*tScale,
                                    ( double ) tXYZ[ 2 ]*tScale );
                        }
                    }
                    break;

               }
               default :
               {
                   MORIS_ERROR( false, "wrong number of dimensions\n");
                   break;
               }
            }

            // end node tag
            std::fprintf( tFile, "$EndNodes\n" );

            // write element tag
            std::fprintf( tFile, "$Elements\n%lu\n",
                    ( long unsigned int ) mNumberOfElements );

            // loop over all elements
            for( auto tElement: mAllElementsOnProc )
            {
                // test if this element is relevant
                if ( tElement->get_owner() == tMyRank && tElement->is_active() )
                {
                    // print element line
                    std::fprintf( tFile, "%lu %s\n",
                            ( long unsigned int ) tElement->get_domain_index()+1,
                            tElement->get_gmsh_string().c_str() );
                }
            }

            // finish element tag
            std::fprintf( tFile, "$EndElements\n" );

            // close file
            std::fclose( tFile );

            // print a debug statement if verbosity is set
            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created GMSH File: %s\n               Writing took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        tFilePath.c_str(),
                        ( double ) tElapsedTime / 1000 );
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::link_twins( )
        {
            // get number of elements of interest
            auto tNumberOfElements = this->get_number_of_elements();

            // loop over all elements of interest
            for( uint k=0; k<tNumberOfElements; ++k )
            {
                // get pointer to Lagrange element
                auto tLagrangeElement = this->get_element( k );

                // link elements
                tLagrangeElement->set_twin( mBSplineMesh->get_element( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::save_to_vtk( const std::string & aFilePath )
        {
            // start timer
            tic tTimer;

            // modify filename
            std::string tFilePath;
            if ( moris::par_size() > 1 )
            {
                auto tFileExt = aFilePath.substr(aFilePath.find_last_of("."),aFilePath.length());
                auto tBasePath = aFilePath.substr(0,aFilePath.find_last_of("."));
                tFilePath = tBasePath + "_" +  std::to_string(moris::par_rank()) + tFileExt;

            }
            else
            {
                tFilePath = aFilePath;
            }

            // open the file
            std::ofstream tFile(tFilePath, std::ios::binary);

            // containers
            //float tFValue = 0;
            //int   tIValue = 0;
            float tFChar = 0;
            int   tIChar = 0;

            tFile << "# vtk DataFile Version 3.0" << std::endl;
            tFile << "GO BUFFS!" << std::endl;
            tFile << "BINARY" << std::endl;
            //tFile << "ASCII" << std::endl;
            luint tNumberOfNodes = mAllBasisOnProc.size();

            // write node data
            tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;

            // ask settings for numner of dimensions
            auto tNumberOfDimensions = mParameters->get_number_of_dimensions();

            if ( tNumberOfDimensions == 2 )
            {
                // loop over all nodes
                for ( luint k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXY = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXY[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXY[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) 0 );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }
            else if ( tNumberOfDimensions == 3 )
            {
                // loop over all nodes
                for ( luint k = 0; k < tNumberOfNodes; ++k )
                {
                    // get coordinate from node
                    const real* tXYZ = mAllBasisOnProc( k )->get_xyz();

                    // write coordinates to mesh
                    tFChar = swap_byte_endian( (float) tXYZ[ 0 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 1 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                    tFChar = swap_byte_endian( (float) tXYZ[ 2 ] );
                    tFile.write( (char*) &tFChar, sizeof(float));
                }
            }

            tFile << std::endl;

            // write element topology
            int tCellType = mAllElementsOnProc( 0 )->get_vtk_type();

            // count number of non padding elements
            luint tNumberOfElements = 0;

            // can only write element data if vtk map exists
            if ( tCellType != 0 )
            {
                int tNumberOfNodesPerElement = swap_byte_endian( (int) mNumberOfBasisPerElement );


                luint tNumberOfAllElementsOnProc = mAllElementsOnProc.size();

                for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        // increment element counter
                        ++tNumberOfElements;
                    }
                }

                // write header for cells
                tFile << "CELLS " << tNumberOfElements << " "
                        << ( mNumberOfBasisPerElement + 1 )*tNumberOfElements  << std::endl;


                // matrix containing node indices
                Mat< luint > tNodes( mNumberOfBasisPerElement, 1 );

                // loop over all elements
                for( luint k=0; k<tNumberOfAllElementsOnProc; ++k )
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tFile.write( (char*) &tNumberOfNodesPerElement, sizeof(int)) ;

                        // ask element for nodes
                        mAllElementsOnProc( k )->get_basis_indices_for_vtk( tNodes );

                        for( uint i=0; i <mNumberOfBasisPerElement; ++i )
                        {
                            tIChar = swap_byte_endian( (int) tNodes( i ) );
                            tFile.write((char *) &tIChar, sizeof(int));
                            //tFile << " " << (int) mAllElementsOnProc( k )->get_basis( i )->get_memory_index();
                        }
                        //tFile << std::endl;
                    }
                }

                tFile << std::endl;

                // write cell types
                tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
                tIChar = swap_byte_endian( tCellType );
                for ( luint k = 0; k < tNumberOfElements; ++k)
                {
                    tFile.write( (char*) &tIChar, sizeof(int));
                    //tFile << tCellType << std::endl;
                }


                // write element data
                tFile << "CELL_DATA " << tNumberOfElements << std::endl;

                // write element ID
                tFile << "SCALARS ELEMENT_ID int" << std::endl;
                tFile << "LOOKUP_TABLE default" << std::endl;
                for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_domain_id() );
                        tFile.write( (char*) &tIChar, sizeof(int));
                    }
                }
                tFile << std::endl;

                // write proc owner
                tFile << "SCALARS ELEMENT_OWNER int" << std::endl;
                tFile << "LOOKUP_TABLE default" << std::endl;
                for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_owner() );
                        tFile.write( (char*) &tIChar, sizeof(float));
                    }
                }
                tFile << std::endl;

                // write level
                tFile << "SCALARS ELEMENT_LEVEL int" << std::endl;
                tFile << "LOOKUP_TABLE default" << std::endl;
                for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_level() );
                        tFile.write( (char*) &tIChar, sizeof(float));
                    }
                }
                tFile << std::endl;

                // write level
                tFile << "SCALARS ELEMENT_CHILD_INDEX int" << std::endl;
                tFile << "LOOKUP_TABLE default" << std::endl;
                for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_child_index() );
                        tFile.write( (char*) &tIChar, sizeof(float));
                    }
                }
                tFile << std::endl;

                // write memory index
                tFile << "SCALARS ELEMENT_MEMORY_INDEX int" << std::endl;
                tFile << "LOOKUP_TABLE default" << std::endl;
                for (moris::uint k = 0; k <  tNumberOfAllElementsOnProc; ++k)
                {
                    if ( mAllElementsOnProc( k )->is_active() )
                    {
                        tIChar = swap_byte_endian( (int) mAllElementsOnProc( k )->get_background_element()->get_memory_index() );
                        tFile.write( (char*) &tIChar, sizeof(float));
                    }
                }
                tFile << std::endl;
            }

            // write node data
            tFile << "POINT_DATA " << tNumberOfNodes << std::endl;

            tFile << "SCALARS NODE_ID int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_domain_id() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;

            tFile << "SCALARS NODE_OWNER int" << std::endl;
            tFile << "LOOKUP_TABLE default" << std::endl;
            for (moris::uint k = 0; k <  tNumberOfNodes; ++k)
            {

                tIChar = swap_byte_endian( (int) mAllBasisOnProc( k )->get_owner() );
                tFile.write( (char*) &tIChar, sizeof(float));
            }
            tFile << std::endl;


            // close the output file
            tFile.close();

            if ( mParameters->is_verbose() )
            {
                // stop timer
                real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"%s Created VTK debug file.\n               Mesh has %lu active and refined Elements and %lu Nodes.\n               Creation took %5.3f seconds.\n\n",
                        proc_string().c_str(),
                        ( long unsigned int ) tNumberOfElements,
                        ( long unsigned int ) tNumberOfNodes,
                        ( double ) tElapsedTime / 1000 );
            }

        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::create_nodes_on_level_zero()
        {
            // ask background mesh for number of dimensions
            Mat< luint > tNumberOfElements
            = mBackgroundMesh->get_number_of_elements_per_direction_on_proc();

            if( mNumberOfDimensions == 2 )
            {
                // get maximum numbers in i direction
                luint tImax = tNumberOfElements( 0, 0 );

                // get maximum numbers in j direction
                luint tJmax = tNumberOfElements( 1, 0 );

                // initialize element counter
                luint tCount = 0;

                // loop over all elements
                for( luint j=0; j<tJmax; ++j )
                {
                    for ( luint i=0; i<tImax; ++i )
                    {
                        mAllCoarsestElementsOnProc( tCount++ )
                            ->create_basis_on_level_zero(
                                mAllElementsOnProc,
                                mNumberOfAllBasis );
                    }
                }
            }
            else if( mNumberOfDimensions == 3 )
            {
                // get maximum numbers in i direction
                luint tImax = tNumberOfElements( 0, 0 );

                // get maximum numbers in j direction
                luint tJmax = tNumberOfElements( 1, 0 );

                // get maximum numbers in k direction
                luint tKmax = tNumberOfElements( 2, 0 );

                // initialize element counter
                luint tCount = 0;

                // loop over all elements
                for( luint k=0; k<tKmax; ++k )
                {
                    for( luint j=0; j<tJmax; ++j )
                    {
                        for ( luint i=0; i<tImax; ++i )
                        {
                            mAllCoarsestElementsOnProc( tCount++ )
                                ->create_basis_on_level_zero(
                                    mAllElementsOnProc,
                                    mNumberOfAllBasis );
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Lagrange_Mesh_Base::update_node_list()
        {
            // tidy up memory
            mNodes.clear();

            // assign memory
            mNodes.resize( mNumberOfUsedNodes, nullptr );

            // initialize counter
            luint tCount = 0;

            for( auto tNode : mAllBasisOnProc )
            {
                if ( tNode->is_used() )
                {
                    mNodes( tCount++ ) = tNode;
                }
            }

            MORIS_ERROR( tCount == mNumberOfUsedNodes, "Number Of used Nodes does not match" );
        }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
