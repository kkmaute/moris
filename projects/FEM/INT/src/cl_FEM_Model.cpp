
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"              //FEM/INT/src

#include "cl_FEM_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Set_User_Info.hpp"

#include "cl_MSI_Equation_Object.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    FEM_Model::FEM_Model(       mtk::Mesh_Manager *                 aMeshManager,
                          const moris_index                       & aMeshPairIndex,
                                moris::Cell< fem::Set_User_Info > & aSetInfo ) : mMeshManager( aMeshManager ),
                                                                                 mMeshPairIndex( aMeshPairIndex )
    {
        // get the number of sets
        uint tNumFemSets = aSetInfo.size();

        // start timer
        tic tTimer1;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 0: initialize
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Get pointers to interpolation and integration mesh
        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMeshManager->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create nodes
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // ask mesh about number of IP nodes on proc
        luint tNumOfIPNodes = tInterpolationMesh->get_num_nodes();

        // create IP node objects
        mIPNodes.resize( tNumOfIPNodes, nullptr );

        for( uint iNode = 0; iNode < tNumOfIPNodes; iNode++ )
        {
            mIPNodes( iNode ) = new fem::Node( &tInterpolationMesh->get_mtk_vertex( iNode ) );
        }

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( "Model: created %u FEM IP nodes in %5.3f seconds.\n\n",
                            ( unsigned int ) tNumOfIPNodes,
                            ( double ) tElapsedTime / 1000 );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create elements
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // start timer
        tic tTimer2;

        // create equation objects
        mFemSets.resize( tNumFemSets, nullptr );     // FIXME try to create them as equation sets

        // get the number of element to create
        mFemClusters.reserve( 100000 ); //FIXME

        //------------------------------------------------------------------------------
        // init the fem set counter
        moris::uint tFemSetCounter = 0;

        // loop over the used fem set
        for( luint iSet = 0; iSet < tNumFemSets; iSet++ )
        {
            mMeshSetToFemSetMap[ aSetInfo( iSet ).get_mesh_index() ] = tFemSetCounter;

            // create a list of clusters
            moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( aSetInfo( iSet ).get_mesh_index() );

            // get global number of clusters
            moris::uint tNumClusters = 0;
            moris::sum_all(tMeshSet->get_num_clusters_on_set(),tNumClusters);

            if (  tNumClusters !=0 )
            {
                // create new fem set
                mFemSets( tFemSetCounter ) = new fem::Set( this,
                                                           tMeshSet,
                                                           aSetInfo( iSet ),
                                                           mIPNodes );
            }
            else
            {
                // FIXME why do we build empty set?
                mFemSets( tFemSetCounter ) = new fem::Set();
            }

            // collect equation objects associated with the block-set
            mFemClusters.append( mFemSets( tFemSetCounter )->get_equation_object_list() );

            // update fem set counter
            tFemSetCounter++;
        }
        mFemClusters.shrink_to_fit();

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

            // print output

            MORIS_LOG_INFO( "Model: created %u FEM elements in %5.3f seconds.\n\n",
                    ( unsigned int ) mFemClusters.size(),
                    ( double ) tElapsedTime / 1000 );
        }
    }

//------------------------------------------------------------------------------

        FEM_Model::~FEM_Model()
        {
            // delete fem nodes
            for( auto tIPNodes : mIPNodes )
            {
                delete tIPNodes;
            }
            mIPNodes.clear();

            // delete the fem sets
            for( auto tFemSet : mFemSets )
            {
                delete tFemSet;
            }
            mFemSets.clear();

            // delete the fem cluster
            mFemClusters.clear();
        }

//------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
