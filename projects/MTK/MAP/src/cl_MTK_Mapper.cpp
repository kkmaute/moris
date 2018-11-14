

#include "assert.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "cl_MSI_Solver_Interface.hpp"

#include "cl_FEM_Element.hpp"
#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_Node.hpp"

#include "op_elemwise_mult.hpp"
#include "op_div.hpp"
#include "fn_dot.hpp"

#include "fn_sum.hpp"
#include "fn_print.hpp"

#include "cl_MDL_Model.hpp"

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

namespace moris
{
    namespace mapper
    {

//------------------------------------------------------------------------------

        Mapper::Mapper( std::shared_ptr< mtk::Mesh > aMesh ) :
                mSourceMesh( aMesh ),
                mTargetMesh( aMesh )
        {

        }

//------------------------------------------------------------------------------

        Mapper::~Mapper()
        {
            // test if model and IWG have been created
            if( mHaveIwgAndModel )
            {
                // delete the fem model
                delete mModel;

                // delete IWG object
                delete mIWG;
            }

            // delete nodes for the filter
            if( mHaveNodes )
            {
                for( Node * tNode : mNodes )
                {
                    delete tNode;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_iwg_and_model()
        {
            if( ! mHaveIwgAndModel )
            {
                // create IWG object
                mIWG = new moris::fem::IWG_L2( );

                // create model
                mModel = new mdl::Model( mTargetMesh.get(), mIWG );

                mHaveIwgAndModel = true;
            }
        }

//-----------------------------------------------------------------------------

        void
        Mapper::perform_mapping(
                const std::string & aSourceLabel,
                const EntityRank    aSourceEntityRank,
                const std::string & aTargetLabel,
                const EntityRank    aTargetEntityRank )
        {
            // get index of source
            moris_index tSourceIndex = mSourceMesh->get_field_ind(
                    aSourceLabel,
                    aSourceEntityRank );

            MORIS_ERROR( tSourceIndex != gNoIndex, "perform_mapping() Source Field not found");

            // get target index
            moris_index tTargetIndex = mTargetMesh->get_field_ind(
                    aTargetLabel,
                    aTargetEntityRank );

            // test if output field has to be initialized
            if( tTargetIndex == gNoIndex )
            {
                // create target field
                tTargetIndex =
                        mTargetMesh->create_scalar_field(
                                aTargetLabel,
                                aTargetEntityRank );

            }

            switch( aSourceEntityRank )
            {
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                case( EntityRank::NODE ) :
                {
                   switch( aTargetEntityRank )
                   {
                       case( EntityRank::BSPLINE_1 ) :
                       case( EntityRank::BSPLINE_2 ) :
                       case( EntityRank::BSPLINE_3 ) :
                       {
                          this->map_node_to_bspline_same_mesh(
                                  tSourceIndex,
                                  tTargetIndex,
                                  aTargetEntityRank );


                           break;
                       }
                       default :
                       {
                           MORIS_ERROR( false, "perform_mapping(): aTargetEntityRank not supported.");
                           break;
                       }
                   }
                   break;
                }
                case( EntityRank::BSPLINE_1 ) :
                case( EntityRank::BSPLINE_2 ) :
                case( EntityRank::BSPLINE_3 ) :
                {
                    switch( aTargetEntityRank )
                    {
                        case( EntityRank::NODE ) :
                        {
                            this->map_bspline_to_node_same_mesh(
                                    tSourceIndex,
                                    aSourceEntityRank,
                                    tTargetIndex );
                            break;
                        }
                        default :
                        {
                            MORIS_ERROR( false, "perform_mapping(): aTargetEntityRank not supported.");
                            break;
                        }
                    }
                    break;
                }
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                default :
                {
                    MORIS_ERROR( false, "perform_mapping(): aSourceEntityRank not supported.");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::map_node_to_bspline_same_mesh(
                const moris_index   aSourceIndex,
                const moris_index   aTargetIndex,
                const EntityRank    aBSplineRank )
        {

            // create the model if it has not been created yet
            this->create_iwg_and_model();

            // get rank of B-Splines
            uint tOrder = entity_rank_to_order( aBSplineRank );

            // set order to 1
            mModel->set_dof_order( tOrder );

            // set weak bcs from field
            mModel->set_weak_bcs_from_nodal_field( aSourceIndex );

            // test if output mesh is HMR
            if( mTargetMesh->get_mesh_type() == MeshType::HMR )
            {
                // perform L2 projection
                mModel->solve( mTargetMesh->get_field( aTargetIndex, aBSplineRank ) );
            }
            else
            {
                // create vector with solution
                Matrix< DDRMat > tSolution;

                // perform L2 projection
                mModel->solve( tSolution );

                // get number of coeffs of
                uint tNumberOfCoeffs = mTargetMesh->get_num_coeffs( tOrder );

                // make sure that solution is correct
                MORIS_ERROR( tNumberOfCoeffs == tSolution.length(),
                                  "perform_mapping() number of coeffs does not match" );

                // copy solution into target
                for( uint k=0; k<tNumberOfCoeffs; ++k )
                {
                    // get ref to value
                    mTargetMesh->get_value_of_scalar_field(
                            aTargetIndex,
                            aBSplineRank,
                            k ) = tSolution( k );
                }

            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::map_bspline_to_node_same_mesh(
                         const moris_index   aSourceIndex,
                         const EntityRank    aBSplineRank,
                         const moris_index   aTargetIndex )
        {

            // get number of nodes
            moris_index tNumberOfNodes = mTargetMesh->get_num_nodes();

            // loop over all nodes
            for( moris_index k=0;  k<tNumberOfNodes; ++k )
            {
                // get weights
                const Matrix< DDRMat > & tTMatrix = mTargetMesh->
                        get_t_matrix_of_node_loc_ind( k, aBSplineRank );

                // get indices
                Matrix< IndexMat > tBSplines = mTargetMesh->
                        get_bspline_inds_of_node_loc_ind( k, aBSplineRank );

                // get number of coefficients
                uint tNumberOfCoeffs = tTMatrix.length();

                // value of node
                real & tValue = mTargetMesh->get_value_of_scalar_field(
                        aTargetIndex,
                        EntityRank::NODE,
                        k );

                // reset value
                tValue = 0.0;

                for( uint i=0; i<tNumberOfCoeffs; ++i )
                {
                    tValue +=  mSourceMesh->get_value_of_scalar_field(
                            aSourceIndex,
                            aBSplineRank,
                            tBSplines( i ) ) * tTMatrix( i );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Mapper::map_node_to_element_same_mesh(
                         const moris_index   aSourceIndex,
                         const moris_index   aTargetIndex )
        {
            // create the model if it has not been created yet
            this->create_iwg_and_model();

            // set weak bcs from field
            mModel->set_weak_bcs_from_nodal_field( aSourceIndex );

            // get number of elements
            uint tNumberOfElements = mTargetMesh->get_num_elems();

            // loop over all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                // get ref to entry in database
                real & tValue = mTargetMesh->get_value_of_scalar_field(
                        aTargetIndex,
                        EntityRank::ELEMENT,
                        e );

                // calculate value
                tValue = mModel->compute_element_average( e );
            }
        }

 //------------------------------------------------------------------------------

        void
        Mapper::create_nodes_for_filter()
        {
            if( ! mHaveNodes )
            {
                // get number of nodes from mesh
                uint tNumberOfNodes = mSourceMesh->get_num_nodes();

                // reserve node container
                mNodes.resize( tNumberOfNodes, nullptr );

                // populate container
                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    mNodes( k ) = new Node( &mSourceMesh->get_mtk_vertex( k ) );
                }

                // link to neighbors
                /*for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    Matrix< IndexMat > tNodeIndices =
                            mSourceMesh->get_entity_connected_to_entity_loc_inds(
                                    k,
                                    EntityRank::NODE,
                                    EntityRank::NODE );

                    uint tNumberOfConnectedNodes = tNodeIndices.length();
                    mNodes( k )->init_neighbor_container( tNumberOfConnectedNodes );

                    for( uint i=0; i<tNumberOfConnectedNodes; ++i )
                    {
                        mNodes( k )->insert_neighbor( mNodes( tNodeIndices( i ) ) );
                    }

                } */

                // set node flag
                mHaveNodes = true;
            }
        }
//------------------------------------------------------------------------------

        void
        Mapper::perform_filter(
                        const std::string & aSourceLabel,
                        const real        & aFilterRadius,
                        Matrix< DDRMat >  & aValues )
        {

            MORIS_ERROR( par_size() == 1,
                    "The filter is not written for parallel. In order do use it, mtk::Mapper needs access to node information from the aura.");

            // fixme: the following two lines only work for HMR
            moris_index tFieldIndex
                = mSourceMesh->get_field_ind(
                        aSourceLabel,
                        EntityRank::NODE );

            const Matrix< DDRMat > & tSourceField =
                    mSourceMesh->get_field( tFieldIndex, EntityRank::NODE );

            // calculate weights if this was not done already
            this->calculate_filter_weights( aFilterRadius );

            // get number of nodes on target
            uint tNumberOfNodes = mNodes.size();

            aValues.set_size( tNumberOfNodes, 1 );

            for( uint k=0; k<tNumberOfNodes; ++k )
            {

                Matrix< IndexMat > & tIndices = mNodes( k )->get_node_indices();

                uint tNumberOfIndices = tIndices.length();

                Matrix< DDRMat > tValues( tNumberOfIndices , 1 );

                for( uint i=0; i<tNumberOfIndices; ++i )
                {
                    tValues( i ) = tSourceField( tIndices( i ) );
                }

                // fill vector with values
                aValues( k ) = dot ( mNodes( k )->get_weights(), tValues );
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::calculate_filter_weights( const real & aFilterRadius )
        {
            if( mFilterRadius != aFilterRadius )
            {
                // remember radius
                mFilterRadius = aFilterRadius;

                // create nodes for the filter
                this->create_nodes_for_filter();

                for( Node * tNode : mNodes )
                {

                    // flag myself
                    tNode->flag();

                    // cell containing neighbors
                    Cell< Node * > tNeighbors;

                    tNode->get_nodes_in_proximity( tNode->get_coords(), aFilterRadius, tNeighbors );

                    uint tNumberOfNeighbors = tNeighbors.size();

                    Matrix< DDRMat > & tWeights = tNode->get_weights();
                    tWeights.set_size( tNumberOfNeighbors, 1 );

                    Matrix< IndexMat > & tIndices = tNode->get_node_indices();
                    tIndices.set_size( tNumberOfNeighbors, 1 );

                    real tMyLevel = tNode->get_level();

                    uint tCount = 0;
                    for( Node * tNeighbor : tNeighbors )
                    {
                        // Kurt's formula with level based average
                        tWeights( tCount )   =
                                ( aFilterRadius - tNeighbor->get_distance() )
                               *  ( tMyLevel + 1.0 ) / ( ( real ) tNeighbor->get_level() + 1.0);

                        // Simple Weight by distance
                        //tWeights( tCount )   =
                        //        ( aFilterRadius - tNeighbor->get_distance() );

                        // save index
                        tIndices( tCount++ ) = tNeighbor->get_index();

                        // unflag neighbors
                        tNeighbor->unflag();
                    }

                    tWeights = tWeights / sum( tWeights );

                    // unflag this node
                    tNode->unflag();
                }
            }

        }
    } /* namespace mtk */
} /* namespace moris */
