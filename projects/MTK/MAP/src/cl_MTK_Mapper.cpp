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

#include "cl_MDL_Model.hpp"

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------

        Mapper::Mapper( mtk::Mesh * aMesh ) :
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
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_iwg_and_model()
        {
            if( ! mHaveIwgAndModel )
            {
                // create IWG object
                mIWG = new moris::fem::IWG_L2();

                // create model
                mModel = new mdl::Model( mTargetMesh, mIWG );

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
    } /* namespace mtk */
} /* namespace moris */
