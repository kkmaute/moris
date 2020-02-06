/*
 * cl_MSI_Solver_Interface.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_

#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Solver_Interface.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
class Dist_Vector;

namespace mtk
{
    class Mesh;
}
namespace mdl
{
    class Model;
}

    namespace MSI
    {
        class MSI_Solver_Interface : public moris::Solver_Interface
        {
        private:
                moris::MSI::Model_Solver_Interface * mMSI = nullptr;
                moris::MSI::Dof_Manager            * mDofMgn = nullptr;

                Dist_Vector                        * mSolutionVector = nullptr;
                Dist_Vector                        * mPrevSolutionVector = nullptr;
                Dist_Vector                        * mExactSolFromFile = nullptr;
                Matrix< DDRMat>  mTime;

                moris::Cell< enum MSI::Dof_Type > mListOfDofTypes;
                Cell< moris::Cell< enum MSI::Dof_Type > > mListOfSecundaryDofTypes;

                mdl::Model * mModel = nullptr;

                bool mIsForward = true;

        public:
            MSI_Solver_Interface( )
            {
                mTime = { {0.0}, {1.0} };
            };

//------------------------------------------------------------------------------
            MSI_Solver_Interface( moris::MSI::Model_Solver_Interface * aMSI ) : mMSI( aMSI ),
                                                                                mDofMgn( mMSI->get_dof_manager() )
            {
                mTime = { {0.0}, {1.0} };

                mMSI->set_solver_interface( this );
            };

//------------------------------------------------------------------------------
            ~MSI_Solver_Interface() {};

//------------------------------------------------------------------------------

            void set_solution_vector( Dist_Vector * aSolutionVector )
            {
                mSolutionVector = aSolutionVector;
            }

//------------------------------------------------------------------------------

            void set_model( mdl::Model * aModel )
            {
                mModel = aModel;
            }

//------------------------------------------------------------------------------

            void set_is_forward( bool aIsForward )
            {
            	mIsForward = aIsForward;
            }


//------------------------------------------------------------------------------

            void get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename );

//------------------------------------------------------------------------------

            void get_residual_vector_for_output( const char* aFilename );

//------------------------------------------------------------------------------

            void write_solution_to_hdf5_file( const char* aFilename );

//------------------------------------------------------------------------------

            void set_solution_vector_prev_time_step( Dist_Vector * aSolutionVector )
            {
                mPrevSolutionVector = aSolutionVector;
            }
//------------------------------------------------------------------------------

            void set_time( const Matrix< DDRMat> & aTime )
            {
                mTime = aTime;
            };

//------------------------------------------------------------------------------

            void perform_mapping( )
            {
            };

//------------------------------------------------------------------------------

            void free_block_memory( const uint aBlockInd )
            {
                mMSI->get_eqn_block( aBlockInd )->free_matrix_memory();
            };

//------------------------------------------------------------------------------

            void initialize_block( const uint aBlockInd,
                                   const bool aIsResidual )
            {
                mMSI->get_eqn_block( aBlockInd )->initialize_set( aIsResidual, mIsForward );
            };

//------------------------------------------------------------------------------

            void initiate_output( const uint aOutputIndex,
                                  const uint aTime );

//------------------------------------------------------------------------------

            void set_requested_dof_types( const moris::Cell< enum MSI::Dof_Type > aListOfDofTypes )
            {
               mListOfDofTypes = aListOfDofTypes;
            };

//------------------------------------------------------------------------------

            void set_secundary_dof_types( const Cell< moris::Cell< enum MSI::Dof_Type > > aListOfDofTypes )
            {
               mListOfSecundaryDofTypes = aListOfDofTypes;
            };

//------------------------------------------------------------------------------

            moris::Cell< enum MSI::Dof_Type > get_requested_dof_types()
            {
                return mListOfDofTypes;
            };

//------------------------------------------------------------------------------

            moris::Cell< moris::Cell< enum MSI::Dof_Type > > get_secundary_dof_types()
            {
                return mListOfSecundaryDofTypes;
            };

//------------------------------------------------------------------------------

            void set_time_levels_for_type( const enum Dof_Type aDofType,
                                           const moris::uint   aNumTimeLevels )
            {
                mDofMgn->set_time_levels_for_type( aDofType, aNumTimeLevels );
            };

//------------------------------------------------------------------------------

            // number of elements blocks on proc
            moris::uint get_num_my_blocks()
            {
                return mMSI->get_num_eqn_blocks();
            };

//------------------------------------------------------------------------------

            // number of elements on proc
            moris::uint get_num_my_elements()
            {
                return mMSI->get_num_eqn_objs();
            };

//------------------------------------------------------------------------------

            moris::uint get_num_my_elements_on_block( uint aBlockInd )
            {
                return mMSI->get_eqn_block( aBlockInd )->get_num_equation_objects();
            };

//------------------------------------------------------------------------------

             moris::uint get_num_my_dofs()
             {
                 return mDofMgn->get_num_owned_adofs();
             };

//------------------------------------------------------------------------------

             moris::uint get_max_num_global_dofs()
             {
                 moris::uint tNumMyDofs        = mDofMgn->get_num_owned_adofs();
                 moris::uint tMaxNumGlobalDofs = mDofMgn->get_num_owned_adofs();

                 // sum up all distributed dofs
                 sum_all( tNumMyDofs, tMaxNumGlobalDofs );

                 return tMaxNumGlobalDofs;
             };

//------------------------------------------------------------------------------
             // local-to-global map
             moris::Matrix< DDSMat > get_my_local_global_map()
             {
                 return mDofMgn->get_local_adof_ids();
             };

//------------------------------------------------------------------------------

             moris::Matrix< DDSMat > get_my_local_global_map( const moris::Cell< enum Dof_Type > & aListOfDofTypes )
             {
                 return mDofMgn->get_local_adof_ids( aListOfDofTypes );
             };

//------------------------------------------------------------------------------
             Matrix< DDSMat > get_my_local_global_overlapping_map( )
             {
                 return mDofMgn->get_local_overlapping_adof_ids();
             };

//------------------------------------------------------------------------------

             Matrix< DDSMat > get_my_local_global_overlapping_map( const moris::Cell< enum Dof_Type > & aListOfDofTypes )
             {
                 return mDofMgn->get_local_overlapping_adof_ids( aListOfDofTypes );
             };

//------------------------------------------------------------------------------

             void get_element_matrix( const moris::uint      & aMyElementInd,
                                            Matrix< DDRMat > & aElementMatrix )
             {
                 mMSI->get_eqn_obj( aMyElementInd )->set_time( mTime );
                 mMSI->get_eqn_obj( aMyElementInd )->get_egn_obj_jacobian( aElementMatrix, mSolutionVector );
             };

//------------------------------------------------------------------------------

             void get_element_matrix( const moris::uint      & aMyBlockInd,
                                      const moris::uint      & aMyElementInd,
                                            Matrix< DDRMat > & aElementMatrix )
             {
                 mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )->set_time( mTime );
                 mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )->get_egn_obj_jacobian( aElementMatrix, mSolutionVector );
             };

//------------------------------------------------------------------------------
             void  get_element_topology( const moris::uint      & aMyElementInd,
                                               Matrix< DDSMat > & aElementTopology )
             {
                 mMSI->get_eqn_obj( aMyElementInd )->get_equation_obj_dof_ids( aElementTopology );
             };

//------------------------------------------------------------------------------
             void  get_element_topology( const moris::uint      & aMyBlockInd,
                                         const moris::uint      & aMyElementInd,
                                               Matrix< DDSMat > & aElementTopology )
             {
                 mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )->get_equation_obj_dof_ids( aElementTopology );
//                 mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )
//                                                   ->get_equation_obj_dof_ids( aElementTopology, mListOfDofTypes, mDofMgn );
             };

//------------------------------------------------------------------------------

             Matrix< DDUMat > get_constr_dof()
             {
                 // Matrix< DDUMat > tLocalConstrIds;// = mDofMgn->get_full_to_free_constraints();
                 return Matrix< DDUMat >(0,0);
             };

//------------------------------------------------------------------------------

             void get_equation_object_rhs( const moris::uint              & aMyElementInd,
                                         Cell< Matrix< DDRMat > > & aElementRHS )
             {
                     mMSI->get_eqn_obj( aMyElementInd )->set_time( mTime );
                     mMSI->get_eqn_obj( aMyElementInd )->get_equation_obj_residual( aElementRHS, mSolutionVector  );
             };

//------------------------------------------------------------------------------

             void get_equation_object_rhs( const moris::uint              & aMyBlockInd,
                                   const moris::uint              & aMyElementInd,
                                         Cell< Matrix< DDRMat > > & aElementRHS )
             {

                     mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )->set_time( mTime );
                     mMSI->get_eqn_block( aMyBlockInd )->get_equation_object_list()( aMyElementInd )->get_equation_obj_residual( aElementRHS, mSolutionVector  );
             };

//------------------------------------------------------------------------------

             mtk::Mesh * get_mesh_pointer_for_multigrid( )
             {
                 return mMSI->get_mesh_pointer_for_multigrid();
             };

//------------------------------------------------------------------------------

             void read_multigrid_maps( const moris::uint               aLevel,
                                       const moris::Matrix< DDSMat > & aExtFineIndices,
                                       const moris::sint               aTypeTimeIdentifier,
                                             moris::Matrix< DDSMat > & aInternalFineIndices)
             {
                 mMSI->get_msi_multigrid_pointer()->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
             };

//------------------------------------------------------------------------------

             const moris::Cell< Matrix< DDSMat > > & get_lists_of_multigrid_identifiers()
             {
                 return mMSI->get_msi_multigrid_pointer()->get_lists_of_multigrid_identifiers();
             };

//------------------------------------------------------------------------------

             const moris::Cell< Matrix< DDUMat > > & get_lists_of_ext_index_multigrid()
             {
                 return mMSI->get_msi_multigrid_pointer()->get_lists_of_ext_index_multigrid();
             };

//------------------------------------------------------------------------------

             const moris::Cell< moris::Cell< Matrix< DDSMat > > > & get_multigrid_map( )
             {
                 return mMSI->get_msi_multigrid_pointer()->get_multigrid_map();
             };

//------------------------------------------------------------------------------

             const moris::Matrix< DDUMat > & get_number_remaining_dofs()
             {
                 return mMSI->get_msi_multigrid_pointer()->get_number_remaining_dofs();
             };

//------------------------------------------------------------------------------

             const Matrix< DDSMat > & get_type_time_identifier_to_type_map()
             {
                 return mDofMgn->get_typetime_identifier_to_type_map();
             };

//------------------------------------------------------------------------------

             moris::sint get_adof_index_for_type( moris::uint aDofType )
             {
                 return mMSI->get_adof_index_for_type( aDofType );;
             };
        };
    }
}

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
