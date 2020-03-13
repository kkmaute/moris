/*
 * cl_Model_Solver_Interface.hpp
 *
 *  Created on: Jul 22, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_MANAGER_HPP_
#define SRC_FEM_CL_EQUATION_MANAGER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"

#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Multigrid.hpp"
#include "cl_MSI_Equation_Set.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    class Dist_Vector;

    namespace mtk
    {
        class Mesh;
    }
    namespace MSI
    {
        class MSI_Solver_Interface;
        class Equation_Model;
        class Model_Solver_Interface
        {
        private:

            //! Parameter list for this model solver interface
            ParameterList mMSIParameterList;

            //! List of equation blocks
            Cell< MSI::Equation_Set * >     mEquationBlocks;

            //! List of equation objects
            moris::Cell< Equation_Object * > mEquationObjectList;

            //! Dof manager object
            Dof_Manager                      mDofMgn;

            //! Pointer to Mesh. Only used in multigrid
            mtk::Mesh                      * mMesh;

            //! Multigrid object pointer
            Multigrid * mMultigrid = nullptr;

            MSI::MSI_Solver_Interface * mSolverInterface;

            friend class MSI_Solver_Interface;

//------------------------------------------------------------------------------

     public:

        Model_Solver_Interface( ){ };

//------------------------------------------------------------------------------

        Model_Solver_Interface( ParameterList                      aMSIParameterList,
                                moris::Cell < Equation_Object* > & aListEqnObj ) : mMSIParameterList( aMSIParameterList ),
                                                                                   mEquationObjectList( aListEqnObj )
        { };

//------------------------------------------------------------------------------

        /**
         * @brief Model solver interface constructor. This function is tested by the test [MSI_Test][MSI_Test_parallel]
         *
         * @param[in] aListEqnObj   List containing all the equation objects.
         * @param[in] aCommTable    Communication table for adofs.
         *
         */
        Model_Solver_Interface(      ParameterList                                       aMSIParameterList,
                                     Cell< MSI::Equation_Set * >                       & aElementBlocks,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs ) : mMSIParameterList( aMSIParameterList ),
                                                                                                          mEquationBlocks( aElementBlocks ),
                                                                                                          mDofMgn( aCommTable, this )
        {
            this->create_equation_object_list();

            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );

            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( mEquationBlocks );
        };

//------------------------------------------------------------------------------

        Model_Solver_Interface(      ParameterList                                       aMSIParameterList,
                                     Cell< MSI::Equation_Set * >                       & aElementBlocks,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs,
                                     mtk::Mesh                                         * aMesh ) : mMSIParameterList( aMSIParameterList ),
                                                                                                   mEquationBlocks( aElementBlocks ),
                                                                                                   mDofMgn( aCommTable, this ),
                                                                                                   mMesh( aMesh )
        {
            this->create_equation_object_list();

            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );

            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( mEquationBlocks );
        };

//------------------------------------------------------------------------------

        Model_Solver_Interface(      ParameterList                                       aMSIParameterList,
                                     std::shared_ptr< MSI::Equation_Model >            & aEquationModel,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs,
                                     mtk::Mesh                                         * aMesh );

//------------------------------------------------------------------------------

        ~Model_Solver_Interface()
        {
            if( mMultigrid != NULL )
            {
                delete mMultigrid;
            }
        };

//------------------------------------------------------------------------------

        void create_equation_object_list()
        {
            moris::uint tNumEquationObj = 0;

            for( luint Ik = 0; Ik < mEquationBlocks.size(); ++Ik )
            {
                tNumEquationObj = tNumEquationObj + mEquationBlocks( Ik )->get_num_equation_objects();
            }

            mEquationObjectList.reserve( tNumEquationObj );

            for( luint Ik = 0; Ik < mEquationBlocks.size(); ++Ik )
            {
                mEquationObjectList.append( mEquationBlocks( Ik )->get_equation_object_list() );
            }
        };

//------------------------------------------------------------------------------

        void finalize()
        {
            for ( luint Ik = 0; Ik < mEquationBlocks.size(); ++Ik )
            {
                mEquationBlocks( Ik )->set_model_solver_interface( this );
            }

            mDofMgn.initialize_pdof_host_list( mEquationObjectList );

            mDofMgn.create_adofs();

            mDofMgn.set_pdof_t_matrix();

            for ( Equation_Object* tElement : mEquationObjectList )
            {
                tElement->create_my_pdof_list();

                tElement->create_my_list_of_adof_ids();

                tElement->set_unique_adof_map();
            }

            if ( mMSIParameterList.get< bool >( "multigrid" ) )
            {
                mMultigrid = new Multigrid( this, mMesh );

                mMultigrid->multigrid_initialize();
            }
        };

//------------------------------------------------------------------------------

        moris::uint get_num_eqn_blocks()
        {
            return mEquationBlocks.size();
        };

//------------------------------------------------------------------------------

        moris::uint get_num_eqn_objs()
        {
            return mEquationObjectList.size();
        };

//------------------------------------------------------------------------------

        moris::uint get_num_eqn_objs_on_block( moris::uint aBlockInd )
        {
            return mEquationBlocks( aBlockInd )->get_num_equation_objects();
        };

//------------------------------------------------------------------------------

        Equation_Set * get_eqn_block( const moris::uint & aMyEquBlockInd )
        {
            return mEquationBlocks( aMyEquBlockInd );
        };

//------------------------------------------------------------------------------

        Equation_Object * get_eqn_obj( const moris::uint & aMyEquObjInd )
        {
            return mEquationObjectList( aMyEquObjInd );
        };

//------------------------------------------------------------------------------

        Dof_Manager * get_dof_manager()
        {
            return & mDofMgn;
        };

//------------------------------------------------------------------------------

        Multigrid * get_msi_multigrid_pointer( )
        {
            return mMultigrid;
        };

//------------------------------------------------------------------------------

        void read_multigrid_maps( const moris::uint               aLevel,
                                  const moris::Matrix< DDSMat > & aExtFineIndices,
                                  const moris::sint               aTypeTimeIdentifier,
                                        moris::Matrix< DDSMat > & aInternalFineIndices)
        {
            mMultigrid->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
        };

//------------------------------------------------------------------------------

        mtk::Mesh * get_mesh_pointer_for_multigrid( )
        {
            return mMesh;
        };

//------------------------------------------------------------------------------

//        boost::variant< bool, sint, real > & set_param( const std::string & aKey )
//        {
//            return mMSIParameterList( aKey );
//        }

//------------------------------------------------------------------------------

        /**
         * @brief Returns the adof index for a given dof type index.
         *
         * @param[in] aDofType    Dof type index.
         */
        moris::sint get_adof_index_for_type( moris::uint aDofType );

//------------------------------------------------------------------------------

        /**
         * @brief determines the maximal underlying adof mesh index based on the used dof types.
         *  Is not limited to underlying BSpline meshes.
         */
        moris::sint get_max_adof_index();

//------------------------------------------------------------------------------

        moris::uint get_time_levels_for_type( const enum Dof_Type aDofType )
        {
            return mDofMgn.get_time_levels_for_type( aDofType );
        };

//------------------------------------------------------------------------------

        void set_solver_interface( MSI::MSI_Solver_Interface * aSolverInterface )
        {
            mSolverInterface = aSolverInterface;
        };

//------------------------------------------------------------------------------

        MSI::MSI_Solver_Interface * get_solver_interface()
        {
            return mSolverInterface;
        };

//------------------------------------------------------------------------------
    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_MANAGER_HPP_ */
