/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Model_Solver_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_EQUATION_MANAGER_HPP_
#define SRC_FEM_CL_EQUATION_MANAGER_HPP_

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"

#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Multigrid.hpp"
#include "cl_MSI_Equation_Set.hpp"

#include "cl_Parameter_List.hpp"

namespace moris
{
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
            Parameter_List mMSIParameterList;

            //! List of equation blocks
            Vector< MSI::Equation_Set* > mEquationSets;

            //! List of equation objects
            Vector< Equation_Object* > mEquationObjectList;

            //! Dof manager object
            Dof_Manager mDofMgn;

            //! Pointer to Mesh. Only used in multigrid
            mtk::Mesh* mMesh = nullptr;

            //! Multigrid object pointer
            Multigrid* mMultigrid = nullptr;

            MSI::MSI_Solver_Interface* mSolverInterface = nullptr;

            std::shared_ptr< MSI::Equation_Model > mEquationModel = nullptr;

            friend class MSI_Solver_Interface;
            friend class Multigrid;

            //------------------------------------------------------------------------------

            void msi_checker();

            //------------------------------------------------------------------------------

            void pdof_host_checker();

            //------------------------------------------------------------------------------

          public:
            Model_Solver_Interface() {};

            //------------------------------------------------------------------------------

            Model_Solver_Interface(
                    const Parameter_List&       aMSIParameterList,
                    Vector< Equation_Object* >& aListEqnObj )
                    : mMSIParameterList( aMSIParameterList )
                    , mEquationObjectList( aListEqnObj ) {};

            //------------------------------------------------------------------------------

            /**
             * @brief Model solver interface constructor. This function is tested by the test [MSI_Test][MSI_Test_parallel]
             *
             * @param[in] aListEqnObj   List containing all the equation objects.
             * @param[in] aCommTable    Communication table for adofs.
             *
             */
            Model_Solver_Interface(
                    const Parameter_List&                                    aMSIParameterList,
                    Vector< MSI::Equation_Set* >&                            aElementBlocks,
                    const Matrix< IdMat >&                                   aCommTable,
                    const moris::map< moris::moris_id, moris::moris_index >& aAdofLocaltoGlobalMap,
                    const moris::uint                                        aNumMaxAdofs )
                    : mMSIParameterList( aMSIParameterList )
                    , mEquationSets( aElementBlocks )
                    , mDofMgn( aCommTable, this )
            {
                this->create_equation_object_list();

                mDofMgn.set_adof_map( &aAdofLocaltoGlobalMap );

                mDofMgn.set_max_num_adofs( aNumMaxAdofs );

                mDofMgn.initialize_pdof_type_list( mEquationSets );
            };

            //------------------------------------------------------------------------------

            Model_Solver_Interface(
                    const Parameter_List&                                    aMSIParameterList,
                    Vector< MSI::Equation_Set* >&                            aElementBlocks,
                    const Matrix< IdMat >&                                   aCommTable,
                    const moris::map< moris::moris_id, moris::moris_index >& aAdofLocaltoGlobalMap,
                    const moris::uint                                        aNumMaxAdofs,
                    mtk::Mesh*                                               aMesh )
                    : mMSIParameterList( aMSIParameterList )
                    , mEquationSets( aElementBlocks )
                    , mDofMgn( aCommTable, this )
                    , mMesh( aMesh )
            {
                this->create_equation_object_list();

                mDofMgn.set_adof_map( &aAdofLocaltoGlobalMap );

                mDofMgn.set_max_num_adofs( aNumMaxAdofs );

                mDofMgn.initialize_pdof_type_list( mEquationSets );
            };

            //------------------------------------------------------------------------------

            Model_Solver_Interface(
                    const Parameter_List&                         aMSIParameterList,
                    const std::shared_ptr< MSI::Equation_Model >& aEquationModel,
                    mtk::Mesh*                                    aMesh );

            //------------------------------------------------------------------------------

            ~Model_Solver_Interface()
            {
                if ( mMultigrid != nullptr )
                {
                    delete mMultigrid;
                }
            };

            //------------------------------------------------------------------------------

            void
            create_equation_object_list();

            //------------------------------------------------------------------------------

            void finalize();

            //------------------------------------------------------------------------------

            moris::uint
            get_num_equation_sets()
            {
                return mEquationSets.size();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_eqn_objs()
            {
                return mEquationObjectList.size();
            };

            //------------------------------------------------------------------------------

            moris::uint
            get_num_eqn_objs_on_block( moris::uint aBlockInd )
            {
                return mEquationSets( aBlockInd )->get_num_equation_objects();
            };

            //------------------------------------------------------------------------------

            Equation_Set*
            get_equation_set( const moris::uint& aMyEquSetInd )
            {
                return mEquationSets( aMyEquSetInd );
            };

            //------------------------------------------------------------------------------

            std::shared_ptr< MSI::Equation_Model >
            get_equation_model()
            {
                return mEquationModel;
            };
            //------------------------------------------------------------------------------

            Equation_Object*
            get_eqn_obj( const moris::uint& aMyEquObjInd )
            {
                return mEquationObjectList( aMyEquObjInd );
            };

            //------------------------------------------------------------------------------

            Dof_Manager*
            get_dof_manager()
            {
                return &mDofMgn;
            };

            //------------------------------------------------------------------------------

            Multigrid*
            get_msi_multigrid_pointer()
            {
                return mMultigrid;
            };

            //------------------------------------------------------------------------------

            void
            read_multigrid_maps(
                    const moris::uint              aLevel,
                    const moris::Matrix< DDSMat >& aExtFineIndices,
                    const moris::sint              aTypeTimeIdentifier,
                    moris::Matrix< DDSMat >&       aInternalFineIndices )
            {
                mMultigrid->read_multigrid_maps( aLevel, aExtFineIndices, aTypeTimeIdentifier, aInternalFineIndices );
            };

            //------------------------------------------------------------------------------

            mtk::Mesh*
            get_mesh_pointer_for_multigrid()
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
            moris_index get_adof_index_for_type( MSI::Dof_Type aDofType );

            //------------------------------------------------------------------------------

            /**
             * @brief determines the maximal underlying adof mesh index based on the used dof types.
             *  Is not limited to underlying BSpline meshes.
             */
            moris::sint get_max_adof_index();

            //------------------------------------------------------------------------------

            moris::uint
            get_time_levels_for_type( const enum Dof_Type aDofType )
            {
                return mDofMgn.get_time_levels_for_type( aDofType );
            };

            //------------------------------------------------------------------------------

            void
            set_solver_interface( MSI::MSI_Solver_Interface* aSolverInterface )
            {
                mSolverInterface = aSolverInterface;
            };

            //------------------------------------------------------------------------------

            MSI::MSI_Solver_Interface*
            get_solver_interface()
            {
                return mSolverInterface;
            };

            //------------------------------------------------------------------------------

            Parameter_List&
            get_msi_parameterlist()
            {
                return mMSIParameterList;
            };

            //-------------------------------------------------------------------------------------------------------

            moris::uint
            get_num_eigen_vectors()
            {
                return mMSIParameterList.get< int >( "number_eigen_vectors" );
            }

            //------------------------------------------------------------------------------
        };
    }    // namespace MSI
}    // namespace moris

#endif /* SRC_FEM_CL_EQUATION_MANAGER_HPP_ */
