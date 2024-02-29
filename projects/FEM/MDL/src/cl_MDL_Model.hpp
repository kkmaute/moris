/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_Model.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_

#include "moris_typedefs.hpp"                       //MRS/COR/src
#include "cl_Vector.hpp"                            //MRS/CNT/src

#include "cl_Map.hpp"                        //MRS/CNT/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    class Library_IO;
    //------------------------------------------------------------------------------
    namespace mtk
    {
        class Mesh_Manager;
        class Field;
    }

    namespace vis
    {
        class Output_Manager;
    }

    namespace fem
    {
        class IWG;
        class Node_Base;
        class Cell;
        class Set;
        class Field_Interpolator;
        class Set_User_Info;
    }

    namespace dla
    {
        class Linear_Solver;
        class Linear_Solver_Algorithm;
    }

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Equation_Model;
        class Design_Variable_Interface;
        enum class Dof_Type;
    }
    namespace tsa
    {
        class Time_Solver;
        class Time_Solver_Algorithm;
    }
    namespace sol
    {
        class SOL_Warehouse;
    }
    namespace mdl
    {
        //------------------------------------------------------------------------------

        class Model
        {
            private:
                // pointer to reference mesh
                std::shared_ptr< mtk::Mesh_Manager > mMeshManager = nullptr;
                moris_index        mMeshPairIndex;

                // pointer to equation model
                std::shared_ptr< MSI::Equation_Model > mEquationModel = nullptr;

                // by default, this value is set to the order of the
                // Lagrange modes
                moris::uint mBSplineIndex = 0;

                // pointer to model solver interface
                MSI::Model_Solver_Interface * mModelSolverInterface = nullptr;

                // pointer to solver interface
                MSI::MSI_Solver_Interface   * mSolverInterface = nullptr;

                // pointer to solver warehouse
                std::shared_ptr< sol::SOL_Warehouse > mSolverWarehouse = nullptr;

                Matrix< DDUMat > mAdofMap;

                // pointer to output manager
                vis::Output_Manager * mOutputManager = nullptr;
                bool mOutputManagerOwned = false;

                // bool for multigrid use
                bool mUseMultigrid = false;

                // bool whether new Ghost sets are being used
                bool mUseNewGhostSets = false;

                // pointer to library for input reading
                std::shared_ptr< Library_IO > mLibrary = nullptr;

                MSI::Design_Variable_Interface * mDesignVariableInterface = nullptr;

                Vector< Vector< ParameterList > > mFEMParameterList;
                Vector< Vector< ParameterList > > mMSIParameterList;
                Vector< Vector< ParameterList > > mSOLParameterList;
                Vector< Vector< ParameterList > > mVISParameterList;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ] aMeshManager   pointer to mesh info
                 * @param[ in ] aBSplineIndex  ???
                 * @param[ in ] aSetInfo       list of set user info
                 * @param[ in ] aMeshPairIndex ???
                 * @param[ in ] aUseMultigrid  bool for multigrid use
                 */
                Model(
                        std::shared_ptr< mtk::Mesh_Manager > aMeshManager,
                        const uint                           aBSplineIndex,
                        Vector< fem::Set_User_Info >  & aSetInfo,
                        const moris_index                    aMeshPairIndex = 0,
                        const bool                           aUseMultigrid  = false );

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ] aMeshManager   pointer to mesh info
                 * @param[ in ] aBSplineIndex  ???
                 * @param[ in ] aSetInfo       cell of set user info
                 * @param[ in ] aDVInterface   pointer to dv interface
                 * @param[ in ] aMeshPairIndex ???
                 * @param[ in ] aUseMultigrid  bool for multigrid use
                 */
                Model(
                        std::shared_ptr< mtk::Mesh_Manager >aMeshManager,
                        const uint                          aBSplineIndex,
                        Vector< fem::Set_User_Info > & aSetInfo,
                        MSI::Design_Variable_Interface    * aDesignVariableInterface,
                        const moris_index                   aMeshPairIndex = 0,
                        const bool                          aUseMultigrid  = false );

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ] aBSplineIndex  ???
                 * @param[ in ] aMeshPairIndex ???
                 */
                Model(
                        std::shared_ptr< Library_IO > aLibrary,
                        const uint                    aBSplineIndex,
                        const moris_index             aMeshPairIndex = 0 );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Model();

                //------------------------------------------------------------------------------

                void free_memory();

                //------------------------------------------------------------------------------
                /**
                 * solve
                 */
                void perform( uint aIndex = 0 );

                //------------------------------------------------------------------------------

                void perform_forward_analysis();

                //------------------------------------------------------------------------------

                Vector< moris::Matrix< DDRMat > > get_IQI_values();

                //------------------------------------------------------------------------------

                void perform_sensitivity_analysis();

                //------------------------------------------------------------------------------

                void set_solver_warehouse_hack( std::shared_ptr< sol::SOL_Warehouse > aSolverWarehouse)
                {
                    mSolverWarehouse = aSolverWarehouse;
                };

                //------------------------------------------------------------------------------
                /**
                 * set MTK performer
                 * @param[ in ] aMTKPerformer the MTK mesh manager
                 */
                void set_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

                //------------------------------------------------------------------------------

                void set_design_variable_interface( MSI::Design_Variable_Interface * aDesignVariableInterface );

                //------------------------------------------------------------------------------
                /**
                 * initialize the set - build FEM, MSI, VIS and SOL
                 */
                void initialize();

                //------------------------------------------------------------------------------
                /**
                 * set dof order
                 * @param[ in ] aOrder an order
                 */
                void set_dof_order( const uint aOrder );

                //------------------------------------------------------------------------------
                /**

                 * get the model solver interface
                 */
                MSI::Model_Solver_Interface * get_model_solver_interface()
                {
                    return mModelSolverInterface;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the solver interface
                 */
                MSI::MSI_Solver_Interface * get_solver_interface()
                {
                    return mSolverInterface;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the fem model
                 * (used for UT)
                 */
                std::shared_ptr< MSI::Equation_Model > get_fem_model()
                {
                    return mEquationModel;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the adof map
                 * @param[ out ] adof map
                 */
                Matrix< DDUMat > & get_adof_map()
                {
                    return mAdofMap;
                }

                //------------------------------------------------------------------------------
                /**
                 * FIXME
                 * sets the weak BCs
                 * @param[ in ] aWeakBCs matrix with weak BCs value
                 */
                void set_weak_bcs( const Matrix< DDRMat > & aWeakBCs );

                //------------------------------------------------------------------------------
                /**
                 * FIXME
                 * sets the weak BCs from a nodal field
                 * @param[ in ] aFieldIndex an index for the field
                 */
                void set_weak_bcs_from_nodal_field( moris_index aFieldIndex );

                //------------------------------------------------------------------------------
                /*
                 * set output manager
                 * @param[ in ] aOutputManager pointer to output manager
                 */
                void set_output_manager( vis::Output_Manager * aOutputManager )
                {
                    mOutputManager = aOutputManager;
                }

                //------------------------------------------------------------------------------
                /*
                 * output solution
                 * @param[ in ] aVisMeshIndex mesh index on with solution is displayed
                 * @param[ in ] aTime         a time at which solution is displayed
                 * @param[ in ] aCloseFile    boolean indicating the closing of the exodus file
                 */
                void output_solution(
                        const uint aVisMeshIndex,
                        const real aTime,
                        const bool aCloseFile );

                //------------------------------------------------------------------------------

                /**
                 * return fields
                 */
                Vector< std::shared_ptr< mtk::Field > > get_mtk_fields();

                //------------------------------------------------------------------------------

                /**
                 * @brief builds an unordered map from the MSI parameter list relating which dof types correspond to which B-spline mesh indices
                 * 
                 * @return std::unordered_map< MSI::Dof_Type, moris_index > 
                 */
                std::unordered_map< enum MSI::Dof_Type, moris_index >
                build_dof_type_to_mesh_index();

                //------------------------------------------------------------------------------

                /**
                 * @brief Turn on usage of new ghost sets. Needed to trigger the FEM model to use the correct sets
                 * 
                 * @param aUseNewGhost 
                 */
                void
                set_use_new_ghost_mesh_sets( bool aUseNewGhost )
                {
                    mUseNewGhostSets = aUseNewGhost;
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief Get the spatial dim
                 *
                 * @return uint 2,3
                 */

                uint get_spatial_dim();
        };
        //------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_MDL_MODEL_HPP_ */

