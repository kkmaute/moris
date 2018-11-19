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
        class Model_Solver_Interface
        {
        private:
            moris::Cell< Equation_Object* > & mEquationObjectList;
            Dof_Manager                       mDofMgn;

            mtk::Mesh                       * mMesh;

        public:
        /**
         * @brief Model solver interface constructor. This function is tested by the test [MSI_Test][MSI_Test_parallel]
         *
         * @param[in] aListEqnObj   List containing all the equation objects.
         * @param[in] aCommTable    Communication table for adofs.
         *
         */
        Model_Solver_Interface(      moris::Cell < Equation_Object* >                  & aListEqnObj,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs ) : mEquationObjectList( aListEqnObj ),
                                                                                                          mDofMgn( aCommTable )
        {
            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );
            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( aListEqnObj );

            mDofMgn.initialize_pdof_host_list( aListEqnObj );

            mDofMgn.create_adofs();

            mDofMgn.set_pdof_t_matrix();

            for ( Equation_Object* tElement : mEquationObjectList )
            {
                tElement->create_my_pdof_list();

                tElement->create_my_list_of_adof_ids();

                tElement->set_unique_adof_map();
            }

            //Multigrid tMultigrid( &mDofMgn );
        };

        Model_Solver_Interface(      moris::Cell < Equation_Object* >                  & aListEqnObj,
                               const Matrix< IdMat >                                   & aCommTable,
                               const moris::map< moris::moris_id, moris::moris_index > & aAdofLocaltoGlobalMap,
                               const moris::uint                                         aNumMaxAdofs,
                                     mtk::Mesh                                         * aMesh ) : mEquationObjectList( aListEqnObj ),
                                                                                                   mDofMgn( aCommTable ),
                                                                                                   mMesh( aMesh )
        {
            mDofMgn.set_adof_map( & aAdofLocaltoGlobalMap );
            mDofMgn.set_max_num_adofs( aNumMaxAdofs );

            mDofMgn.initialize_pdof_type_list( aListEqnObj );

            mDofMgn.initialize_pdof_host_list( aListEqnObj );

            mDofMgn.create_adofs();

            mDofMgn.set_pdof_t_matrix();

            for ( Equation_Object* tElement : mEquationObjectList )
            {
                tElement->create_my_pdof_list();

                tElement->create_my_list_of_adof_ids();

                tElement->set_unique_adof_map();
            }

            Multigrid tMultigrid( this, mMesh );
        };

//------------------------------------------------------------------------------
        ~Model_Solver_Interface(){};

//------------------------------------------------------------------------------
        moris::uint get_num_eqn_objs()
        {
            return mEquationObjectList.size();
        };

//------------------------------------------------------------------------------
        Dof_Manager * get_dof_manager()
        {
            return & mDofMgn;
        };

//------------------------------------------------------------------------------

        void
        get_equation_obj_jacobian( const moris::uint      & aEqnObjInd,
                                         Matrix< DDRMat > & aEqnObjMatrix,
                                         Dist_Vector      * aSolutionVector)
        {
            mEquationObjectList( aEqnObjInd )->get_egn_obj_jacobian( aEqnObjMatrix, aSolutionVector );
        };

//------------------------------------------------------------------------------
        void get_equation_obj_residual( const moris::uint      & aEqnObjInd,
                                              Matrix< DDRMat > & aEqnObjRHS,
                                              Dist_Vector      * aSolutionVector )
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_residual( aEqnObjRHS, aSolutionVector  );
        };

//------------------------------------------------------------------------------
        void get_equation_obj_dof_ids( const moris::uint      & aEqnObjInd,
                                             Matrix< DDSMat > & aElementTopology )
        {
            mEquationObjectList( aEqnObjInd )->get_equation_obj_dof_ids( aElementTopology );
        };

//------------------------------------------------------------------------------
    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_MANAGER_HPP_ */
