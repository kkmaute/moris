/*
 * cl_MSI_Equation_Set.hpp
 *
 *  Created on: Apr 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_
#define SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_

#include "assert.h"
#include "cl_Communication_Tools.hpp"               //FEM/INT/src
#include "cl_Map.hpp"

namespace moris
{
namespace mtk
{
   class Cell;
}
    namespace MSI
    {
    class Model_Solver_Interface;
    class Equation_Object;
    enum class Dof_Type;
//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Equation_Set
    {
    private:


    protected:
        Cell< MSI::Equation_Object* > mEquationObjList;

        Matrix< DDRMat > mResidual;
        Matrix< DDRMat > mJacobian;

        moris::Cell< moris::Cell< Matrix< DDRMat > > > mJacobians;
        moris::Cell< Matrix< DDRMat > > mResiduals;

        // maps for the master and slave dof type
        moris::Matrix< DDSMat > mMasterDofTypeMap;
        moris::Matrix< DDSMat > mSlaveDofTypeMap;

        // map of master and slave dof types for assembly
        moris::Matrix< DDSMat > mDofAssemblyMap;  //FIXME delete
        moris::Matrix< DDSMat > mDofAssemblyMap_2;
        uint                    mTotalDof;
        moris::Matrix< DDSMat > mNumDofMap;
        Cell< moris::map< enum MSI::Dof_Type, moris::uint > > mRequestedTypeToIndexMap;

        bool mJacobianExist = false;
        bool mResidualExist = false;

        Matrix< DDRMat >mTime;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type > mEqnObjDofTypeList; // List of dof types of this equation obj

        Model_Solver_Interface * mModelSolverInterface = nullptr;

        friend class MSI::Equation_Object;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Equation_Set( ){};

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        virtual ~Equation_Set(){};

//------------------------------------------------------------------------------

//        void delete_pointers();

//-------------------------------------------------------------------------------------------------

        void free_matrix_memory()
        {
            if ( mJacobianExist )
            {
                mJacobian.resize( 0, 0 );

                mJacobianExist = false;
            }
            if ( mResidualExist )
            {
                mResidual.resize( 0, 0 );

                mResidualExist = false;
            }
        };

//-------------------------------------------------------------------------------------------------

        void initialize_set()
        {
        };

//------------------------------------------------------------------------------

        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Set::finalize(), not implemented");
        };

//------------------------------------------------------------------------------

        void get_dof_types( moris::Cell< enum Dof_Type > & aDofType )
        {
            aDofType = mEqnObjDofTypeList;
        }

//------------------------------------------------------------------------------

        void set_model_solver_interface( Model_Solver_Interface * aModelSolverInterface)
        {
            mModelSolverInterface = aModelSolverInterface;
        };

//------------------------------------------------------------------------------

        Model_Solver_Interface * get_model_solver_interface()
        {
            return mModelSolverInterface;
        };

//------------------------------------------------------------------------------

        uint get_num_equation_objects()
        {
            return mEquationObjList.size();
        };

//------------------------------------------------------------------------------

        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mEquationObjList;
        };

//------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mEqnObjDofTypeList;
        }
//------------------------------------------------------------------------------

    void create_dof_assembly_map_EQ();

//------------------------------------------------------------------------------

    moris::Matrix< DDSMat > & get_dof_assembly_map_2()
    {
        return mDofAssemblyMap_2;
    }

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_ */
