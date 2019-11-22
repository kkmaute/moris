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



        // maps for the master and slave dof type
        moris::Matrix< DDSMat > mMasterDofTypeMap;
        moris::Matrix< DDSMat > mSlaveDofTypeMap;

        // map of master and slave dof types for assembly
        Cell< moris::Matrix< DDSMat > > mDofAssemblyMap;

        Cell< moris::map< enum MSI::Dof_Type, moris::uint > > mRequestedTypeToIndexMap;

        bool mJacobianExist = false;
        bool mResidualExist = false;

        Matrix< DDRMat >mTime;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type > mEqnObjDofTypeList; // List of dof types of this equation obj

        Model_Solver_Interface * mModelSolverInterface = nullptr;

        bool mIsEmptySet = false;    //FIXME this flag is a hack. find better solution

        friend class MSI::Equation_Object;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

    public:

        sint get_dof_index_for_type( enum MSI::Dof_Type aDofType,
                                     mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
//                    MORIS_ASSERT( mMasterDofTypeMap( static_cast< int >( aDofType ) ) != -1, "Set::get_dof_index_for_type(), dof type does not exist in map " );

                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(), "Set::get_dof_index_for_type(), dof type does not exist in map " );

                    return mMasterDofTypeMap( static_cast< int >( aDofType ) );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
//                    MORIS_ASSERT( mSlaveDofTypeMap( static_cast< int >( aDofType ) ) != -1, "Set::get_dof_index_for_type(), dof type does not exist in map " );

                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(), "Set::get_dof_index_for_type(), dof type does not exist in map " );
                    sint tSlaveIndex = mSlaveDofTypeMap( static_cast< int >( aDofType ) );

                    if ( tSlaveIndex == -1 )
                    {
                        return tSlaveIndex;
                    }
                    else
                    {
                        moris::sint tMaxMasterIndex = mMasterDofTypeMap.max();

                        MORIS_ASSERT( tMaxMasterIndex != -1, "Set::get_dof_index_for_type(), mMasterDofTypeMap is seems to be empty " );

                        return tSlaveIndex + tMaxMasterIndex + 1;
                    }
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return 0;
                }
            }
        }

        sint get_dof_index_for_type_1( enum MSI::Dof_Type aDofType,
                                       mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
//                    MORIS_ASSERT( mMasterDofTypeMap( static_cast< int >( aDofType ) ) != -1, "Set::get_dof_index_for_type(), dof type does not exist in map " );
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(), "Set::get_dof_index_for_type(), dof type does not exist in map " );

                    return mMasterDofTypeMap( static_cast< int >( aDofType ) );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
//                    MORIS_ASSERT( mSlaveDofTypeMap( static_cast< int >( aDofType ) ) != -1, "Set::get_dof_index_for_type(), dof type does not exist in map " );
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(), "Set::get_dof_index_for_type(), dof type does not exist in map " );

                    return mSlaveDofTypeMap( static_cast< int >( aDofType ) );
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return 0;
                }
            }
        }

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

        virtual void initialize_set()
        {
            MORIS_ERROR(false, "initialize_set(), not implemented for virtual memeber function");
        };

//-------------------------------------------------------------------------------------------------

        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Set::finalize(), not implemented");
        };

//-------------------------------------------------------------------------------------------------

        Matrix< DDRMat > & get_residual()
        {
            return mResidual;
        };

//-------------------------------------------------------------------------------------------------

        Matrix< DDRMat > & get_jacobian()
        {
            return mJacobian;
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

        moris::uint get_num_unique_dof_types()
        {
            return mEqnObjDofTypeList.size();
        }


//------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_ */
