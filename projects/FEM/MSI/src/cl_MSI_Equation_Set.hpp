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

#include "cl_MTK_Enums.hpp"
#include "cl_MSI_Design_Variable_Interface.hpp"



namespace moris
{
    namespace mtk
    {
       class Cell;
       class Set;
    }
    namespace vis
    {
        enum class Output_Type;
        enum class Field_Type;
    }
    namespace MSI
    {
    class Model_Solver_Interface;
    class Equation_Object;
    class Design_Variable_Interface;
    enum class Dof_Type;
    enum class Dv_Type;
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

        // lists of master and slave groups of dof types
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mMasterDofTypes;
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > mSlaveDofTypes;

        // maps for the master and slave dof type
        moris::Matrix< DDSMat > mMasterDofTypeMap;
        moris::Matrix< DDSMat > mSlaveDofTypeMap;

        // map of master and slave dof types for assembly
        Cell< moris::Matrix< DDSMat > > mResDofAssemblyMap;
        Cell< moris::Matrix< DDSMat > > mJacDofAssemblyMap;

        // lists of master and slave groups of dv types
        moris::Cell< moris::Cell< enum GEN_DV > > mMasterDvTypes;
        moris::Cell< moris::Cell< enum GEN_DV > > mSlaveDvTypes;

        // maps for the master and slave dv type
        moris::Matrix< DDSMat > mMasterDvTypeMap;
        moris::Matrix< DDSMat > mSlaveDvTypeMap;

        // map of master and slave dv types for assembly
        Cell< moris::Matrix< DDSMat > > mDvAssemblyMap;

        Cell< moris::map< enum MSI::Dof_Type, moris::uint > > mRequestedTypeToIndexMap;

        bool mJacobianExist = false;
        bool mResidualExist = false;

        Matrix< DDRMat >mTime;

        // map of the element active dof types
        moris::Cell< enum MSI::Dof_Type > mEqnObjDofTypeList; // List of dof types of this equation obj
        moris::Cell< enum GEN_DV >  mUniqueDvTypeList; // List of dv types of this equation obj

        Model_Solver_Interface * mModelSolverInterface = nullptr;

        bool mIsEmptySet = false;    //FIXME this flag is a hack. find better solution

        Matrix< DDRMat > * mSetElementalValues;
        Matrix< DDRMat > * mSetNodalValues;
        moris::real      * mSetGlobalValues;

        Matrix< DDUMat > mSetNodalCounter;

        MSI::Design_Variable_Interface * mDesignVariableInterface = nullptr;

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
        /**
         * get set index for a dv type
         */
        sint get_dv_index_for_type( enum GEN_DV       aDvType,
                                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type(), dv type does not exist in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type(), dv type does not exist in map." );

                    // get the set index for dv type
                    sint tSlaveIndex = mSlaveDvTypeMap( static_cast< int >( aDvType ) );

                    // if index is -1
                    if ( tSlaveIndex == -1 )
                    {
                        return tSlaveIndex;
                    }
                    else
                    {
                        // get the max set index for dv types
                        moris::sint tMaxMasterIndex = mMasterDvTypeMap.max();

                        // check if mMasterDvTypeMap is set
                        MORIS_ASSERT( tMaxMasterIndex != -1,
                                      "Equation_Set::get_dv_index_for_type - mMasterDvTypeMap is empty." );

                        // return set index for dv type
                        return tSlaveIndex + tMaxMasterIndex + 1;
                    }
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dv_index_for_type - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get set index for dv type
         */
        sint get_dv_index_for_type_1( enum GEN_DV       aDvType,
                                      mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type - dv type does not exist in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type - dv type does not exist in map." );

                    // return set index for dv type
                    return mSlaveDvTypeMap( static_cast< int >( aDvType ) );
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dv_index_for_type - can only be MASTER or SLAVE.");
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
            this->free_memory();
        };

//-------------------------------------------------------------------------------------------------

        virtual void initialize_set( const bool aIsResidual )
        {
            MORIS_ERROR(false, "initialize_set(), not implemented for virtual memeber function");
        };
//-------------------------------------------------------------------------------------------------

        virtual void free_memory()
        {
            MORIS_ERROR(false, "free_memory(), not implemented for virtual memeber function");
        };

//-------------------------------------------------------------------------------------------------

        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Set::finalize(), not implemented");
        };

//-------------------------------------------------------------------------------------------------

         virtual void set_Dv_interface( MSI::Design_Variable_Interface * aDesignVariableInterface )
         {
           MORIS_ERROR( false, "");
         }

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

        moris::Cell< enum GEN_DV > & get_unique_dv_type_list()
        {
            return mUniqueDvTypeList;
        }

//------------------------------------------------------------------------------

        moris::uint get_num_unique_dof_types()
        {
            return mEqnObjDofTypeList.size();
        }

//------------------------------------------------------------------------------

        moris::uint get_num_unique_dv_types()
        {
            return mUniqueDvTypeList.size();
        }

//------------------------------------------------------------------------------

        virtual moris::Cell < enum MSI::Dof_Type > get_requested_dof_types()
        {
            MORIS_ERROR(false, "get_requested_dof_types(), not implemented for virtual member function");
            return moris::Cell< enum MSI::Dof_Type >( 0 );
        };

//------------------------------------------------------------------------------

        virtual moris::Cell< moris::Cell< enum MSI::Dof_Type > > get_secundary_dof_types()
        {
            MORIS_ERROR(false, "get_secundary_dof_types(), not implemented for virtual memeber function");
            return moris::Cell< moris::Cell< enum MSI::Dof_Type > >(0);
        };

//------------------------------------------------------------------------------

        virtual void set_visualization_set( const uint              aMeshIndex,
                                                  moris::mtk::Set * aVisMeshSet,
                                            const bool              aOnlyPrimayCells )
        {
            MORIS_ASSERT( false, "set_visualization_set(), not implemented for base class" );
        }

//------------------------------------------------------------------------------
        virtual void compute_quantity_of_interest( const uint            aMeshIndex,
                                                   Matrix< DDRMat >      * aElementFieldValues,
                                                   Matrix< DDRMat >      * aNodalFieldValues,
                                                   moris::real           * aGlobalScalar,
                                                   enum vis::Output_Type   aOutputType,
                                                   enum vis::Field_Type    aFieldType)
        {
            MORIS_ASSERT( false, "compute_quantity_of_interest(), not implemented for base class" );
        }

//------------------------------------------------------------------------------
        /**
         * get dof type map
         * @param[ in ] aIsMaster enum for master or slave
         */
        Matrix< DDSMat > & get_dof_type_map( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDofTypeMap;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDofTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return mMasterDofTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------
        /**
         * get dv type map
         * @param[ in ] aIsMaster enum for master or slave
         */
        Matrix< DDSMat > & get_dv_type_map( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDvTypeMap;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDvTypeMap;
                }
                default:
                {
                    MORIS_ERROR(false, "Set::get_dv_type_map - can only be MASTER or SLAVE");
                    return mMasterDvTypeMap;
                }
            }
        }

//------------------------------------------------------------------------------

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_ */
