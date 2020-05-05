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
#include "cl_FEM_Enums.hpp"
#include "cl_GEN_Pdv_Enums.hpp"


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
    class Equation_Model;
    enum class Dof_Type;

    //------------------------------------------------------------------------------
    /**
     * \brief Equation Set class
     */
    class Equation_Set
    {
    private:

    protected:
        moris::Cell< MSI::Equation_Object* > mEquationObjList;

        moris::Cell< Matrix< DDRMat > > mResidual;
        Matrix< DDRMat >                mJacobian;
        moris::Cell< Matrix< DDRMat > > mQI;
        moris::Cell< Matrix< DDRMat > > mdRdp;
        moris::Cell< moris::Cell< Matrix< DDRMat > > > mdQIdp;

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
        moris::Cell< moris::Cell< enum PDV > > mMasterDvTypes;
        moris::Cell< moris::Cell< enum PDV > > mSlaveDvTypes;

        // maps for the master and slave dv type
        moris::Matrix< DDSMat > mMasterDvTypeMap;
        moris::Matrix< DDSMat > mSlaveDvTypeMap;

        // FIXME map of master and slave dv types for assembly
        Cell< moris::Matrix< DDSMat > > mDvAssemblyMap;

        moris::Cell< moris::Cell< enum moris::fem::IQI_Type > > mRequestedIQITypes;
        moris::Cell< moris::Cell< moris_index > >               mRequestedIQITypeAssemblyMap;

        // initialization flag for jacobian, residual, QI, dRdp, dQIdp
        bool mJacobianExist = false;
        bool mResidualExist = false;
        bool mQIExist       = false;
        bool mdRdpExist     = false;
        bool mdQIdpExist    = false;

        Matrix< DDRMat > mTime;

        // unique list of dof and dv types
        moris::Cell< enum MSI::Dof_Type > mUniqueDofTypeList;
        moris::Cell< enum PDV >        mUniqueDvTypeList;

        // pointer to the model solver interface
        Model_Solver_Interface * mModelSolverInterface = nullptr;

        // FIXME pointer to the GEN MSI interface
        MSI::Design_Variable_Interface * mDesignVariableInterface = nullptr;

        bool mIsEmptySet = false;    //FIXME this flag is a hack. find better solution

        Matrix< DDRMat > * mSetElementalValues;
        Matrix< DDRMat > * mSetNodalValues;
        moris::real      * mSetGlobalValues;

        Matrix< DDUMat > mSetNodalCounter;

        uint tNumRHS = 1;

        MSI::Equation_Model * mEquationModel = nullptr;

        friend class MSI::Equation_Object;
        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Time_Sideset;
        friend class Element;

    public:

//------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        Equation_Set(){};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        virtual ~Equation_Set(){};

//------------------------------------------------------------------------------

        void set_equation_model(  MSI::Equation_Model * aEquationModel )
        {
            mEquationModel = aEquationModel;
        }

//------------------------------------------------------------------------------
        /**
         * get dof type list
         * @param[ in ] aIsMaster enum for master or slave
         */
        moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDofTypes;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDofTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_type_list - can only be MASTER or SLAVE");
                    return mMasterDofTypes;
                }
            }
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
                    MORIS_ERROR(false, "Equation_Set::get_dof_type_map - can only be MASTER or SLAVE");
                    return mMasterDofTypeMap;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get dof index for type
         * (return consecutive indices for master/slave,
         * i.e. index for slave starts at max index for master)
         * @param[ in ] aDofType  enum for dof type
         * @param[ in ] aIsMaster enum for master or slave
         * @param[ out ] sint     consecutive index for dof type
         */
        sint get_dof_index_for_type( enum MSI::Dof_Type     aDofType,
                                          mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if master dof type in list
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(),
                                  "Equation_Set::get_dof_index_for_type - master dof type does not exist in map." );

                    // get the index for master dof type from map
                    sint tMasterIndex = mMasterDofTypeMap( static_cast< int >( aDofType ) );

//                    // check if master dof type assigned in map
//                    MORIS_ASSERT( tMasterIndex != -1,
//                                  "Equation_Set::get_dof_index_for_type - master dof type not assigned in map." );

                    // return master dof type index
                    return tMasterIndex;
                }

                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if slave dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(),
                                  "Equation_Set::get_dof_index_for_type - slave dof type does not exist in map." );

                    // get the index for slave dof type from map
                    sint tSlaveIndex = mSlaveDofTypeMap( static_cast< int >( aDofType ) );

//                    // check if slave dof type assigned in map
//                    MORIS_ASSERT( tSlaveIndex != -1,
//                                  "Equation_Set::get_dof_index_for_type - slave dof type not assigned in map." );


                    if ( tSlaveIndex == -1 )
                    {
                        return tSlaveIndex;
                    }
                    else
                    {
                        // get maximum index from master map
                        moris::sint tMaxMasterIndex = mMasterDofTypeMap.max();

                        // check if master map was assigned
                        MORIS_ASSERT( tMaxMasterIndex != -1,
                                      "Equation_Set::get_dof_index_for_type - mMasterDofTypeMap is empty." );

                        // return slave dof type index
                        return tSlaveIndex + tMaxMasterIndex + 1;
                    }
                }

                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_index_for_type - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get dof index for type 1
         * (return non-consecutive indices for master/slave,
         *  i.e. index for slave restarts at 0)
         * @param[ in ]  aDofType  enum for dof type
         * @param[ in ]  aIsMaster enum for master or slave
         * @param[ out ] sint      non-consecutive index for dof type
         */
        sint get_dof_index_for_type_1( enum MSI::Dof_Type     aDofType,
                                            mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if master dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mMasterDofTypeMap.numel(),
                                  "Equation_Set::get_dof_index_for_type_1 - master dof type does not exist in map." );

//                    // check if master dof type assigned in map
//                    MORIS_ASSERT( mMasterDofTypeMap( static_cast< int >( aDofType ) ) != -1,
//                                  "Equation_Set::get_dof_index_for_type_1 - master dof type does not assigned in map." );

                    // return index for master dof type
                    return mMasterDofTypeMap( static_cast< int >( aDofType ) );
                }

                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if slave dof type in map
                    MORIS_ASSERT( static_cast< uint >( aDofType ) < mSlaveDofTypeMap.numel(),
                                  "Equation_Set::get_dof_index_for_type_1 - slave dof type does not exist in map." );

//                    // check if slave dof type assigned in map
//                    MORIS_ASSERT( mSlaveDofTypeMap( static_cast< int >( aDofType ) ) != -1,
//                                  "Equation_Set::get_dof_index_for_type_1 - slave dof type not assigned in map." );

                    // return index for slave dof type
                    return mSlaveDofTypeMap( static_cast< int >( aDofType ) );
                }

                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dof_index_for_type_1 - can only be MASTER or SLAVE." );
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get dv type list
         * @param[ in ] aIsMaster enum for master or slave
         */
        moris::Cell< moris::Cell< PDV > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    return mMasterDvTypes;
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    return mSlaveDvTypes;
                }
                default:
                {
                    MORIS_ERROR( false, "Equation_Set::get_dv_type_list - can only be MASTER or SLAVE.");
                    return mMasterDvTypes;
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
                    MORIS_ERROR(false, "Equation_Set::get_dv_type_map - can only be MASTER or SLAVE");
                    return mMasterDvTypeMap;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * get dv index for type
         * (return consecutive indices for master/slave,
         * i.e. index for slave starts at max index for master)
         * @param[ in ] aDvType  enum for dv type
         * @param[ in ] aIsMaster enum for master or slave
         * @param[ out ] sint     consecutive index for dv type
         */
        sint get_dv_index_for_type( enum PDV            aDvType,
                                         mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type - master dv type does not exist in map." );

//                    // check if dv type assigned in map
//                    MORIS_ASSERT( mMasterDvTypeMap( static_cast< int >( aDvType ) ),
//                                  "Equation_Set::get_dv_index_for_type - master dv type not assigned in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }

                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if dv type exists in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type(), slave dv type does not exist in map." );

//                    // check if dv type assigned in map
//                    MORIS_ASSERT( mSlaveDvTypeMap( static_cast< int >( aDvType ) ),
//                                  "Equation_Set::get_dv_index_for_type - slave dv type not assigned in map." );

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
         * get dv index for type 1
         * (return non-consecutive indices for master/slave,
         *  i.e. index for slave restarts at 0)
         * @param[ in ]  aDvType   enum for dv type
         * @param[ in ]  aIsMaster enum for master or slave
         * @param[ out ] sint      non-consecutive index for dv type
         */
        sint get_dv_index_for_type_1( enum PDV            aDvType,
                                           mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ):
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mMasterDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type_1 - master dv type does not exist in map." );

//                    // check if dv type is set in map
//                    MORIS_ASSERT( mMasterDvTypeMap( static_cast< int >( aDvType ) ),
//                                  "Equation_Set::get_dv_index_for_type_1 - master dv type not assigned in map." );

                    // return set index for dv type
                    return mMasterDvTypeMap( static_cast< int >( aDvType ) );
                }
                case( mtk::Master_Slave::SLAVE ):
                {
                    // check if dv type is set in map
                    MORIS_ASSERT( static_cast< uint >( aDvType ) < mSlaveDvTypeMap.numel(),
                                  "Equation_Set::get_dv_index_for_type_1 - slave dv type does not exist in map." );

//                    // check if dv type is set in map
//                    MORIS_ASSERT( mSlaveDvTypeMap( static_cast< int >( aDvType ) ),
//                                  "Equation_Set::get_dv_index_for_type_1 - slave dv type not assigned in map." );

                    // return set index for dv type
                    return mSlaveDvTypeMap( static_cast< int >( aDvType ) );
                }
                default:
                {
                    MORIS_ERROR(false, "Equation_Set::get_dv_index_for_type_1 - can only be MASTER or SLAVE.");
                    return 0;
                }
            }
        }

//-------------------------------------------------------------------------------------------------
        /**
         * free matrix memory,
         * i.e. free residual and jacobian matrices
         */
        void free_matrix_memory()
        {
            // if the jacobian matrix was created
            if ( mJacobianExist )
            {
                // resize it to 0x0
                mJacobian.resize( 0, 0 );

                // reset the exist flag
                mJacobianExist = false;
            }

            // if the residual cell of matrices was created
            if ( mResidualExist )
            {
                // resize each matrix to 0x0
                for( auto & tResidual : mResidual )
                {
                    tResidual.resize( 0, 0 );
                }
                mResidual.clear();

                // reset the exist flag
                mResidualExist = false;
            }

            // free additional memory
            this->free_memory();
        };

//-------------------------------------------------------------------------------------------------
        /**
         * initialize set
         */
        virtual void initialize_set( const bool aIsResidual,
                                     const bool aIsForward )
        {
            MORIS_ERROR( false, "Equation_Set::initialize_set - not implemented for virtual member function");
        };

//-------------------------------------------------------------------------------------------------
        /**
         * free memory
         */
        virtual void free_memory()
        {
            MORIS_ERROR( false, "Equation_Set::free_memory - not implemented for virtual member function");
        };

//-------------------------------------------------------------------------------------------------
        /**
         * finalize
         */
        virtual void finalize( MSI::Model_Solver_Interface * aModelSolverInterface )
        {
            MORIS_ERROR(false,"Equation_Set::finalize - not implemented for msi base class.");
        };

//-------------------------------------------------------------------------------------------------
        /**
         * set GEN/MSI interface
         * @param[ in ] aDesignVariableInterface a GEN/MSI interface pointer
         */
        virtual void set_dv_interface( MSI::Design_Variable_Interface * aDesignVariableInterface )
        {
            MORIS_ERROR( false, "Equation_Set::set_dv_interface - not implemented for msi base class." );
        };

//-------------------------------------------------------------------------------------------------
        /**
         * get residual
         */
        Cell< Matrix< DDRMat > > & get_residual()
        {
            return mResidual;
        };

//------------------------------------------------------------------------------
        /**
         * get residual dof assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_res_dof_assembly_map()
        {
            return mResDofAssemblyMap;
        }

//-------------------------------------------------------------------------------------------------
        /**
         * get QI
         */
        Cell< Matrix< DDRMat > > & get_QI()
        {
            return mQI;
        };

//-------------------------------------------------------------------------------------------------
        /**
         * get the element on the set
         * param[ out ] aElementType element type for the set
         */
        virtual enum fem::Element_Type get_element_type() const
        {
            MORIS_ERROR( false, "Equation_Set::get_element_type - not implemented for virtual member function");
            return fem::Element_Type::UNDEFINED;
        }

//------------------------------------------------------------------------------
        /**
         * get QI assembly map
         */
        moris_index get_QI_assembly_index( const enum Phase_Type    aPhaseType,
                                            const enum fem::IQI_Type aIQIType );

//-------------------------------------------------------------------------------------------------
        /**
         * get jacobian
         */
        Matrix< DDRMat > & get_jacobian()
        {
            return mJacobian;
        };

//------------------------------------------------------------------------------
        /**
         * get jacobian dof assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_jac_dof_assembly_map()
        {
            return mJacDofAssemblyMap;
        }

//-------------------------------------------------------------------------------------------------
        /**
         * get dRdp
         */
        moris::Cell< Matrix< DDRMat > > & get_dRdp()
        {
            return mdRdp;
        };

//-------------------------------------------------------------------------------------------------
        /**
         * get dQIdp
         */
        moris::Cell< moris::Cell< Matrix< DDRMat > > > & get_dQIdp()
        {
            return mdQIdp;
        };

//------------------------------------------------------------------------------
        /**
         * get dv assembly map
         */
        moris::Cell< moris::Matrix< DDSMat > > & get_dv_assembly_map()
        {
            return mDvAssemblyMap;
        }

//-------------------------------------------------------------------------------------------------
        /**
         * get number of right hand side
         */
        uint get_num_rhs()
        {
           return tNumRHS;
        };

//------------------------------------------------------------------------------
        /**
         * set model solver interface
         * @param[ in ] aModelSolverInterface a model solver interface pointer
         */
        void set_model_solver_interface( Model_Solver_Interface * aModelSolverInterface )
        {
            mModelSolverInterface = aModelSolverInterface;
        };

//------------------------------------------------------------------------------
        /**
         * get model solver interface
         * @param[ out ] aModelSolverInterface a model solver interface pointer
         */
        Model_Solver_Interface * get_model_solver_interface()
        {
            return mModelSolverInterface;
        };

//------------------------------------------------------------------------------
        /**
         * get number of equation objects
         */
        uint get_num_equation_objects()
        {
            return mEquationObjList.size();
        };

//------------------------------------------------------------------------------
        /**
         * get list of equation object pointers
         */
        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mEquationObjList;
        };

//------------------------------------------------------------------------------
        /**
         * get unique dof type list
         * @param[ out ] mUniqueDofTypeList a unique list of dof type
         */
        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mUniqueDofTypeList;
        }

//------------------------------------------------------------------------------
        /**
         * get number of unique dof types
         */
        moris::uint get_num_unique_dof_types()
        {
            return mUniqueDofTypeList.size();
        }

//------------------------------------------------------------------------------
        /**
         * get unique dv type list
         * @param[ out ] mUniqueDvTypeList a unique list of dv type
         */
        moris::Cell< enum PDV > & get_unique_dv_type_list()
        {
            return mUniqueDvTypeList;
        }

//------------------------------------------------------------------------------
        /**
         * get number of unique dv types
         */
        moris::uint get_num_unique_dv_types()
        {
            return mUniqueDvTypeList.size();
        }

//------------------------------------------------------------------------------
        /**
         * get requested dof types
         */
        moris::Cell < enum MSI::Dof_Type > get_requested_dof_types();

//------------------------------------------------------------------------------
        /**
         * get secundary dof types
         */
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > get_secundary_dof_types();

//------------------------------------------------------------------------------
        /**
         * set requested IQI types
         */
        void set_requested_IQI_types( const moris::Cell< moris::Cell< enum fem::IQI_Type > > & aRequestedIQITypes );

//------------------------------------------------------------------------------
        /**
          * get requested IQI types
          */
        const moris::Cell< moris::Cell< enum fem::IQI_Type > > & get_requested_IQI_types();

//------------------------------------------------------------------------------
        /**
          * create requested IQI type map
          */
        void create_requested_IQI_type_map();

//------------------------------------------------------------------------------
        /**
         * get requested dv types
         */
        moris::Cell < enum PDV > get_requested_dv_types();

//------------------------------------------------------------------------------
        /**
         * set visualization set
         * @param[ in ] aMeshIndex
         * @param[ in ] aVisMeshSet
         * @param[ in ] aOnlyPrimayCells
         */
        virtual void set_visualization_set( const uint              aMeshIndex,
                                                  moris::mtk::Set * aVisMeshSet,
                                            const bool              aOnlyPrimayCells )
        {
            MORIS_ASSERT( false, "Equation_Set::set_visualization_set(), not implemented for base class" );
        }

//------------------------------------------------------------------------------
        /**
         * compute quantity of interest
         * @param[ in ] aMeshIndex
         * @param[ in ] aElementFieldValues
         * @param[ in ] aNodalFieldValues
         * @param[ in ] aGlobalScalar
         * @param[ in ] aOutputType
         * @param[ in ] aFieldType
         */
        virtual void compute_quantity_of_interest( const uint              aMeshIndex,
                                                   Matrix< DDRMat >      * aElementFieldValues,
                                                   Matrix< DDRMat >      * aNodalFieldValues,
                                                   moris::real           * aGlobalScalar,
                                                   enum vis::Output_Type   aOutputType,
                                                   enum vis::Field_Type    aFieldType)
        {
            MORIS_ASSERT( false, "Equation_Set::compute_quantity_of_interest - not implemented for base class." );
        }

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_SET_HPP_ */
