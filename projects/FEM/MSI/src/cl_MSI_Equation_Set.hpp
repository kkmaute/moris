/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Set.hpp
 *
 */

#ifndef SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_
#define SRC_FEM_CL_MSI_EQUATION_BLOCK_HPP_

#include "assert.h"
#include "cl_Communication_Tools.hpp"    //FEM/INT/src
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
    }    // namespace mtk
    namespace fem
    {
        class Field;
    }

    namespace vis
    {
        enum class Output_Type;
        enum class Field_Type;
    }    // namespace vis

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

            moris::Cell< Matrix< DDRMat > >                mResidual;
            Matrix< DDRMat >                               mJacobian;
            moris::Cell< Matrix< DDRMat > >                mQI;
            moris::Cell< Matrix< DDRMat > >                mdRdp;
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
            moris::Cell< moris::Cell< enum PDV_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< enum PDV_Type > > mSlaveDvTypes;

            // maps for the master and slave dv type
            moris::Matrix< DDSMat > mMasterDvTypeMap;
            moris::Matrix< DDSMat > mSlaveDvTypeMap;

            // lists of master and slave groups of field types
            moris::Cell< moris::Cell< enum mtk::Field_Type > > mMasterFieldTypes;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > mSlaveFieldTypes;

            // maps for the master and slave field type
            moris::Matrix< DDSMat > mMasterFieldTypeMap;
            moris::Matrix< DDSMat > mSlaveFieldTypeMap;

            // map of master and slave mat pdv types for assembly
            Cell< moris::Matrix< DDSMat > >                      mPdvMatAssemblyMap;
            moris::Matrix< DDSMat >                              mPdvMatAssemblyVector;
            std::map< std::pair< moris_index, PDV_Type >, uint > mPdvGeoAssemblyMap;
            moris::Matrix< DDSMat >                              mPdvGeoAssemblyVector;
            bool                                                 mPdvGeoAssemblyFlag = false;

            // Map from requested IQI Name to index.
            // I do not know if this is slow because the map is called per gauss point.
            // However as long as we go by name I do not see another way.
            moris::map< std::string, moris_index > mRequestedIQINamesAssemblyMap;

            // initialization flag for jacobian, residual, QI, dRdp, dQIdp
            bool mJacobianExist = false;
            bool mResidualExist = false;
            bool mQIExist       = false;
            bool mdRdpMatExist  = false;
            bool mdQIdpMatExist = false;

            // bool for time continuity
            bool mIsStaggered = false;

            Matrix< DDRMat > mTime;

            // unique list of dof and dv types
            moris::Cell< moris::Cell< enum MSI::Dof_Type > >   mUniqueDofTypeListMasterSlave;
            moris::Cell< moris::Cell< enum PDV_Type > >        mUniqueDvTypeListMasterSlave;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > mUniqueFieldTypeListMasterSlave;

            // unique list of dof and dv types. Master and Slave are combined
            moris::Cell< enum MSI::Dof_Type >   mUniqueDofTypeList;
            moris::Cell< enum PDV_Type >        mUniqueDvTypeList;
            moris::Cell< enum mtk::Field_Type > mUniqueFieldTypeList;

            // pointer to the model solver interface
            Model_Solver_Interface* mModelSolverInterface = nullptr;

            // FIXME pointer to the GEN MSI interface
            MSI::Design_Variable_Interface* mDesignVariableInterface = nullptr;

            bool mIsEmptySet = false;    // FIXME this flag is a hack. find better solution

            Matrix< DDRMat >* mSetElementalValues = nullptr;
            Matrix< DDRMat >* mSetNodalValues     = nullptr;
            Matrix< DDRMat >* mSetGlobalValues    = nullptr;

            uint tNumRHS = 1;

            MSI::Equation_Model* mEquationModel = nullptr;

            //! actual pdof values. Cells are for different multi-vectors
            moris::Cell< Matrix< DDRMat > > mPdofValues;

            //! previous pdof values
            moris::Cell< Matrix< DDRMat > > mPreviousPdofValues;

            //! adjoint pdof values
            moris::Cell< Matrix< DDRMat > > mAdjointPdofValues;

            //! adjoint pdof values
            moris::Cell< Matrix< DDRMat > > mPreviousAdjointPdofValues;

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

            void
            set_equation_model( MSI::Equation_Model* aEquationModel )
            {
                mEquationModel = aEquationModel;
            }

            //------------------------------------------------------------------------------

            MSI::Equation_Model*
            get_equation_model()
            {
                return mEquationModel;
            }

            //------------------------------------------------------------------------------
            /**
             * get dof type list
             * @param[ in ] aIsMaster enum for master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > >& get_dof_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dof type map
             * @param[ in ] aIsMaster enum for master or slave
             */
            Matrix< DDSMat >& get_dof_type_map(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dof index for type
             * (return consecutive indices for master/slave,
             * i.e. index for slave starts at max index for master)
             * @param[ in ] aDofType  enum for dof type
             * @param[ in ] aIsMaster enum for master or slave
             * @param[ out ] sint     consecutive index for dof type
             */
            sint get_dof_index_for_type(
                    enum MSI::Dof_Type aDofType,
                    mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dof index for type 1
             * (return non-consecutive indices for master/slave,
             *  i.e. index for slave restarts at 0)
             * @param[ in ]  aDofType  enum for dof type
             * @param[ in ]  aIsMaster enum for master or slave
             * @param[ out ] sint      non-consecutive index for dof type
             */
            sint get_dof_index_for_type_1(
                    enum MSI::Dof_Type aDofType,
                    mtk::Master_Slave  aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dv type list
             * @param[ in ] aIsMaster enum for master or slave
             */
            const moris::Cell< moris::Cell< PDV_Type > >& get_dv_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dv type map
             * @param[ in ] aIsMaster enum for master or slave
             */
            const Matrix< DDSMat >& get_dv_type_map(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dv index for type
             * (return consecutive indices for master/slave,
             * i.e. index for slave starts at max index for master)
             * @param[ in ] aDvType   enum for dv type
             * @param[ in ] aIsMaster enum for master or slave
             * @param[ out ] sint     consecutive index for dv type
             */
            sint get_dv_index_for_type(
                    enum PDV_Type     aDvType,
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get dv index for type 1
             * (return non-consecutive indices for master/slave,
             *  i.e. index for slave restarts at 0)
             * @param[ in ]  aDvType   enum for dv type
             * @param[ in ]  aIsMaster enum for master or slave
             * @param[ out ] sint      non-consecutive index for dv type
             */
            sint get_dv_index_for_type_1(
                    enum PDV_Type     aDvType,
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get field type list
             * @param[ in ] aIsMaster enum for master or slave
             */
            const moris::Cell< moris::Cell< mtk::Field_Type > >& get_field_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get fieldtype map
             * @param[ in ] aIsMaster enum for master or slave
             */
            const Matrix< DDSMat >& get_field_type_map(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get field index for type 1
             * (return non-consecutive indices for master/slave,
             *  i.e. index for slave restarts at 0)
             * @param[ in ]  aFieldType   enum for field type
             * @param[ in ]  aIsMaster enum for master or slave
             * @param[ out ] sint      non-consecutive index for field type
             */
            sint get_field_index_for_type_1(
                    enum mtk::Field_Type aFieldType,
                    mtk::Master_Slave    aIsMaster = mtk::Master_Slave::MASTER );

            //-------------------------------------------------------------------------------------------------
            /**
             * free matrix memory,
             * i.e. free residual and jacobian matrices
             */
            void free_matrix_memory();

            //-------------------------------------------------------------------------------------------------
            /**
             * initialize set
             */
            virtual void
            initialize_set(
                    const bool aIsStaggered,
                    const bool aTimeContinuityOnlyFlag )
            {
                MORIS_ERROR( false, "Equation_Set::initialize_set - not implemented for virtual member function" );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * free memory
             */
            virtual void
            free_memory()
            {
                MORIS_ERROR( false, "Equation_Set::free_memory - not implemented for virtual member function" );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * finalize
             */
            virtual void
            finalize( MSI::Model_Solver_Interface* aModelSolverInterface )
            {
                MORIS_ERROR( false, "Equation_Set::finalize - not implemented for msi base class." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * set GEN/MSI interface
             * @param[ in ] aDesignVariableInterface a GEN/MSI interface pointer
             */
            virtual void
            set_dv_interface( MSI::Design_Variable_Interface* aDesignVariableInterface )
            {
                MORIS_ERROR( false, "Equation_Set::set_dv_interface - not implemented for msi base class." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get residual
             */
            Cell< Matrix< DDRMat > >&
            get_residual()
            {
                return mResidual;
            }

            //------------------------------------------------------------------------------
            /**
             * get residual dof assembly map
             */
            moris::Cell< moris::Matrix< DDSMat > >&
            get_res_dof_assembly_map()
            {
                return mResDofAssemblyMap;
            };

            //-------------------------------------------------------------------------------------------------
            /**
             * get QI
             */
            Cell< Matrix< DDRMat > >&
            get_QI()
            {
                return mQI;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get the element on the set
             * param[ out ] aElementType element type for the set
             */
            virtual enum fem::Element_Type
            get_element_type() const
            {
                MORIS_ERROR( false, "Equation_Set::get_element_type - not implemented for virtual member function" );
                return fem::Element_Type::UNDEFINED;
            }

            //------------------------------------------------------------------------------
            /**
             * get QI assembly map
             */
            moris_index get_QI_assembly_index( const std::string& aIQIName );

            //-------------------------------------------------------------------------------------------------
            /**
             * get jacobian
             */
            Matrix< DDRMat >&
            get_jacobian()
            {
                return mJacobian;
            }

            //------------------------------------------------------------------------------
            /**
             * get jacobian dof assembly map
             */
            moris::Cell< moris::Matrix< DDSMat > >&
            get_jac_dof_assembly_map()
            {
                return mJacDofAssemblyMap;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get dRdp
             */
            moris::Cell< Matrix< DDRMat > >&
            get_drdp()
            {
                return mdRdp;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get dRdpMat
             */
            Matrix< DDRMat >&
            get_drdpmat()
            {
                return mdRdp( 0 );
            }

            //------------------------------------------------------------------------------
            /**
             * get dRdpMat pdv assembly map
             */
            moris::Cell< moris::Matrix< DDSMat > >&
            get_mat_pdv_assembly_map()
            {
                return mPdvMatAssemblyMap;
            }

            //------------------------------------------------------------------------------
            /**
             * get dRdpMat pdv assembly vector
             */
            Matrix< DDSMat >&
            get_mat_pdv_assembly_vector()
            {
                return mPdvMatAssemblyVector;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get dRdpGeo
             */
            Matrix< DDRMat >&
            get_drdpgeo()
            {
                return mdRdp( 1 );
            }

            //------------------------------------------------------------------------------
            /**
             * get dRdpGeo pdv assembly map
             */
            std::map< std::pair< moris_index, PDV_Type >, uint >&
            get_geo_pdv_assembly_map()
            {
                return mPdvGeoAssemblyMap;
            }

            //------------------------------------------------------------------------------
            /**
             * get dRdpGeo pdv assembly flag
             */
            bool
            get_geo_pdv_assembly_flag()
            {
                return mPdvGeoAssemblyFlag;
            }

            //------------------------------------------------------------------------------
            /**
             * get dRdpGeo pdv assembly vector
             */
            Matrix< DDSMat >&
            get_geo_pdv_assembly_vector()
            {
                return mPdvGeoAssemblyVector;
            }

            //-----------------------------------------------------------------------------------------
            /**
             * get dQIdp
             */
            moris::Cell< moris::Cell< Matrix< DDRMat > > >&
            get_dqidp()
            {
                return mdQIdp;
            }

            //-----------------------------------------------------------------------------------------
            /**
             * get dQIdp for material pdv
             */
            moris::Cell< Matrix< DDRMat > >&
            get_dqidpmat()
            {
                return mdQIdp( 0 );
            }

            //-----------------------------------------------------------------------------------------
            /**
             * get dQIdp for geometry pdv
             */
            moris::Cell< Matrix< DDRMat > >&
            get_dqidpgeo()
            {
                return mdQIdp( 1 );
            }

            //-----------------------------------------------------------------------------------------
            /**
             * get number of requested IQI for SA on set
             * @param[ out ] uint number of requested IQI for SA on set
             */
            virtual uint
            get_number_of_requested_IQIs()
            {
                MORIS_ERROR( false, "not implemented for base class." );
                return 0;
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * get number of right hand side
             */
            uint
            get_num_rhs()
            {
                return tNumRHS;
            }

            //------------------------------------------------------------------------------
            /**
             * set model solver interface
             * @param[ in ] aModelSolverInterface a model solver interface pointer
             */
            void
            set_model_solver_interface( Model_Solver_Interface* aModelSolverInterface )
            {
                mModelSolverInterface = aModelSolverInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * get model solver interface
             * @param[ out ] aModelSolverInterface a model solver interface pointer
             */
            Model_Solver_Interface*
            get_model_solver_interface()
            {
                return mModelSolverInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of equation objects
             */
            uint
            get_num_equation_objects()
            {
                return mEquationObjList.size();
            }

            //------------------------------------------------------------------------------
            /**
             * get list of equation object pointers
             */
            Cell< MSI::Equation_Object* >&
            get_equation_object_list()
            {
                return mEquationObjList;
            };

            //------------------------------------------------------------------------------
            /**
             * get unique dof type list
             * @param[ out ] mUniqueDofTypeList a unique list of dof type
             */
            const moris::Cell< enum MSI::Dof_Type >&
            get_unique_dof_type_list()
            {
                return mUniqueDofTypeList;
            }

            //------------------------------------------------------------------------------
            /**
             * get unique dof type list
             * @param[ out ] mUniqueDofTypeList a unique list of dof type
             */
            const moris::Cell< moris::Cell< enum MSI::Dof_Type > >&
            get_unique_master_slave_dof_type_list()
            {
                return mUniqueDofTypeListMasterSlave;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of unique dof types
             */
            moris::uint
            get_num_unique_dof_types()
            {
                return mUniqueDofTypeList.size();
            }

            //------------------------------------------------------------------------------
            /**
             * get unique dv type list
             * @param[ out ] mUniqueDvTypeList a unique list of dv type
             */
            const moris::Cell< enum PDV_Type >&
            get_unique_dv_type_list()
            {
                return mUniqueDvTypeList;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of unique dv types
             */
            moris::uint
            get_num_unique_dv_types()
            {
                return mUniqueDvTypeList.size();
            }

            //------------------------------------------------------------------------------
            /**
             * get unique field type list
             * @param[ out ] mUniqueFieldTypeList a unique list of field type
             */
            const moris::Cell< enum mtk::Field_Type >&
            get_unique_field_type_list()
            {
                return mUniqueFieldTypeList;
            }

            //------------------------------------------------------------------------------
            /**
             * get number of unique field types
             */
            moris::uint
            get_num_unique_field_types()
            {
                return mUniqueFieldTypeList.size();
            }

            //------------------------------------------------------------------------------
            /**
             * get requested dof types
             */
            const moris::Cell< enum MSI::Dof_Type >& get_requested_dof_types();

            //------------------------------------------------------------------------------
            /**
             * get secondary dof types
             */
            const moris::Cell< enum MSI::Dof_Type >& get_secondary_dof_types();

            //------------------------------------------------------------------------------
            /**
             * create requested IQI type map
             */
            void create_requested_IQI_type_map();

            //------------------------------------------------------------------------------
            /**
             * get requested dv types
             */
            moris::Cell< enum PDV_Type > get_requested_dv_types();

            //------------------------------------------------------------------------------
            /**
             * set visualization set
             * @param[ in ] aMeshIndex
             * @param[ in ] aVisMeshSet
             * @param[ in ] aOnlyPrimayCells
             */
            virtual void
            set_visualization_set(
                    const uint       aMeshIndex,
                    moris::mtk::Set* aVisMeshSet,
                    const bool       aOnlyPrimayCells )
            {
                MORIS_ASSERT( false, "Equation_Set::set_visualization_set(), not implemented for base class" );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest global
             * @param[ in ] aMeshIndex   mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues matrix to be filled with QI global values
             * @param[ in ] aQINames     list of QI names to compute
             */
            virtual void
            compute_quantity_of_interest_global(
                    const uint                        aMeshIndex,
                    Matrix< DDRMat >*                 aFieldValues,
                    const moris::Cell< std::string >& aQINames )
            {
                MORIS_ASSERT( false, "Equation_Set::compute_quantity_of_interest_global - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest nodal
             * @param[ in ] aMeshIndex   mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues matrix to be filled with QI nodal values
             * @param[ in ] aQINames     list of QI names to compute
             */
            virtual void
            compute_quantity_of_interest_nodal(
                    const uint                        aMeshIndex,
                    Matrix< DDRMat >*                 aFieldValues,
                    const moris::Cell< std::string >& aQINames )
            {
                MORIS_ASSERT( false, "Equation_Set::compute_quantity_of_interest_nodal - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest elemental
             * @param[ in ] aMeshIndex   mesh index to defined IG mesh to use
             * @param[ in ] aFieldValues matrix to be filled with QI elemental values
             * @param[ in ] aQINames     list of QI names to compute
             */
            virtual void
            compute_quantity_of_interest_elemental(
                    const uint                        aMeshIndex,
                    Matrix< DDRMat >*                 aFieldValues,
                    const moris::Cell< std::string >& aQINames )
            {
                MORIS_ASSERT( false, "Equation_Set::compute_quantity_of_interest_elemental - not implemented for base class." );
            }

            //------------------------------------------------------------------------------

            virtual void
            populate_fields(
                    moris::Cell< std::shared_ptr< fem::Field > >& aFieldToPopulate,
                    moris::Cell< std::string > const &            aFieldIQINames )
            {
                MORIS_ERROR( false, "populate_fields(), no child implementation." );
            }

            //------------------------------------------------------------------------------

            virtual std::string
            get_set_name()
            {

                MORIS_ERROR( false, "get_set_name(), not implemented for base class." );
                return "";
            }

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    }    // namespace MSI
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_EQUATION_SET_HPP_ */
