/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Model.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_
#define PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CON/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

#include "cl_Library_IO.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
        class Dist_Map;
    }    // namespace sol
    //------------------------------------------------------------------------------
    namespace mtk
    {
        class Field;
    }
    //------------------------------------------------------------------------------

    namespace MSI
    {
        class Model_Solver_Interface;
        class MSI_Solver_Interface;
        class Equation_Set;
        class Equation_Object;
        class Design_Variable_Interface;
        enum class Dof_Type;
        //------------------------------------------------------------------------------

        class Equation_Model
        {
          protected:
            // list of equation sets
            moris::Cell< MSI::Equation_Set* > mFemSets;

            // list of equation objects
            moris::Cell< MSI::Equation_Object* > mFemClusters;

            // map from mesh set indices to fem set indices
            map< std::tuple< moris_index, bool, bool >, moris_index > mMeshSetToFemSetMap;

            // distributed solution vectors for current and previous time slabs
            sol::Dist_Vector* mSolutionVector                = nullptr;
            sol::Dist_Vector* mPrevSolutionVector            = nullptr;
            sol::Dist_Vector* mAdjointSolutionVector         = nullptr;
            sol::Dist_Vector* mPreviousAdjointSolutionVector = nullptr;

            sol::Dist_Map* mdQIdpMap = nullptr;

            moris::Cell< moris::Matrix< DDRMat > > mGlobalIQIVal;

            sol::Dist_Vector* mImplicitdQidp = nullptr;
            sol::Dist_Vector* mExplicitdQidp = nullptr;
            sol::Dist_Vector* mdQIdp         = nullptr;

            // matrices for current and previous time slabs
            Matrix< DDRMat > mTime;
            Matrix< DDRMat > mPrevTime;

            MSI::Design_Variable_Interface* mDesignVariableInterface = nullptr;

            bool mIsForwardAnalysis             = true;
            bool mIsOffDiagonalTimeContribution = false;

            moris::sint mNumSensitivityAnalysisRHS = -1;

            //------------------------------------------------------------------------------
            // Dummy Variables
            moris::Cell< std::string > mDummy;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Equation_Model(){};

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            virtual ~Equation_Model();

            //------------------------------------------------------------------------------

            virtual void free_memory() = 0;

            //------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            moris::sint get_num_rhs();

            //------------------------------------------------------------------------------
            /**
             * get equation sets for test
             */
            moris::Cell< MSI::Equation_Set* >&
            get_equation_sets()
            {
                return mFemSets;
            }

            //------------------------------------------------------------------------------
            /**
             * get equation objects
             */
            moris::Cell< MSI::Equation_Object* >&
            get_equation_objects()
            {
                return mFemClusters;
            }

            //------------------------------------------------------------------------------
            /**
             * MTK set to fem set index map
             */
            map< std::tuple< moris_index, bool, bool >, moris_index >&
            get_mesh_set_to_fem_set_index_map()
            {
                return mMeshSetToFemSetMap;
            }

            //------------------------------------------------------------------------------
            /**
             * finalize the equation sets
             * @param[ in ] aModelSolverInterface pointer to a model solver interface
             */
            virtual void
            finalize_equation_sets( MSI::Model_Solver_Interface* aModelSolverInterface )
            {
                MORIS_ERROR( false, "Equation_Model::finalize_equation_sets - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * set solution vector
             * @param[ in ] aSolutionVector distributed solution vector
             */
            void
            set_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * get solution vector
             * @param[ out ] aSolutionVector distributed solution vector
             */
            sol::Dist_Vector*
            get_solution_vector()
            {
                return mSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * set previous solution vector
             * @param[ in ] aSolutionVector previous distributed solution vector
             */
            void
            set_previous_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mPrevSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */
            sol::Dist_Vector*
            get_previous_solution_vector()
            {
                return mPrevSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * set sensitivity solution vector
             * @param[ in ] aSensitivitySolutionVector distributed solution vector for sensitivity
             */
            void
            set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mAdjointSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * set previous adjoint solution vector
             * @param[ in ] aSensitivitySolutionVector distributed solution vector for sensitivity
             */
            void
            set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mPreviousAdjointSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * get adjoint solution vector
             * @param[ out ] aSolutionVector adjoint distributed solution vector
             */
            sol::Dist_Vector*
            get_adjoint_solution_vector()
            {
                return mAdjointSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * get previous adjoint solution vector
             * @param[ out ] aSolutionVector previous adjoint distributed solution vector
             */
            sol::Dist_Vector*
            get_previous_adjoint_solution_vector()
            {
                return mPreviousAdjointSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * returns the implicit dQidu
             * @param[ out ] mImplicitdQidu returns a pointer to dQidu
             */
            sol::Dist_Vector*
            get_implicit_dQidp()
            {
                return mImplicitdQidp;
            }

            //------------------------------------------------------------------------------
            /**
             * returns the explicit dQidu
             * @param[ out ] mExplicitdQidu returns a pointer to dQidu
             */
            sol::Dist_Vector*
            get_explicit_dQidp()
            {
                return mExplicitdQidp;
            }

            //------------------------------------------------------------------------------
            /**
             * returns the dQIdp
             * @param[ out ] mdQidp returns a pointer to dQIdp
             */
            sol::Dist_Vector* get_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * set time for current time slab
             * @param[ in ] aTime matrix for time in current time slab
             */
            void
            set_time( Matrix< DDRMat >& aTime )
            {
                mTime = aTime;
            }

            //------------------------------------------------------------------------------
            /**
             * get time for current time slab
             * @param[ out ] mTime matrix for time in current time slab
             */
            Matrix< DDRMat >&
            get_time()
            {
                return mTime;
            }

            //------------------------------------------------------------------------------
            /**
             * set time for previous time slab
             * @param[ in ] aPrevTime matrix for time in previous time slab
             */
            void
            set_previous_time( Matrix< DDRMat >& aPrevTime )
            {
                mPrevTime = aPrevTime;
            }

            //------------------------------------------------------------------------------

            /**
             * get time for previous time slab
             * @param[ out ] mPrevTime matrix for time in previous time slab
             */
            Matrix< DDRMat >&
            get_previous_time()
            {
                return mPrevTime;
            }

            //------------------------------------------------------------------------------
            /**
             * set pointer to design variable interface
             * @param[ in ] aDesignVariableInterface pointer to design variable interface
             */
            void
            set_design_variable_interface( MSI::Design_Variable_Interface* aDesignVariableInterface )
            {
                mDesignVariableInterface = aDesignVariableInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * get pointer to design variable interface
             * @param[ out ] mDesignVariableInterface pointer to design variable interface
             */
            MSI::Design_Variable_Interface*
            get_design_variable_interface()
            {
                return mDesignVariableInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            virtual void
            set_requested_IQI_names( const moris::Cell< std::string >& aRequestedIQINames )
            {
                MORIS_ERROR( false, "Equation_Model::set_requested_IQI_names - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * get requested IQI names
             */
            virtual const moris::Cell< std::string >&
            get_requested_IQI_names()
            {
                MORIS_ERROR( false, "Equation_Model::get_requested_IQI_names - not implemented for base class." );
                return mDummy;
            }

            //------------------------------------------------------------------------------
            /**
             * set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            virtual void
            create_IQI_map()
            {
                MORIS_ERROR( false, "Equation_Model::create_requested_IQI_map - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * indicated that this equation model is used for the sensitivity analysis
             */
            void
            set_is_sensitivity_analysis()
            {
                this->reset();

                mIsForwardAnalysis = false;
            }

            //------------------------------------------------------------------------------
            /**
             * indicated that this equation model is used for the forward analysis
             */
            void
            set_is_forward_analysis()
            {
                this->reset();

                mIsForwardAnalysis = true;
            }

            //------------------------------------------------------------------------------
            /**
             * resets member variables of the equation object
             */
            virtual void reset(){};

            //------------------------------------------------------------------------------
            /**
             * report on assembly information
             */
            virtual void report_on_assembly(){};

            //------------------------------------------------------------------------------
            /**
             * returns if this is a forward analysis
             * @param[ out ] mIsForwardAnalysis bool true if forward analysis
             */
            bool
            get_is_forward_analysis()
            {
                return mIsForwardAnalysis;
            }

            //------------------------------------------------------------------------------
            /**
             * indicated that
             */
            void
            set_is_adjoint_off_diagonal_time_contribution(
                    const bool aIsOffDiagonalTimeContribution )
            {
                mIsOffDiagonalTimeContribution = aIsOffDiagonalTimeContribution;
            }

            //------------------------------------------------------------------------------
            /**
             * returns if this
             * @param[ out ] mIsForwardAnalysis
             */
            bool
            get_is_adjoint_off_diagonal_time_contribution()
            {
                return mIsOffDiagonalTimeContribution;
            }

            //------------------------------------------------------------------------------
            /**
             * initialize explicit and implicit dQidp
             */
            void initialize_explicit_and_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * compute implicit dQidp
             */
            void compute_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * compute explicit dQidp
             */
            void compute_explicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * compute explicit and implicit dQidp
             */
            void compute_explicit_and_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * initialize QI
             */
            virtual void
            initialize_IQIs()
            {
                MORIS_ERROR( 0, "Equation_Model::initialize_IQIs, not implemenetd in the base class" );
            }

            //------------------------------------------------------------------------------
            /**
             * compute QI
             */
            void compute_IQIs();

            //------------------------------------------------------------------------------
            /**
             * get QI global values
             * @param[ out ] mGlobalIQIVal cell filled with global QI values
             */
            moris::Cell< moris::Matrix< DDRMat > >&
            get_IQI_values()
            {
                return mGlobalIQIVal;
            }

            /**
             * Scale the IQIs according to user input. Default does nothing, scaling is done in child class.
             */
            virtual void normalize_IQIs();

            //------------------------------------------------------------------------------
            /**
             * get integration xyz active flags
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            get_integration_xyz_active_flags(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const moris::Cell< PDV_Type >& aPdvTypes,
                    Matrix< DDSMat >&              aIsActiveDv )
            {
                MORIS_ERROR( false,
                        "Equation_Model::get_integration_xyz_active_flags - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aXYZPdvIds  matrix to fill with pdv ids
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            get_integration_xyz_pdv_ids(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const moris::Cell< PDV_Type >& aPdvTypes,
                    Matrix< DDSMat >&              aXYZPdvIds )
            {
                MORIS_ERROR( false,
                        "Equation_Model::get_integration_xyz_pdv_ids - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             * @param[ in ] aXYZPdvIds   matrix to fill with pdv ids
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            get_integration_xyz_pdv_active_flags_and_ids(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const moris::Cell< PDV_Type >& aRequestedPdvTypes,
                    Matrix< DDSMat >&              aIsActiveDv,
                    Matrix< DDSMat >&              aXYZPdvIds )
            {
                MORIS_ERROR( false,
                        "Equation_Model::get_integration_xyz_pdv_active_flags_and_ids - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * reset integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices list of node indices to reset
             */
            virtual void
            reset_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >& aNodeIndices )
            {
                MORIS_ERROR( false,
                        "Equation_Model::reset_integration_xyz_pdv_assembly_indices - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * set integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            set_integration_xyz_pdv_assembly_index(
                    moris_index   aNodeIndex,
                    enum PDV_Type aPdvType,
                    moris_index   aXYZPdvAssemblyIndex )
            {
                MORIS_ERROR( false,
                        "Equation_Model::set_integration_xyz_pdv_assembly_index - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * get integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            get_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const moris::Cell< PDV_Type >& aRequestedPdvTypes,
                    Matrix< DDSMat >&              aXYZPdvAssemblyIndices )
            {
                MORIS_ERROR( false,
                        "Equation_Model::get_integration_xyz_pdv_assembly_indices - not implemented for base class." );
            }

            //------------------------------------------------------------------------------

            virtual void
            populate_fields()
            {
                MORIS_ERROR( false,
                        "Equation_Model::populate_fields - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * return fields
             */
            virtual moris::Cell< std::shared_ptr< mtk::Field > > get_fields() = 0;

            //------------------------------------------------------------------------------

            /**
             * initialize the FEM model from parameter lists + create the interpolation nodes & FEM sets
             * @param[ in ] aLibrary       a file path for property functions
             */
            virtual void
            initialize_from_inputfile( std::shared_ptr< Library_IO > aLibrary ) = 0;

            //------------------------------------------------------------------------------

            /**
             * @brief Set the dof type to Bspline mesh index map
             *
             * @param aDofTypeToBsplineMeshIndex
             */
            virtual void
            set_dof_type_to_Bspline_mesh_index( std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex ) = 0;
            //------------------------------------------------------------------------------

            /**
             * @brief set flag whether to use new ghost sets
             *
             * @param aUseNewGhostSets
             */
            virtual void
            set_use_new_ghost_sets( bool aUseNewGhostSets ) = 0;

            //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace MSI */
} /* namespace moris */

#endif /* PROJECTS_FEM_MDL_SRC_CL_MSI_EQUATION_MODEL_HPP_ */
