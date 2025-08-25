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

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "GEN_Data_Types.hpp"

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
            Vector< MSI::Equation_Set* > mFemSets;

            // list of equation objects
            Vector< MSI::Equation_Object* > mFemClusters;

            // map from mesh set indices to fem set indices
            map< std::tuple< moris_index, bool, bool >, moris_index > mMeshSetToFemSetMap;

            // distributed solution vectors for current and previous time slabs
            sol::Dist_Vector* mSolutionVector                = nullptr;
            sol::Dist_Vector* mPrevSolutionVector            = nullptr;
            sol::Dist_Vector* mAdjointSolutionVector         = nullptr;
            sol::Dist_Vector* mPreviousAdjointSolutionVector = nullptr;
            sol::Dist_Vector* mEigenSolutionVector           = nullptr;

            sol::Dist_Map* mdQIdpMap = nullptr;

            Vector< moris::Matrix< DDRMat > > mGlobalIQIVal;

            sol::Dist_Vector* mImplicitdQidp = nullptr;
            sol::Dist_Vector* mExplicitdQidp = nullptr;
            sol::Dist_Vector* mdQIdp         = nullptr;

            // matrices for current and previous time slabs
            Matrix< DDRMat > mTime;
            Matrix< DDRMat > mPrevTime;

            std::shared_ptr< MSI::Design_Variable_Interface > mDesignVariableInterface = nullptr;

            bool mIsForwardAnalysis             = true;
            bool mIsAdjointSensitivityAnalysis  = true;
            bool mIsOffDiagonalTimeContribution = false;

            moris::sint mNumSensitivityAnalysisRHS = -1;

            //------------------------------------------------------------------------------

            std::shared_ptr< Vector< real > > mEigenValues = nullptr;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Equation_Model() {};

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            virtual ~Equation_Model();

            //------------------------------------------------------------------------------

            virtual void free_memory() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief get equation sets for test
             */
            moris::sint get_num_rhs();

            //------------------------------------------------------------------------------
            /**
             * @brief get equation sets for test
             */
            Vector< MSI::Equation_Set* >&
            get_equation_sets()
            {
                return mFemSets;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get equation objects
             */
            Vector< MSI::Equation_Object* >&
            get_equation_objects()
            {
                return mFemClusters;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief MTK set to fem set index map
             */
            map< std::tuple< moris_index, bool, bool >, moris_index >&
            get_mesh_set_to_fem_set_index_map()
            {
                return mMeshSetToFemSetMap;
            }

            //------------------------------------------------------------------------------

            virtual Matrix< DDSMat >
            get_XYZ_local_pdv_assembly_map() const = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief finalize the equation sets
             * @param[ in ] aModelSolverInterface pointer to a model solver interface
             */
            virtual void
            finalize_equation_sets( MSI::Model_Solver_Interface* aModelSolverInterface ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief set solution vector
             * @param[ in ] aSolutionVector distributed solution vector
             */
            void
            set_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get solution vector
             * @returns aSolutionVector distributed solution vector
             */
            sol::Dist_Vector*
            get_solution_vector()
            {
                return mSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set previous solution vector
             * @param[ in ] aSolutionVector previous distributed solution vector
             */
            void
            set_previous_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mPrevSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get previous solution vector
             * @returns aSolutionVector previous distributed solution vector
             */
            sol::Dist_Vector*
            get_previous_solution_vector()
            {
                return mPrevSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set eigen vectors
             * @param[ in ] aEigenVector distributed solution vector
             */
            void
            set_eigen_solution_vector( sol::Dist_Vector* aEigenSolutionVector )
            {
                mEigenSolutionVector = aEigenSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief Set the eigen values object
             *
             * @param aEigenValues
             */

            void set_eigen_values( std::shared_ptr< Vector< real > > aEigenValues )
            {
                mEigenValues = std::move( aEigenValues );
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get previous solution vector
             *
             * @returns aSolutionVector previous distributed solution vector
             */
            sol::Dist_Vector*
            get_eigen_solution_vector()
            {
                return mEigenSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * get previous solution vector
             * @param[ out ] aSolutionVector previous distributed solution vector
             */

            std::shared_ptr< Vector< real > >&
            get_eigen_values()
            {
                return mEigenValues;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set sensitivity solution vector
             * @param[ in ] aSensitivitySolutionVector distributed solution vector for sensitivity
             */
            void
            set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mAdjointSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set previous adjoint solution vector
             * @param[ in ] aSensitivitySolutionVector distributed solution vector for sensitivity
             */
            void
            set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
            {
                mPreviousAdjointSolutionVector = aSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get adjoint solution vector
             * @returns aSolutionVector adjoint distributed solution vector
             */
            sol::Dist_Vector*
            get_adjoint_solution_vector()
            {
                return mAdjointSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get previous adjoint solution vector
             * @returns aSolutionVector previous adjoint distributed solution vector
             */
            sol::Dist_Vector*
            get_previous_adjoint_solution_vector()
            {
                return mPreviousAdjointSolutionVector;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief returns the implicit dQidu
             * @returns  a pointer to dQidu
             */
            sol::Dist_Vector*
            get_implicit_dQidp()
            {
                return mImplicitdQidp;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief returns the explicit dQidu
             * @returns a pointer to dQidu
             */
            sol::Dist_Vector*
            get_explicit_dQidp()
            {
                return mExplicitdQidp;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief returns the dQIdp
             * @returns a pointer to dQIdp
             */
            sol::Dist_Vector* get_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * @brief set time for current time slab
             * @param[ in ] aTime matrix for time in current time slab
             */
            void
            set_time( Matrix< DDRMat >& aTime )
            {
                mTime = aTime;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get time for current time slab
             * @returns mTime matrix for time in current time slab
             */
            Matrix< DDRMat >&
            get_time()
            {
                return mTime;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set time for previous time slab
             * @param[ in ] aPrevTime matrix for time in previous time slab
             */
            void
            set_previous_time( Matrix< DDRMat >& aPrevTime )
            {
                mPrevTime = aPrevTime;
            }

            //------------------------------------------------------------------------------

            /**
             * @brief get time for previous time slab
             * @returns mPrevTime matrix for time in previous time slab
             */
            Matrix< DDRMat >&
            get_previous_time()
            {
                return mPrevTime;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set pointer to design variable interface
             * @param[ in ] aDesignVariableInterface pointer to design variable interface
             */
            void
            set_design_variable_interface( std::shared_ptr< MSI::Design_Variable_Interface > aDesignVariableInterface )
            {
                mDesignVariableInterface = aDesignVariableInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief get pointer to design variable interface
             * @returns pointer to design variable interface
             */
            std::shared_ptr< const MSI::Design_Variable_Interface >
            get_design_variable_interface()
            {
                return mDesignVariableInterface;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            virtual void set_requested_IQI_names( const Vector< std::string >& aRequestedIQINames ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief get requested IQI names
             */
            virtual const Vector< std::string >&
            get_requested_IQI_names() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief set requested IQI names
             * @param[ in ] aRequestedIQINames List of requested IQI names
             */
            virtual void
            create_IQI_map() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief indicated that this equation model is used for the sensitivity analysis
             */
            void
            set_sensitivity_analysis_type( bool tIsAdjointSensitivityAnalysis )
            {
                this->reset();

                mIsForwardAnalysis            = false;
                mIsAdjointSensitivityAnalysis = tIsAdjointSensitivityAnalysis;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief indicated that this equation model is used for the forward analysis
             */
            void
            set_is_forward_analysis()
            {
                this->reset();

                mIsForwardAnalysis = true;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief resets member variables of the equation object
             */
            virtual void reset() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief report on assembly information
             */
            virtual void report_on_assembly() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief returns if this is a forward analysis
             * @returns mIsForwardAnalysis bool true if forward analysis
             */
            bool
            is_forward_analysis() const
            {
                return mIsForwardAnalysis;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief returns if adjoint sensitivity analysis is used; if not it is direct
             * @returns mIsForwardAnalysis bool true if forward analysis
             */
            bool
            is_adjoint_sensitivity_analysis() const
            {
                return mIsAdjointSensitivityAnalysis;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief indicated that
             */
            void
            set_is_adjoint_off_diagonal_time_contribution(
                    const bool aIsOffDiagonalTimeContribution )
            {
                mIsOffDiagonalTimeContribution = aIsOffDiagonalTimeContribution;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief returns if this
             * @returns mIsForwardAnalysis
             */
            bool
            get_is_adjoint_off_diagonal_time_contribution() const
            {
                return mIsOffDiagonalTimeContribution;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief initialize explicit and implicit dQidp
             */
            void initialize_explicit_and_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * compute implicit dQidp
             */
            //            void compute_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * compute explicit dQidp
             */
            //            void compute_explicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * @brief compute explicit and implicit dQidp
             */
            void compute_explicit_and_implicit_dQIdp();

            //------------------------------------------------------------------------------
            /**
             * @brief initialize QI
             */
            virtual void initialize_IQIs() = 0;

            /**
             * @brief computes the "requested" IQIs which are set in GEN parameter list
             */
            void compute_IQIs();

            //------------------------------------------------------------------------------
            /**
             * @brief get QI global values
             *
             * @returns mGlobalIQIVal cell filled with global QI values
             */
            Vector< moris::Matrix< DDRMat > >&
            get_IQI_values()
            {
                return mGlobalIQIVal;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief Scale the IQIs according to user input. Default does nothing, scaling is done in child class.
             */
            virtual void normalize_IQIs() {};

            //------------------------------------------------------------------------------
            /**
             * @brief get integration xyz active flags
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void get_integration_xyz_active_flags(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const Vector< gen::PDV_Type >& aPdvTypes,
                    Matrix< DDSMat >&              aIsActiveDv ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aXYZPdvIds  matrix to fill with pdv ids
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void get_integration_xyz_pdv_ids(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const Vector< gen::PDV_Type >& aPdvTypes,
                    Matrix< DDSMat >&              aXYZPdvIds ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief get integration xyz pdv ids
             * @param[ in ] aNodeIndices list of node indices
             * @param[ in ] aPdvTypes    list of pdv types
             * @param[ in ] aIsActiveDv  matrix to fill with 0/1 when pdv is active
             *                           ( tNumNodeIndices x tNumPdvTypes )
             * @param[ in ] aXYZPdvIds   matrix to fill with pdv ids
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void get_integration_xyz_pdv_active_flags_and_ids(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const Vector< gen::PDV_Type >& aRequestedPdvTypes,
                    Matrix< DDSMat >&              aIsActiveDv,
                    Matrix< DDSMat >&              aXYZPdvIds ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief reset integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices list of node indices to reset
             */
            virtual void reset_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >& aNodeIndices ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief set integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            set_integration_xyz_pdv_assembly_index(
                    moris_index        aNodeIndex,
                    enum gen::PDV_Type aPdvType,
                    moris_index        aXYZPdvAssemblyIndex ) = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief get integration xyz pdv local cluster assembly indices
             * @param[ in ] aNodeIndices           list of node indices
             * @param[ in ] aPdvTypes              list of pdv types
             * @param[ in ] aXYZPdvAssemblyIndices matrix to fill with assembly indices for pdv
             *                           ( tNumNodeIndices x tNumPdvTypes )
             */
            virtual void
            get_integration_xyz_pdv_assembly_indices(
                    const Matrix< IndexMat >&      aNodeIndices,
                    const Vector< gen::PDV_Type >& aRequestedPdvTypes,
                    Matrix< DDSMat >&              aXYZPdvAssemblyIndices ) = 0;

            //------------------------------------------------------------------------------

            virtual void
            populate_fields() = 0;

            /**
             * Updates the stored fields.
             */
            virtual void update_fields() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief return all fields
             */
            virtual Vector< std::shared_ptr< mtk::Field > >
            get_fields() = 0;

            //------------------------------------------------------------------------------
            /**
             * @brief Method to update the fem sets that need to be reinitialized in every newton iteration
             */
            virtual void update_equation_sets() {};

            //------------------------------------------------------------------------------
            /**
             * @brief initialize the FEM model from parameter lists + create the interpolation nodes & FEM sets
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

#endif /* PROJECTS_FEM_MDL_SRC_CL_MSI_MODEL_HPP_ */
