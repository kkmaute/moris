#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_

#include "cl_Matrix.hpp"
#include "../additional/cl_GEN_Matrix_Base_Utilities.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        enum class Phase_Table_Structure
        {
            EXP_BASE_2
        };

        //--------------------------------------------------------------------------------------------------------------

        class Phase_Table
        {

        private:
            moris::moris_index mNumPhases; // number of phases
            Cell<std::string> mPhaseNames; // phase names
            enum Phase_Table_Structure mPhaseTableStructure; // enum defining phase table structure
            moris::Matrix<moris::IndexMat> mPhaseTable; // the actual phase table

            /**
             * Private constructor to eliminate redundant code
             *
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(std::string aStructure, Cell<std::string> aPhaseNames);

        public:
            /**
             * Constructor with explicitly defined phase table
             *
             * @param aPhaseTable Explicit phase table
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(Matrix<IndexMat> aPhaseTable, std::string aStructure, Cell<std::string> aPhaseNames = {});

            /**
             * Constructor for phase table built on structure
             *
             * @param aNumPhi Number of fields
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(uint aNumPhi, std::string aStructure, Cell<std::string> aPhaseNames = {});

            /**
             * Get the number of phases
             *
             * @return Number of phases
             */
            moris::moris_index get_num_phases();

            /**
             * Gets the sign of a given phase and geometry
             *
             * @param aPhaseIndex Phase index
             * @param aGeometryIndex Geometry index
             * @return Phase table value
             */
            moris::moris_index
            get_phase_sign_of_given_phase_and_geometry(moris::moris_index aPhaseIndex,
                                                       moris::moris_index aGeometryIndex);

            /**
             * Get phase index based on entity phase info
             *
             * @param aEntityPhaseInfo Phase info
             * @return Phase index
             */
            moris::moris_index get_phase_index(moris::Matrix< moris::IndexMat > const & aEntityPhaseInfo);

            /**
             * Gets the name of a requested phase
             *
             * @param aPhaseIndex The index of the requested phase
             * @return Phase name
             */
            std::string const & get_phase_name(moris::moris_index const & aPhaseIndex);

        private:
            /**
             * Checks the structure of the stored phase table
             *
             * @return If the phase table is correct based on the underlying structure
             */
            bool check_phase_table_structure();

            /**
             * Set all of the phases to have default names (p_i)
             */
            void set_default_phase_names();

        };
    }
}

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_ */
