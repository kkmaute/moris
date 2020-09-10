#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {

        // ----------------------------------------------------------------------------------

        enum class Phase_Table_Structure
        {
            EXP_BASE_2,
            INVALID
        };

        // ----------------------------------------------------------------------------------

        class Phase_Table
        {

        public:
            // ----------------------------------------------------------------------------------

            /**
             * Constructor with explicitly defined phase table
             *
             * @param aPhaseTable Explicit phase table
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table( Matrix<IndexMat>      aPhaseTable, 
                         Phase_Table_Structure aStructure, 
                         Cell<std::string>     aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            /**
             * Constructor for phase table built on structure
             *
             * @param aNumPhi Number of fields
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table( uint                       aNumPhi,
                         enum Phase_Table_Structure aStructure, 
                         Cell<std::string>          aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            /**
             * Constructor with explicitly defined phase table
             *
             * @param aPhaseTable Explicit phase table
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table( Matrix<IndexMat>  aPhaseTable, 
                         std::string       aStructure, 
                         Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            Phase_Table( Phase_Table_Structure aStructure, 
                         uint                  aNumPhases, 
                         Cell<std::string>     aPhaseNames);

            // ----------------------------------------------------------------------------------

            /**
             * Constructor for phase table built on structure
             *
             * @param aNumPhi Number of fields
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table( uint              aNumPhi, 
                         std::string       aStructure,
                         Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            /**
             * Get the number of phases
             *
             * @return Number of phases
             */
            moris_index 
            get_num_phases();

            // ----------------------------------------------------------------------------------

            /**
             * Gets the sign of a given phase and geometry
             *
             * @param aPhaseIndex Phase index
             * @param aGeometryIndex Geometry index
             * @return Phase table value
             */
            moris_index
            get_phase_sign_of_given_phase_and_geometry( moris_index aPhaseIndex,
                                                        moris_index aGeometryIndex);

            // ----------------------------------------------------------------------------------

            /**
             * Get phase index based on entity phase info
             *
             * @param aEntityPhaseInfo Phase info
             * @return Phase index
             */
            moris_index 
            get_phase_index( Matrix<IndexMat> const & aEntityPhaseInfo );

            // ----------------------------------------------------------------------------------

            /**
             * Gets the name of a requested phase
             *
             * @param aPhaseIndex The index of the requested phase
             * @return Phase name
             */
            std::string const & 
            get_phase_name( moris_index const & aPhaseIndex);

            // ----------------------------------------------------------------------------------

        private:
            Cell<std::string>     mPhaseNames;              // phase names
            Phase_Table_Structure mPhaseTableStructure; // enum defining phase table structure
            Matrix<IndexMat>      mPhaseTable; // the actual phase table
            moris_index           mNumPhases;              // number of phases

            /**
             * Checks the structure of the stored phase table
             *
             * @return If the phase table is correct based on the underlying structure
             */
            bool 
            check_phase_table_structure();

            // ----------------------------------------------------------------------------------

            /**
             * Set all of the phases to have default names (p_i)
             */
            void
            set_default_phase_names();

            // ----------------------------------------------------------------------------------

            /**
             * Gets the phase table structure enum from the corresponding string
             *
             * @param aStringStructure string defining the phase table structure
             * @return Phase_Table_Structure enum
             */
            enum Phase_Table_Structure 
            get_phase_table_structure( std::string aStringStructure );

            // ----------------------------------------------------------------------------------
        };
    } // namespace ge
} // namespace moris

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_ */
