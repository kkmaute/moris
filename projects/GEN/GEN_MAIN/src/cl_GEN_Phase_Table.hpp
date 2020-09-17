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

            /**
             * @brief Allocate the phase table, if no num bulk-phases and map provided
             * assume that it is a standard 2^n table
             * @param aNumGeometries Number of geometries
             * @param aNumBulkPhase Number of bulk phases 
             * @param aGeomIndexToBulkPhase Geometric index to bulk phase map
             * 
             */
            Phase_Table( moris::uint      const & aNumGeometries,
                         moris::uint      const & aNumBulkPhase = MORIS_UINT_MAX,
                         Matrix<IndexMat> const & aGeomIndexToBulkPhase = {{}});

            // ----------------------------------------------------------------------------------

            /**
             * Constructor with explicitly defined phase table
             *
             * @param aPhaseTable Explicit phase table
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table( Matrix<IndexMat> aPhaseTable,
                         Phase_Table_Structure aStructure,
                         Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            /**
             * Constructor for phase table built on structure
             *
             * @param aNumPhi Number of fields
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(uint aNumPhi,
                        enum Phase_Table_Structure aStructure,
                        Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            /**
             * Constructor with explicitly defined phase table
             *
             * @param aPhaseTable Explicit phase table
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(Matrix<IndexMat> aPhaseTable,
                        std::string aStructure,
                        Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------

            Phase_Table(Phase_Table_Structure aStructure,
                        uint aNumPhases,
                        Cell<std::string> aPhaseNames);

            // ----------------------------------------------------------------------------------

            /**
             * Constructor for phase table built on structure
             *
             * @param aNumPhi Number of fields
             * @param aStructure Phase table structure
             * @param aPhaseNames (optional) Phase names
             */
            Phase_Table(uint aNumPhi,
                        std::string aStructure,
                        Cell<std::string> aPhaseNames = {});

            // ----------------------------------------------------------------------------------
            /**
             * @brief set the index to bulk phase map
             */ 
            void
            set_index_to_bulk_phase_map(Matrix<IndexMat> const & aIndexToBulkPhase);
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
            get_phase_sign_of_given_phase_and_geometry(moris_index aPhaseIndex,
                                                       moris_index aGeometryIndex);

            // ----------------------------------------------------------------------------------

            /**
             * Get phase index based on entity phase info
             *
             * @param aEntityPhaseInfo Phase info
             * @return Phase index
             */
            moris_index
            get_phase_index(Matrix<IndexMat> const &aEntityPhaseInfo);

            // ----------------------------------------------------------------------------------

            /**
             * Gets the name of a requested phase
             *
             * @param aPhaseIndex The index of the requested phase
             * @return Phase name
             */
            std::string const &
            get_phase_name(moris_index const &aPhaseIndex);

            // ----------------------------------------------------------------------------------

            /*!
            * @brief Print information for setting up phase table
            */
            void
            print();

        private:
            //fixme:remove
            Cell<std::string> mPhaseNames; // phase names
            Matrix<IndexMat> mPhaseTable;  // the phase table
                                           //fixme: end

            moris_index mNumPhases;               // number of bulk phases
            moris_index mNumGeometries;       // total number of geometries
            Matrix<IndexMat> mGeomValToBulkPhase; // Geometric value to Bulk phase

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
            get_phase_table_structure(std::string aStringStructure);

            // ----------------------------------------------------------------------------------

            Matrix<IndexMat>
            get_geometry_to_phase_index();
        };
    } // namespace ge
} // namespace moris

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_ */
