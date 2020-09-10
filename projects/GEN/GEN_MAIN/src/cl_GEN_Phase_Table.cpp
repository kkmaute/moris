#include "cl_GEN_Phase_Table.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            Phase_Table_Structure aStructure, 
            uint                  aNumPhases, 
            Cell<std::string>     aPhaseNames)
        : mPhaseTableStructure(aStructure),
        mNumPhases(aNumPhases)
        {
            // Phase names
            if(aPhaseNames.size() == 0 && aNumPhases > 0)
            {
                this->set_default_phase_names();
            }
            else
            {
                mPhaseNames = aPhaseNames;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            Matrix<IndexMat>      aPhaseTable, 
            Phase_Table_Structure aStructure, 
            Cell<std::string>     aPhaseNames)
        : Phase_Table(aStructure, aPhaseTable.n_rows(), aPhaseNames)
        {
            mPhaseTable = aPhaseTable;
            MORIS_ASSERT(this->check_phase_table_structure(), "Data structure does not adhere to the guidelines see wiki pdf on multi_phase for an explanation");
            MORIS_ASSERT(mPhaseNames.size() == mPhaseTable.n_rows(), "Dimension mismatch between phase names and phase table");
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            uint                       aNumPhi, 
            enum Phase_Table_Structure aStructure, 
            Cell<std::string>          aPhaseNames)
        : Phase_Table(aStructure, std::pow(2, aNumPhi), aPhaseNames)
        {
            if (aNumPhi > 0)
            {
                // Determine structure for building
                switch(mPhaseTableStructure) {
                    case (Phase_Table_Structure::EXP_BASE_2):
                    {
                        // Allocate phase table
                        mPhaseTable = moris::Matrix<moris::IndexMat>(mNumPhases, aNumPhi);
                        moris::moris_index tAlternator = mNumPhases;
                        moris::moris_index tCount = 0;
                        moris::moris_index tVal = 0;
                        for (uint iC = 0; iC < aNumPhi; iC++) {
                            tAlternator = tAlternator / 2;
                            for (moris::moris_index iR = 0; iR < mNumPhases; iR++) {
                                mPhaseTable(iR, iC) = tVal;
                                if (tCount < tAlternator - 1) {
                                    tCount++;
                                } else {
                                    if (tVal == 0) {
                                        tVal = 1;
                                    } else {
                                        tVal = 0;
                                    }
                                    tCount = 0;
                                }
                            }
                            tVal = 0;
                            tCount = 0;
                        }
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR(false, "Phase table structure has not been implemented for generation using a number of fields.");
                    }
                }

                MORIS_ASSERT(this->check_phase_table_structure(), "Data structure does not adhere to the guidelines see wiki pdf on multi_phase for an explanation");
                MORIS_ASSERT(mPhaseNames.size() == mPhaseTable.n_rows(), "Dimension mismatch between phase names and phase table");
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            Matrix<IndexMat>  aPhaseTable, 
            std::string       aStructure,
            Cell<std::string> aPhaseNames)
        : Phase_Table(aPhaseTable, get_phase_table_structure(aStructure), aPhaseNames)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            uint              aNumPhi,
            std::string       aStructure, 
            Cell<std::string> aPhaseNames)
        : Phase_Table(aNumPhi, get_phase_table_structure(aStructure), aPhaseNames)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Phase_Table::set_default_phase_names()
        {
            mPhaseNames = Cell<std::string>(mNumPhases,"    ");
            std::string tBase = "p_";

            for(moris::moris_index i = 0; i < mNumPhases; i++)
            {
                mPhaseNames(i) = tBase + std::to_string(i);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index Phase_Table::get_num_phases()
        {
            return mPhaseTable.n_rows();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index
        Phase_Table::get_phase_sign_of_given_phase_and_geometry(
            moris::moris_index aPhaseIndex,
            moris::moris_index aGeometryIndex)
        {
            return mPhaseTable(aPhaseIndex,aGeometryIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index Phase_Table::get_phase_index(
            moris::Matrix< moris::IndexMat > const & aEntityPhaseInfo)
        {

            switch(mPhaseTableStructure)
            {
                case(Phase_Table_Structure::EXP_BASE_2):
                {
                    MORIS_ASSERT(aEntityPhaseInfo.n_cols() == mPhaseTable.n_cols(), "Need information about every phase for this entity because using 2^n phase rule");
                    moris::moris_index i = 0;
                    for(moris::size_t j = 0; j<mPhaseTable.n_cols(); j++)
                    {
                        MORIS_ASSERT(aEntityPhaseInfo(0,j) == 0 || aEntityPhaseInfo(0,j) == 1 ,"Phase not 1 or 0. Note: 1 corresponds to a positive and 0 to a negative");

                        i += mNumPhases/std::pow(2,j+1) *  aEntityPhaseInfo(0,j);
                    }
                    return i;
                }

                default:
                {
                    MORIS_ASSERT(0,"Unhandled phase table structure");
                    return 0;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        std::string const & Phase_Table::get_phase_name(
            moris::moris_index const & aPhaseIndex )
        {
            MORIS_ASSERT(aPhaseIndex<(moris::moris_index)mPhaseNames.size(),"Phase index out of bounds");
            return mPhaseNames(aPhaseIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Phase_Table::check_phase_table_structure()
        {
            switch(mPhaseTableStructure)
            {
                case (Phase_Table_Structure::EXP_BASE_2):
                    {
                    moris::Matrix< moris::IndexMat > tRow(1, mPhaseTable.n_cols());
                    moris::moris_index tIndex = 0;
                    for (moris::size_t iR = 0; iR < mPhaseTable.n_rows(); iR++ )
                    {
                        tRow = mPhaseTable.get_row(iR);

                        tIndex = get_phase_index(tRow);

                        if (tIndex!=(moris::moris_index)iR)
                        {
                            return false;
                        }
                    }
                    return true;
                    }

                default:
                {
                    return false;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table_Structure Phase_Table::get_phase_table_structure(
            std::string aStructure )
        {
            // Phase table structure
            if (aStructure == "exp_base_2")
            {
                return Phase_Table_Structure::EXP_BASE_2;
            }
            else
            {
                MORIS_ERROR(false, aStructure.append(" is not a valid phase table structure.").c_str());
                return Phase_Table_Structure::INVALID;
            }
        }

    }
}
