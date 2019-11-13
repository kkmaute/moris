/*
 * cl_XTK_Phase_Table.hpp
 *
 *  Created on: Oct 24, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_
#define PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_



// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes
//#include "cl_GEN_Enums.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_GEN_Matrix_Base_Utilities.hpp"
#include "cl_Matrix.hpp"

//Phase_Table_Structure

namespace moris
{
namespace ge
{

    class GEN_Phase_Table
    {
    public:
        GEN_Phase_Table(moris::Matrix< moris::IndexMat > const & aPhaseTable,
                    Cell<std::string> const & aPhaseNames,
                    enum Phase_Table_Structure const & aStructure = Phase_Table_Structure::EXP_BASE_2)
    {
            MORIS_ASSERT(aPhaseNames.size()==aPhaseTable.n_rows(),"Dimension mismatch between phase names and phase table");
            mPhaseTable = aPhaseTable.copy();
            mNumPhases = std::pow(2,mPhaseTable.n_cols());
            mPhaseTableStructure = aStructure;
            MORIS_ASSERT(this->check_phase_table_structure(),"Data structure does not adhere to the guidelines see wiki pdf on multi_phase for an explanation");
    }


        GEN_Phase_Table()
        {

        }

        GEN_Phase_Table(moris::moris_index aNumPhi,
                    enum Phase_Table_Structure const & aStructure = Phase_Table_Structure::EXP_BASE_2,
                    Cell<std::string> const aPhaseNames = {})
        {
            mPhaseTableStructure = aStructure;
            switch(aStructure)
            {
                case(Phase_Table_Structure::EXP_BASE_2):
                {
                    mNumPhases = std::pow(2,aNumPhi);


                    // Allocate phase table
                    mPhaseTable = moris::Matrix< moris::IndexMat >(mNumPhases,aNumPhi);
                    moris::moris_index tAlternator = mNumPhases;
                    moris::moris_index tCount = 0;
                    moris::moris_index tVal = 0;
                    for(moris::moris_index iC = 0; iC<aNumPhi; iC++)
                    {
                        tAlternator = tAlternator/2;
                        for(moris::moris_index iR = 0; iR<mNumPhases; iR ++)
                        {
                            mPhaseTable(iR,iC) = tVal;
                            if(tCount < tAlternator-1)
                            {

                                tCount ++;
                            }
                            else
                            {
                                if(tVal == 0)
                                {
                                    tVal = 1;
                                }
                                else
                                {
                                    tVal = 0;
                                }
                                tCount = 0;
                            }
                        }
                        tVal = 0;
                        tCount = 0;
                    }
                }
            }

            // Set default names
            if(aPhaseNames.size() == 0)
            {

                mPhaseNames = Cell<std::string>(mNumPhases,"    ");
                std::string tBase = "p_";

                for(moris::moris_index i = 0; i<mNumPhases; i++)
                {
                    mPhaseNames(i) = tBase + std::to_string(i);
                }

            }


        }

        moris::moris_index get_num_phases()
        {
            return mPhaseTable.n_rows();
        }

        moris::moris_index
        get_phase_sign_of_given_phase_and_geometry(moris::moris_index aPhaseIndex,
                                                   moris::moris_index aGeometryIndex)
        {
            return mPhaseTable(aPhaseIndex,aGeometryIndex);
        }

        moris::moris_index get_phase_index(moris::Matrix< moris::IndexMat > const & aEntityPhaseInfo)
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
                    break;
                }
            }
        }

        std::string const & get_phase_name(moris::moris_index const & aPhaseIndex)
        {
            MORIS_ASSERT(aPhaseIndex<(moris::moris_index)mPhaseNames.size(),"Phase index out of bounds");
            return mPhaseNames(aPhaseIndex);
        }


        // For test purposes.
        moris::Matrix< moris::IndexMat > const &
        get_phase_table_data()
        {
            return mPhaseTable;
        }
    private:
        moris::moris_index mNumPhases;
        Cell<std::string> mPhaseNames;
        enum Phase_Table_Structure mPhaseTableStructure;
        moris::Matrix< moris::IndexMat > mPhaseTable;

        bool
        check_phase_table_structure()
        {

            switch(mPhaseTableStructure)
            {
                case(Phase_Table_Structure::EXP_BASE_2):
                    {
                    moris::Matrix< moris::IndexMat > tRow(1,mPhaseTable.n_cols());
                    moris::moris_index tIndex = 0;
                    for(moris::size_t iR = 0; iR<mPhaseTable.n_rows(); iR++ )
                    {
                        tRow = mPhaseTable.get_row(iR);

                        tIndex = get_phase_index(tRow);

                        if(tIndex!=(moris::moris_index)iR)
                        {
                            return false;
                        }
                    }
                    return true;
                    }

                default:
                {
                    std::cout<<"Unhandled phase table structure";
                    return false;
                    break;
                }
            }
        }


    };
}
}

#endif /* PROJECTS_GEN_SRC_RIPPED_ADDITIONAL_CL_GEN_PHASE_TABLE_HPP_ */
