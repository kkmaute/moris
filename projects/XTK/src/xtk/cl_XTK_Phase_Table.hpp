/*
 * cl_XTK_Phase_Table.hpp
 *
 *  Created on: Oct 24, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_PHASE_TABLE_HPP_
#define SRC_XTK_CL_XTK_PHASE_TABLE_HPP_

#include "assert/fn_xtk_assert.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix.hpp"

//Phase_Table_Structure
#include "xtk/cl_XTK_Enums.hpp"


namespace xtk
{

    template<typename Integer, typename Integer_Matrix>
    class Phase_Table
    {
    public:
        Phase_Table(moris::Mat_New<Integer, Integer_Matrix> const & aPhaseTable,
                    Cell<std::string> const & aPhaseNames,
                    enum Phase_Table_Structure const & aStructure = Phase_Table_Structure::EXP_BASE_2)
    {
            XTK_ASSERT(aPhaseNames.size()==aPhaseTable.n_rows(),"Dimension mismatch between phase names and phase table");
            mPhaseTable = aPhaseTable.copy();
            mNumPhases = std::pow(2,mPhaseTable.n_cols());
            mPhaseTableStructure = aStructure;
            XTK_ASSERT(this->check_phase_table_structure(),"Data structure does not adhere to the guidelines see wiki pdf on multi_phase for an explanation");
    }


        Phase_Table()
        {

        }

        Phase_Table(Integer aNumPhi,
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
                    mPhaseTable = moris::Mat_New<Integer, Integer_Matrix>(mNumPhases,aNumPhi);
                    Integer tAlternator = mNumPhases;
                    Integer tCount = 0;
                    Integer tVal = 0;
                    for(Integer iC = 0; iC<aNumPhi; iC++)
                    {
                        tAlternator = tAlternator/2;
                        for(Integer iR = 0; iR<mNumPhases; iR ++)
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

                for(Integer i = 0; i<mNumPhases; i++)
                {
                    mPhaseNames(i) = tBase + std::to_string(i);
                }

            }


        }

        Integer get_num_phases()
        {
            return mPhaseTable.n_rows();
        }

        Integer get_phase_index(moris::Mat_New<Integer, Integer_Matrix> const & aEntityPhaseInfo)
        {

            switch(mPhaseTableStructure)
            {
                case(Phase_Table_Structure::EXP_BASE_2):
                    {
                    XTK_ASSERT(aEntityPhaseInfo.n_cols() == mPhaseTable.n_cols(), "Need information about every phase for this entity because using 2^n phase rule");
                    Integer i = 0;
                    for(Integer j = 0; j<mPhaseTable.n_cols(); j++)
                    {
                        XTK_ASSERT(aEntityPhaseInfo(0,j) == 0 || aEntityPhaseInfo(0,j) == 1 ,"Phase not 1 or 0. Note: 1 corresponds to a positive and 0 to a negative");

                        i += mNumPhases/std::pow(2,j+1) *  aEntityPhaseInfo(0,j);
                    }

                    return i;

                    }

                default:
                {
                    XTK_ERROR<<"Unhandled phase table structure";
                    return 0;
                    break;
                }
            }
        }

        std::string const & get_phase_name(Integer const & aPhaseIndex)
        {
            XTK_ASSERT(aPhaseIndex<mPhaseNames.size(),"Phase index out of bounds");
            return mPhaseNames(aPhaseIndex);
        }


        // For test purposes.
        moris::Mat_New<Integer, Integer_Matrix> const &
        get_phase_table_data()
        {
            return mPhaseTable;
        }
    private:
        Integer mNumPhases;
        Cell<std::string> mPhaseNames;
        enum Phase_Table_Structure mPhaseTableStructure;
        moris::Mat_New<Integer, Integer_Matrix> mPhaseTable;

        bool
        check_phase_table_structure()
        {

            switch(mPhaseTableStructure)
            {
                case(Phase_Table_Structure::EXP_BASE_2):
                    {
                    moris::Mat_New<Integer, Integer_Matrix> tRow(1,mPhaseTable.n_cols());
                    Integer tIndex = 0;
                    for(Integer iR = 0; iR<mPhaseTable.n_rows(); iR++ )
                    {
                        tRow = mPhaseTable.get_row(iR);

                        tIndex = get_phase_index(tRow);

                        if(tIndex!=iR)
                        {
                            return false;
                        }
                    }
                    return true;
                    }

                default:
                {
                    XTK_ERROR<<"Unhandled phase table structure";
                    return false;
                    break;
                }
            }
        }


    };
}

#endif /* SRC_XTK_CL_XTK_PHASE_TABLE_HPP_ */
