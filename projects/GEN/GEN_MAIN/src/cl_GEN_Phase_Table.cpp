#include "cl_GEN_Phase_Table.hpp"
#include "fn_linspace.hpp"


namespace moris
{
    namespace ge
    {
        Phase_Table::Phase_Table( moris::uint      const & aNumGeometries,
                                  moris::uint      const & aNumBulkPhases,
                                  Matrix<IndexMat> const & aGeomIndexToBulkPhase)
                                  :mNumGeometries(aNumGeometries)

        {   
            // fill with default 2n values
            if(aNumBulkPhases == MORIS_UINT_MAX)
            {
                mNumPhases = std::pow(2,mNumGeometries);
                mGeomValToBulkPhase = moris::linspace(0,mNumPhases-1,mNumPhases);
            }
            else
            {
              mNumPhases = aNumBulkPhases;
              mGeomValToBulkPhase= aGeomIndexToBulkPhase;
            }

        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(Phase_Table_Structure aStructure, uint aNumPhases, Cell<std::string> aPhaseNames)
        :mNumPhases(aNumPhases)
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
            MORIS_ASSERT(mPhaseNames.size() == mPhaseTable.n_rows(), "Dimension mismatch between phase names and phase table");
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
                uint                  aNumPhi,
                Phase_Table_Structure aStructure,
                Cell<std::string>     aPhaseNames)
        : Phase_Table(aStructure, std::pow(2, aNumPhi), aPhaseNames)
        {
 
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(
            Matrix<IndexMat> aPhaseTable, 
             std::string aStructure,
             Cell<std::string> aPhaseNames)
        : Phase_Table(aPhaseTable, get_phase_table_structure(aStructure), aPhaseNames)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Phase_Table::Phase_Table(uint aNumPhi, std::string aStructure, Cell<std::string> aPhaseNames)
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

        void
        Phase_Table::set_index_to_bulk_phase_map(Matrix<IndexMat> const & aIndexToBulkPhase)
        {
            MORIS_ASSERT(aIndexToBulkPhase.numel() == std::pow(2,aIndexToBulkPhase.numel()),"aIndexToBulkPhase needs to be of length 2^N_geom.");

            mGeomValToBulkPhase = aIndexToBulkPhase;
            mNumPhases = mGeomValToBulkPhase.max();
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index Phase_Table::get_num_phases()
        {
            return mNumPhases;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index
        Phase_Table::get_phase_sign_of_given_phase_and_geometry(moris::moris_index aPhaseIndex,
                moris::moris_index aGeometryIndex)
        {
            MORIS_ERROR(0,"REMOVED");
            return mPhaseTable(aPhaseIndex,aGeometryIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::moris_index 
        Phase_Table::get_phase_index(moris::Matrix<moris::IndexMat> const & aEntityPhaseInfo)
        {

            // MORIS_ASSERT(aEntityPhaseInfo.n_cols() == mPhaseTable.n_cols(), "Need information about every phase for this entity because using 2^n phase rule");
            moris::moris_index i = 0;
            for(moris::sint j = 0; j<mNumGeometries; j++)
            {
                MORIS_ASSERT(aEntityPhaseInfo(0,j) == 0 || aEntityPhaseInfo(0,j) == 1 ,"Phase not 1 or 0. Note: 1 corresponds to a positive and 0 to a negative");

                i += std::pow(2,mNumGeometries)/std::pow(2,j+1) *  aEntityPhaseInfo(0,j);
            }

            return mGeomValToBulkPhase(i);
        }

        //--------------------------------------------------------------------------------------------------------------

        std::string const & Phase_Table::get_phase_name(moris::moris_index const & aPhaseIndex)
        {
            MORIS_ASSERT(aPhaseIndex<(moris::moris_index)mPhaseNames.size(),"Phase index out of bounds");
            return mPhaseNames(aPhaseIndex);
        }

        //--------------------------------------------------------------------------------------------------------------
        void
        Phase_Table::print()
        {
            std::cout<<"Phase Table Info:"<<std::endl;
            std::cout<<"  Number of Geometries:  " <<mNumGeometries<<std::endl;
            std::cout<<"  Number of Bulk Phases: "<<mNumPhases<<std::endl;

            Matrix<IndexMat> tPhaseIndexTable = this->get_geometry_to_phase_index();

            // print the header 
            std::cout<<std::setw(8)<<"i" <<" | "<<std::setw(8)<<"Bp"<<" | ";
            for(moris::sint iG  = 0; iG < mNumGeometries; iG++)
            {
                std::cout<<std::setw(8)<<"G_" + std::to_string(iG) << " | ";
            }
            std::cout<<std::endl;

            for(moris::uint r  = 0; r < tPhaseIndexTable.n_rows(); r++)
            {            
                std::cout<<std::setw(8)<< r <<" | "<<std::setw(8)<<mGeomValToBulkPhase(r)<<" | ";

                for(moris::uint c  = 0; c < tPhaseIndexTable.n_cols(); c++)
                {
                    if(tPhaseIndexTable(r,c) > 0)
                    {
                        std::cout<<std::setw(8)<< "+" << " | ";
                    }
                    else
                    {
                        std::cout<<std::setw(8)<< "-" << " | ";
                    }
                }
                std::cout<<std::endl;
            }

            std::cout<<"  i   -> phase index  "<<std::endl;
            std::cout<<"  Bp  -> bulk phase index  "<<std::endl;
            std::cout<<"  G_j -> jth geometry  "      <<std::endl;
            std::cout<<"  +   -> greater than threshold  "<<std::endl;
            std::cout<<"  -   -> less than threshold  "<<std::endl;

        }

        //--------------------------------------------------------------------------------------------------------------
  
        Phase_Table_Structure 
        Phase_Table::get_phase_table_structure(std::string aStructure)
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
        
        //--------------------------------------------------------------------------------------------------------------
  
        Matrix<IndexMat>
        Phase_Table::get_geometry_to_phase_index()
        {
            moris_index tNumRows = std::pow(2,mNumGeometries);
            moris::Matrix< moris::IndexMat > tPhaseTable(tNumRows,mNumGeometries);
            moris::moris_index tAlternator = tNumRows;
            moris::moris_index tCount = 0;
            moris::moris_index tVal = 0;
            for(sint iC = 0; iC < mNumGeometries; iC++)
            {
                tAlternator = tAlternator/2;
                for(moris::moris_index iR = 0; iR<tNumRows; iR ++)
                {
                    tPhaseTable(iR,iC) = tVal;
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

            return tPhaseTable;
                
        }
        
        //--------------------------------------------------------------------------------------------------------------
  
    }
}
