/*
 * cl_XTK_DxDp.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_SENSITIVITY_HPP_
#define SRC_XTK_CL_XTK_SENSITIVITY_HPP_

#include <unordered_map>

#include "xtk_typedefs.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
// XTKL: Containers
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"


namespace xtk
{

class Sensitivity
{
public:

    Sensitivity(){};

    Sensitivity(moris::moris_index aNumDesignVars,
                moris::moris_index aNumNodes,
                std::string aDxDpBaseName,
                std::string aDxDPADVIndsBaseName,
                std::string aDxDPNADVIndsBaseName,
                bool aOutputSparsely = false)
    {
        initialize(aNumDesignVars,
                   aNumNodes,
                   aDxDpBaseName,
                   aDxDPADVIndsBaseName,
                   aDxDPNADVIndsBaseName,
                   aOutputSparsely);
    }


    /*
     * Add the dxdp for a given node aNodeIndex. If the node has not been given sensitivity then add it to the map
     * If it has sensitivity but the compare flag is true then make sure these values are the same.
     */
    void
    add_node_sensitivity(moris::moris_index const & aNodeIndex,
                         moris::Matrix< moris::IndexMat > const & aDesignVarIndices,
                         moris::Matrix< moris::DDRMat > const & aDxDp,
                         bool aCompare = false)
    {
        auto tIterator = mDxDpIndices.begin();
        if(node_has_sensitivity(aNodeIndex,tIterator))
        {
            if(aCompare)
            {

            }
        }

        else
        {
            add_node_sensitivities_to_map(aNodeIndex,aDesignVarIndices,aDxDp);

        }
    }


    moris::Cell<moris::Matrix< moris::DDRMat >> const &
    get_sensitivity_data()
    {
        return  mDxDp;
    }

    moris::Matrix< moris::DDRMat > const &
    get_sensitivity_data(moris::moris_index aADVIndex)
    {
        return  mDxDp(aADVIndex);
    }

    moris::Matrix< moris::IndexMat > const &
    get_node_dxdp_map(moris::moris_index aNodeIndex)
    {
        auto tIterator = mDxDpIndices.begin();
        if(node_has_sensitivity(aNodeIndex,tIterator))
        {
            return (tIterator->second);
        }

        else
        {
            std::cout<<"Requesting node dxdp map on node that does not have sensitivity";
            return (tIterator->second);
        }
    }

    std::unordered_map<moris::moris_index, moris::Matrix< moris::IndexMat >> const &
    get_full_dxdp_map()
    {
        return mDxDpIndices;
    }


    /*
     * Returns the field name for the number of adv indices for a given node
     */
    std::string const & get_num_adv_ind_field_name()
    {
        return mDxDPNumIndName;
    }

    /*
     * Returns the field names for the adv indices
     */
    moris::Cell<std::string> const & get_adv_ind_field_name()
    {
        return mDxDpIndiceNames;
    }


    std::string const & get_field_name(moris::moris_index & aADVIndex)
    {
        return mDxDpNames(aADVIndex);
    }

    moris::Cell<std::string> const &
    get_all_field_names() const
    {
        return mDxDpNames;
    }
    /*
     * Resize out the extra size
     */
    void
    commit_sensitivities()
    {
         for(moris::moris_index i = 0; i< mNumDesignVars; i++)
         {
             mDxDp(i).resize(mDxDpCounter(0,i),3);
         }
    }


    /*
     * Returns the if the sensitivities should be outputted sparsely or densely.
     */

    bool
    output_sparesly()
    {
        return mOutputSparsely;
    }

    /*
     * Returns the maximum number of advs on a given nodes (only useful if sparse)
     */
    moris::moris_index
    get_max_num_advs()
    {
        return mMaxADVs;
    }

private:
    bool mOutputSparsely;
    moris::moris_index mNumDesignVars;

    // Maximum number of design variables a node depends on
    moris::moris_index mMaxADVs;

    // A vector of string names
    moris::Cell<std::string> mDxDpNames;

    // Sparse field names
    std::string mDxDPNumIndName;
    moris::Cell<std::string> mDxDpIndiceNames;

    // Contains the node indices corresponding to the rows in mDxDp (sparsity map)
    // Key - Node Index
    // Value - 2xnD where nD is the number of design variables the node depends on
    //         Row1 - Design Variable Index
    //         Row2 - Row Index in dXdP of the given design variable index
    bool mDenseIndices;
    std::unordered_map<moris::moris_index, moris::Matrix< moris::IndexMat >> mDxDpIndices;

    // Sensitivity Data
    moris::Cell<moris::Matrix< moris::DDRMat >> mDxDp;

    // Counts the locations used for resizing
    moris::Matrix< moris::IndexMat > mDxDpCounter;

    void
    initialize(moris::moris_index aNumDesignVars,
               moris::moris_index aNumNodes,
               std::string aDxDpBaseName,
               std::string aDxDPADVIndsBaseName,
               std::string aDxDPNADVIndsBaseName,
               bool aOutputSparsely)
    {
        mMaxADVs = 0;
        mNumDesignVars = aNumDesignVars;
        // Allocate space for names and declare names
        mOutputSparsely = aOutputSparsely;
        mDxDpNames = moris::Cell<std::string>(mNumDesignVars);
        for(moris::moris_index i = 0; i<mNumDesignVars;i++)
        {
            mDxDpNames(i) = aDxDpBaseName + std::to_string(i);
        }

        if(mOutputSparsely)
        {
            mDxDPNumIndName  = aDxDPNADVIndsBaseName;
            mDxDpIndiceNames = moris::Cell<std::string>(mNumDesignVars);
            for(moris::moris_index i = 0; i<mNumDesignVars;i++)
            {
                mDxDpIndiceNames(i) = aDxDPADVIndsBaseName + std::to_string(i);
            }
        }

        // Allocate for DxDp data
        mDxDpCounter = moris::Matrix< moris::IndexMat >(1,mNumDesignVars,0); // Needs to be filled with 0s
        mDxDp = moris::Cell<moris::Matrix< moris::DDRMat >>(mNumDesignVars);
        for(moris::moris_index i = 0; i<mNumDesignVars; i++)
        {
            mDxDp(i) = moris::Matrix< moris::DDRMat >(aNumNodes, 3); // 3 because its a vector (x,y,z)
        }
    }

    template<typename Iterator>
    bool
    node_has_sensitivity(moris::moris_index const & aNodeIndex,
                         Iterator & aIterator)
    {
        bool tHasSens = false;

        aIterator = mDxDpIndices.find(aNodeIndex);

        if(aIterator != mDxDpIndices.end())
        {
            tHasSens = true;
        }
        return tHasSens;
    }


    void
    add_node_sensitivities_to_map(moris::moris_index const & aNodeIndex,
                                  moris::Matrix< moris::IndexMat > const & aDesignVarIndices,
                                  moris::Matrix< moris::DDRMat > const & aDxDp)
    {
        // Allocate space
        moris::moris_index tNumDesignVars = aDesignVarIndices.n_cols();
        if(tNumDesignVars > mMaxADVs)
        {
            mMaxADVs = tNumDesignVars;
        }

        moris::Matrix< moris::IndexMat > tDxDpIndices(2,tNumDesignVars);
        tDxDpIndices.fill(0.0);
        conservative_copy(aDesignVarIndices, tDxDpIndices);

        // Add to map
        mDxDpIndices[aNodeIndex] = tDxDpIndices;

        // Add Data to DxDp and Increment Counts

        for(moris::moris_index i = 0; i<tNumDesignVars; i++)
        {
            moris::moris_index const & tDVIndex = aDesignVarIndices(0,i);

            if(!mOutputSparsely)
            {
                moris::moris_index & tCount = mDxDpCounter(0,tDVIndex);
                replace_row(i, aDxDp, tCount,mDxDp(tDVIndex));
                tDxDpIndices(1,i) = tCount;
                tCount++;
            }

            else
            {
                moris::moris_index & tCount = mDxDpCounter(0,i);
                replace_row(i, aDxDp, tCount,mDxDp(i));
                tDxDpIndices(1,i) = tCount;
                tCount++;
            }
        }
    }


};
}



#endif /* SRC_XTK_CL_XTK_SENSITIVITY_HPP_ */
