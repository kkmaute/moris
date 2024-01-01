/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Output_Options.hpp
 *
 */

#ifndef UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_
#define UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_

#include "cl_Vector.hpp"

namespace xtk
{

// This class is used for unzipping the mesh in XTK model.
struct Output_Options
{
public:
// Add node sets from background mesh
bool mAddNodeSets;

// Add side sets
bool mAddSideSets;

// Add phase field
bool mAddPhaseField;

// Tell the mesh to locally index entities
bool mInternalUseFlag;

// Split the interface block from the other blocks (needed for exodus writing on background hex8)
bool mSeparateInterfaceBlock;

bool mHaveInterface = true;

// split the background side sets into interface and non interface parts
bool mSplitBackgroundSideSet = true;

// Specify that the mesh has phase information
bool mHasPhaseInfo;

// add parallel fields
bool mAddParallelFields;

// Appendix for sets indicating material phase
std::string mMaterialAppendix;

// Appendix for sets indicating interface
std::string mInterfaceAppendix;

// Specify the phase field name
std::string mPhaseFieldName;

// Sensitivity Options
bool mPackageDxDpSparsely;
bool mPackageDxDpDensely;
std::string mDxDpName;
std::string mDxDpIndicesName;
std::string mDxDpNumIndicesName;

// Other fields to add to the mesh
moris::Vector<std::string> mRealNodeExternalFieldNames;
moris::Vector<std::string> mIntNodeExternalFieldNames;

moris::Vector<std::string> mRealElementExternalFieldNames;
moris::Vector<std::string> mIntElementExternalFieldNames;

// Add cluster information to STK integration Mesh
bool mAddClusters;

Output_Options():
    mAddNodeSets(true),
    mAddSideSets(true),
    mAddPhaseField(false),
    mInternalUseFlag(false),
    mSeparateInterfaceBlock(true),
    mHasPhaseInfo(true),
    mAddParallelFields(false),
    mMaterialAppendix("_mat_"),
    mInterfaceAppendix("_i"),
    mPhaseFieldName("phase"),
    mPackageDxDpSparsely(true),
    mPackageDxDpDensely(false),
    mDxDpName("dxdp_"),
    mDxDpIndicesName("dxdp_inds_"),
    mDxDpNumIndicesName("dxdp_ninds"),
    mRealNodeExternalFieldNames({}),
    mIntNodeExternalFieldNames({}),
    mRealElementExternalFieldNames({}),
    mIntElementExternalFieldNames({}),
    mAddClusters(false),
    mOutputAllPhases(true)
{

}

// Ask whether I should output a given phase
bool output_phase(size_t const & aPhaseIndex) const
{
    if(mOutputAllPhases || mPhasesToOutput(aPhaseIndex) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}
/*
* Modify which phases are to be outputted
* aNumPhase -number of possible phae indices
* aPhasesToOutput - Vector of phase indices to output
*/
void change_phases_to_output(size_t const & aNumPhases,
                             moris::Vector<size_t> const & aPhasesToOutput)
{
    MORIS_ASSERT(mOutputAllPhases, "Phases have already been added, please only call this function once");
    mPhasesToOutput = moris::Vector<size_t>(aNumPhases,0);

    mNumPhasesToOutput = aPhasesToOutput.size();

    for(size_t i = 0; i<aPhasesToOutput.size(); i++)
    {
        mPhasesToOutput(aPhasesToOutput(i)) = 1;
    }

    // Don't output all phases
    mOutputAllPhases = false;
}

moris::uint
num_phases_to_output() const
{
    return mNumPhasesToOutput;
}

bool
output_all_phases() const
{
    return mOutputAllPhases;
}

private:

bool                mOutputAllPhases;
moris::uint         mNumPhasesToOutput;
moris::Vector<size_t> mPhasesToOutput;

};
}

#endif /* UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_ */

