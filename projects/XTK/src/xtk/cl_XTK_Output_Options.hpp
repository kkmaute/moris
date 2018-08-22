/*
 * cl__XTK_Output_Options.hpp
 *
 *  Created on: Oct 12, 2017
 *      Author: ktdoble
 */

#ifndef UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_
#define UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{

// This class is used for unzipping the mesh in XTK model.
template<typename Integer>
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

// Split the interface block from the other blocks (needed for exodus writing)
bool mSeparateInterface;

// Specify that the mesh has phase information
bool mHasPhaseInfo;

// Appendix for sets indicating material phase
std::string mMaterialAppendix;

// Appendix for sets indicating interface
std::string mInterfaceAppendix;

// Specify the phase field name
std::string mPhaseFieldName;

// Sensitivity Options
bool mPackageDxDpSparsely;
std::string mDxDpName;
std::string mDxDpIndicesName;
std::string mDxDpNumIndicesName;

// Other fields to add to the mesh
Cell<std::string> mRealNodeExternalFieldNames;
Cell<std::string> mIntNodeExternalFieldNames;

Cell<std::string> mRealElementExternalFieldNames;
Cell<std::string> mIntElementExternalFieldNames;


Output_Options():
    mAddNodeSets(true),
    mAddSideSets(true),
    mAddPhaseField(false),
    mInternalUseFlag(false),
    mSeparateInterface(false),
    mHasPhaseInfo(true),
    mMaterialAppendix("_mat_"),
    mInterfaceAppendix("_i"),
    mPhaseFieldName("phase"),
    mPackageDxDpSparsely(true),
    mDxDpName("dxdp_"),
    mDxDpIndicesName("dxdp_inds_"),
    mDxDpNumIndicesName("dxdp_ninds"),
    mRealNodeExternalFieldNames({}),
    mIntNodeExternalFieldNames({}),
    mRealElementExternalFieldNames({}),
    mIntElementExternalFieldNames({}),
    mOutputAllPhases(true)
{

}


// Ask whether I should output a given phase
bool output_phase(Integer const & aPhaseIndex) const
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
* aPhasesToOutput - Cell of phase indices to output
*/
void change_phases_to_output(Integer const & aNumPhases,
                             Cell<Integer> const & aPhasesToOutput)
{
    XTK_ASSERT(mOutputAllPhases, "Phases have already been added, please only call this function once");
    mPhasesToOutput = Cell<Integer>(aNumPhases,0);

    for(Integer i = 0; i<aPhasesToOutput.size(); i++)
    {
        mPhasesToOutput(aPhasesToOutput(i)) = 1;
    }

    // Don't output all phases
    mOutputAllPhases = false;
}

private:

bool mOutputAllPhases;
Cell<Integer> mPhasesToOutput;

};
}




#endif /* UNIT_TEST_SRC_XTK_CL_XTK_OUTPUT_OPTIONS_HPP_ */
