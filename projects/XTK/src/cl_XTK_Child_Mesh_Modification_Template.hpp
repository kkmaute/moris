/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Child_Mesh_Modification_Template.hpp
 *
 */

#ifndef SRC_XTK_CL_XTK_CHILD_MESH_MODIFICATION_TEMPLATE_HPP_
#define SRC_XTK_CL_XTK_CHILD_MESH_MODIFICATION_TEMPLATE_HPP_

#include "cl_Matrix.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_Mesh_Enums.hpp"
namespace xtk
{
class Mesh_Modification_Template
{
public:

    Mesh_Modification_Template()
    {

    }

    // Constructor when the template does not have node inheritance information
    Mesh_Modification_Template( moris::moris_index                       aParentElemInd,
                                moris::size_t                            aElemToReplace,
                                moris::Matrix< moris::IndexMat > const & aNodeInds,
                                moris::Matrix< moris::IndexMat > const & aParentEdgeInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentEdgeRanks,
                                moris::Matrix< moris::IndexMat > const & aParentFaceInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentFaceRanks,
                                enum TemplateType                        aTemplateType,
                                moris::size_t                            aPermutationId = 0):
                                    mParentNodeInds(0,0),
                                    mParentNodeRanks(0,0),
                                    mSpatialDimension(0),
                                    mElementTopology(CellTopology::INVALID)
    {
        mElemIndToReplace = aElemToReplace;
        mParentElemInd    = aParentElemInd;
        mNodeInds         = aNodeInds.copy();
        mParentEdgeInds   = aParentEdgeInds.copy();
        mParentEdgeRanks  = aParentEdgeRanks.copy();
        mParentFaceInds   = aParentFaceInds.copy();
        mParentFaceRanks  = aParentFaceRanks.copy();
        template_catalog(aPermutationId,aTemplateType);
    }

    // Constructor for when the template has node inheritance information
    Mesh_Modification_Template( moris::moris_index                       aParentElemInd,
                                moris::size_t                            aElemToReplace,
                                moris::Matrix< moris::IndexMat > const & aNodeInds,
                                moris::Matrix< moris::IndexMat > const & aParentNodeInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentNodeRanks,
                                moris::Matrix< moris::IndexMat > const & aParentEdgeInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentEdgeRanks,
                                moris::Matrix< moris::IndexMat > const & aParentFaceInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentFaceRanks,
                                enum TemplateType                        aTemplateType,
                                moris::size_t                            aPermutationId = 0)
    {
        MORIS_ASSERT(aTemplateType == TemplateType::REGULAR_SUBDIVISION_HEX8 || aTemplateType == TemplateType::REGULAR_SUBDIVISION_QUAD4 || aTemplateType ==TemplateType::TET_4 || aTemplateType == TemplateType::TRI_3,
                     "Entered node inheritance constructor for template without node inheritance. Currently, only the regular subdivision hex8 and tet templates have node inheritance");

        mTemplateType	  = aTemplateType;
        mElemIndToReplace = aElemToReplace;
        mParentElemInd    = aParentElemInd;
        mNodeInds         = aNodeInds.copy();
        mParentNodeInds   = aParentNodeInds.copy();
        mParentNodeRanks  = aParentNodeRanks.copy();
        mParentEdgeInds   = aParentEdgeInds.copy();
        mParentEdgeRanks  = aParentEdgeRanks.copy();
        mParentFaceInds   = aParentFaceInds.copy();
        mParentFaceRanks  = aParentFaceRanks.copy();
        template_catalog(aPermutationId,aTemplateType);
    }

    // Constructor for when the 2D template has node inheritance information
    // Separate constructor could not be implemented because non-inheritance constructor has same parameters
    /*Mesh_Modification_Template( moris::moris_index                       aParentElemInd,
                                moris::size_t                            aElemToReplace,
                                moris::Matrix< moris::IndexMat > const & aNodeInds,
                                moris::Matrix< moris::IndexMat > const & aParentNodeInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentNodeRanks,
                                moris::Matrix< moris::IndexMat > const & aParentEdgeInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentEdgeRanks,
                                moris::Matrix< moris::IndexMat > const & aParentFaceInds,
                                moris::Matrix< moris::DDSTMat >  const & aParentFaceRanks,
                                enum TemplateType                        aTemplateType,
                                moris::size_t                            aPermutationId = 0)
    {
        MORIS_ASSERT(aTemplateType == TemplateType::REGULAR_SUBDIVISION_QUAD4 || aTemplateType ==TemplateType::TRI_3,
                     "Entered node inheritance constructor for template without node inheritance. Currently, only the regular subdivision quad4 and tri3 templates have node inheritance in 2D.");

        mElemIndToReplace = aElemToReplace;
        mParentElemInd    = aParentElemInd;
        mNodeInds         = aNodeInds.copy();
        mParentNodeInds   = aParentNodeInds.copy();
        mParentNodeRanks  = aParentNodeRanks.copy();
        mParentEdgeInds   = aParentEdgeInds.copy();
        mParentEdgeRanks  = aParentEdgeRanks.copy();
        mParentFaceInds   = aParentFaceInds.copy();
        mParentFaceRanks  = aParentFaceRanks.copy();
        template_catalog(aPermutationId,aTemplateType);
    }*/

    // Template type
    enum TemplateType mTemplateType;

    // Number of elements to replace and number of new elements in the template
    moris::size_t                    mNumNewElem;
    moris::size_t                    mNumElemToReplace;
    moris::size_t                    mElemIndToReplace;

    // Parent Entity's parent information
    // This is the information relative to the parent this template is created from
    moris::moris_index               mParentElemInd;
    moris::Matrix< moris::IndexMat > mParentNodeInds;
    moris::Matrix< moris::DDSTMat  > mParentNodeRanks;
    moris::Matrix< moris::IndexMat > mParentEdgeInds;
    moris::Matrix< moris::DDSTMat  > mParentEdgeRanks;
    moris::Matrix< moris::IndexMat > mParentFaceInds;
    moris::Matrix< moris::DDSTMat  > mParentFaceRanks;

    // spatial dimension of template
    moris::uint mSpatialDimension;

    // topology of children elements
    enum CellTopology mElementTopology;

    // Node indices in the template
    moris::Matrix< moris::IndexMat > mNodeInds;

    // Element to Node Connectivity
    moris::Matrix< moris::IndexMat > mNewElementToNode;

    // Parent's of an entity ordered by ordinal relative to
    // (these are the inheritance of the new elements created by this template)
    moris::Matrix< moris::DDSTMat >  mNewParentEdgeRanks;
    moris::Matrix< moris::IndexMat > mNewParentEdgeOrdinals;
    moris::Matrix< moris::DDSTMat >  mNewParentFaceRanks;
    moris::Matrix< moris::IndexMat > mNewParentFaceOrdinals;
    moris::Matrix< moris::DDSTMat >  mNewElementInterfaceSides;

    // Node inheritance information
    bool mHasNodeInheritance = false;
    moris::Matrix< moris::DDSTMat >  mNewNodeParentRanks;
    moris::Matrix< moris::IndexMat > mNewNodeParentOrdinals;

    // Reindex flag (meaning this template has been reindexed by the Child Mesh)
    bool mIsReindexed = false;
private:

    // Note everything below this line is just template data and selecting the correct template
    void
    template_catalog(moris::size_t const & aPermutationId,
                     enum TemplateType aTemplateType)
    {
        switch(aTemplateType)
        {
        case(TemplateType::REGULAR_SUBDIVISION_HEX8):
        {
            hex_8_reg_sub_template();
            break;
        }
        case(TemplateType::TET_4):
         {
             tet4_template();
             break;
         }
        case(TemplateType::REGULAR_SUBDIVISION_QUAD4):
        {
           	quad_4_reg_sub_template();
           	break;
        }
           case(TemplateType::TRI_3):
        {
           	tri3_template();
           	break;
        }
        case(TemplateType::CONFORMAL_TRI3):
        {
            conformal_tri_template(aPermutationId);
            break;
        }

        case(TemplateType::HIERARCHY_TET4_3N):
         {
            hierarchy_tet4_3N(aPermutationId);
             break;
         }
        case(TemplateType::HIERARCHY_TET4_4Na):
         {
            hierarchy_tet4_4na(aPermutationId);
             break;
         }
        case(TemplateType::HIERARCHY_TET4_4Nb):
         {
            hierarchy_tet4_4nb(aPermutationId);
             break;
         }
        case(TemplateType::HIERARCHY_TET4_4Nc):
         {
            hierarchy_tet4_4nc(aPermutationId);
             break;
         }
        case(TemplateType::BISECTED_TET4):
        {
            bisected_tet(aPermutationId);
            break;
        }
        case(TemplateType::HIERARCHY_TET4_2):
        {
            hierarchy_tet4_2(aPermutationId);
            break;
        }
        default :
            std::cout<<"Template not found in the catalog"<<std::endl;
            break;
        }
    }

    void
    hex_8_reg_sub_template()
    {
        MORIS_ASSERT(mNodeInds.n_cols() == 15, "For a Hex8 regular subdivision template, there must be 15 node inds.");

        mNewElementToNode = moris::Matrix< moris::IndexMat >({
                                                      {0, 8, 1,  14},
                                                      {1, 8, 5,  14},
                                                      {4, 5, 8,  14},
                                                      {0, 4, 8,  14},
                                                      {1, 9, 2,  14},
                                                      {2, 9, 6,  14},
                                                      {5, 6, 9,  14},
                                                      {1, 5, 9,  14},
                                                      {2, 10, 3, 14},
                                                      {2, 6, 10, 14},
                                                      {6, 7, 10, 14},
                                                      {3, 10, 7, 14},
                                                      {0, 3, 11, 14},
                                                      {3, 7, 11, 14},
                                                      {4, 11, 7, 14},
                                                      {0, 11, 4, 14},
                                                      {0, 1, 12, 14},
                                                      {1, 2, 12, 14},
                                                      {2, 3, 12, 14},
                                                      {0, 12, 3, 14},
                                                      {4, 13, 5, 14},
                                                      {5, 13, 6, 14},
                                                      {6, 13, 7, 14},
                                                      {4, 7, 13, 14}});

        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;
        mNumNewElem = 24;
        mNumElemToReplace = 1;
        mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {1, 2, 2, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {2, 2, 1, 3, 3, 3}, {1, 2, 2, 3, 3, 3}});
        mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 0, 0, 0, 0, 0}, {0, 0, 9, 0, 0, 0}, {4, 0, 0, 0, 0, 0}, {8, 0, 0, 0, 0, 0}, {1, 1, 1, 0, 0, 0}, {1, 1, 10, 0, 0, 0}, {5, 1, 1, 0, 0, 0}, {9, 1, 1, 0, 0, 0}, {2, 2, 2, 0, 0, 0}, {10, 2, 2, 0, 0, 0}, {6, 2, 2, 0, 0, 0}, {2, 2, 11, 0, 0, 0}, {3, 3, 3, 0, 0, 0}, {11, 3, 3, 0, 0, 0}, {3, 3, 7, 0, 0, 0}, {3, 3, 8, 0, 0, 0}, {0, 4, 4, 0, 0, 0}, {1, 4, 4, 0, 0, 0}, {2, 4, 4, 0, 0, 0}, {4, 4, 3, 0, 0, 0}, {5, 5, 4, 0, 0, 0}, {5, 5, 5, 0, 0, 0}, {5, 5, 6, 0, 0, 0}, {7, 5, 5, 0, 0, 0}});
        mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}, {3, 3, 3, 2}});
        mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 4}, {0, 0, 0, 4}, {0, 0, 0, 4}, {0, 0, 0, 4}, {0, 0, 0, 5}, {0, 0, 0, 5}, {0, 0, 0, 5}, {0, 0, 0, 5}});
        mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >(24,1,std::numeric_limits<moris::size_t>::max());

        // Mark this template as having node inheritance
        mHasNodeInheritance    = true;
        mNewNodeParentRanks    = moris::Matrix< moris::DDSTMat >({{0,0,0,0,0,0,0,0,2,2,2,2,2,2,3}});
        mNewNodeParentOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3,4,5,6,7,0,1,2,3,4,5,0}});
    }

    void
    tet4_template()
    {
        mSpatialDimension = 3;
        mElementTopology          = CellTopology::TET4;
        mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,3}});
        mNumNewElem = 1;
        mNumElemToReplace = 0;
        mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 1, 1, 1, 1, 1}});
        mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 3, 4, 5}});
        mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{2, 2, 2, 2}});
        mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 3}});
        mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >(1,1,std::numeric_limits<moris::size_t>::max());

    }

    void
	quad_4_reg_sub_template()
    {
        MORIS_ASSERT(mNodeInds.n_cols() == 5, "For a Quad4 regular subdivision template, there must be 5 node inds.");
        mSpatialDimension = 2;
        mElementTopology          = CellTopology::TRI3;
        mNewElementToNode = {{0, 1, 4},
                             {1, 2, 4},
                             {2, 3, 4},
                             {3, 0, 4}};

        mNumNewElem 	  = 4;
        mNumElemToReplace = 1;
        mNewParentEdgeRanks    = {{1, 3, 3},
                                  {1, 3, 3},
                                  {1, 3, 3},
                                  {1, 3, 3}};
        mNewParentEdgeOrdinals = {{0, 0, 0},
                                  {1, 0, 0},
                                  {2, 0, 0},
                                  {3, 0, 0}};
        mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >(4, 1, std::numeric_limits<moris::size_t>::max());

        // Mark template as having node inheritance
        mHasNodeInheritance    = true;
        mNewNodeParentRanks    = {{0, 0, 0, 0, 3}};
        mNewNodeParentOrdinals = {{0, 1, 2, 3, 0}};

        mElementTopology = CellTopology::TRI3;
    }

    void
	tri3_template()
    {
        mSpatialDimension = 2;
        mElementTopology          = CellTopology::TRI3;
        mNewElementToNode 	  = {{0, 1, 2}};
        mNumNewElem 		  = 1;
        mNumElemToReplace 	  = 0;
        mNewParentEdgeRanks 	  = {{1, 1, 1}};
        mNewParentEdgeOrdinals 	  = {{0, 1, 2}};
        mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >(1,1,std::numeric_limits<moris::size_t>::max());
    }

    void
    conformal_tri_template(moris::size_t const & aPermutation)
    {

        mSpatialDimension = 2;
        moris::size_t tMax = std::numeric_limits<moris::size_t>::max();

        switch(aPermutation)
        {
            case 1:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0, 3, 2},{3,1,4},{3,4,2}};
                mNumNewElem               = 3;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1, 3, 1},{1, 1, 3},{3, 1, 3}};
                mNewParentEdgeOrdinals    = {{0, 0, 2},{0, 1, 0},{0, 1, 0}};
                mNewElementInterfaceSides = {{tMax},{2},{0}};
                break;
            }
            case 2:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0,3,4},{3,2,4},{3,1,2}};
                mNumNewElem               = 3;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1,3,1},{3,1,3},{1,1,3}};
                mNewParentEdgeOrdinals    = {{0,0,2},{0,2,0},{0,1,0}};
                mNewElementInterfaceSides = {{1},{2},{tMax}};
                break;
            }
            case 3:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0,1,3},{0,3,4},{3,2,4}};
                mNumNewElem               = 3;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1,1,3},{3,3,1},{1,1,3}};
                mNewParentEdgeOrdinals    = {{0,1,0},{0,0,2},{1,2,0}};
                mNewElementInterfaceSides = {{tMax},{1},{2}};
                break;
            }
            case 10:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0,3,2},{3,1,2}};
                mNumNewElem               = 2;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1,3,1},{1,1,3}};
                mNewParentEdgeOrdinals    = {{0,0,2},{0,1,0}};
                mNewElementInterfaceSides = {{1},{2}};
                break;
            }
            case 11:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0,1,3},{0,3,2}};
                mNumNewElem               = 2;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1,1,3},{3,1,1}};
                mNewParentEdgeOrdinals    = {{0,1,0},{0,1,2}};
                mNewElementInterfaceSides = {{2},{0}};
                break;
            }
            case 12:
            {
                mElementTopology          = CellTopology::TRI3;
                mNewElementToNode         = {{0,1,3},{1,2,3}};
                mNumNewElem               = 2;
                mNumElemToReplace         = 1;
                mNewParentEdgeRanks       = {{1,3,1},{1,1,3}};
                mNewParentEdgeOrdinals    = {{0,0,2},{1,2,0}};
                mNewElementInterfaceSides = {{1},{2}};
                break;
            }
            default:
                MORIS_ERROR(0,"Invalid conformal tri 3 permutation");
                break;
        }
    }

    void
    hierarchy_tet4_3N(moris::size_t const & aPermutation)
    {

        mSpatialDimension = 3;
        mElementTopology          = CellTopology::TET4;
        switch(aPermutation)
        {
        case(320):
        {
            // Permutation 320
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 2, 3, 0, 2}, {0, 3, 3, 0, 0, 2}, {1, 2, 3, 0, 2, 2}, {4, 5, 1, 0, 3, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 2, 3}, {0, 0, 0, 3}, {0, 2, 0, 3}, {0, 2, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }

        case(32):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 3, 0, 3, 0}, {2, 2, 2, 3, 3, 0}, {5, 3, 2, 3, 0, 0}, {1, 4, 5, 3, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {3, 0, 0, 2}, {0, 0, 0, 2}, {3, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});

            break;
        }

        case(203):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 2, 3}, {3, 0, 0, 2, 2, 3}, {4, 0, 0, 2, 3, 3}, {5, 1, 4, 2, 2, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 3, 0}, {2, 0, 0, 0}, {0, 3, 0, 0}, {2, 3, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(251):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 5, 2, 3, 2}, {1, 1, 1, 3, 3, 2}, {4, 5, 1, 3, 2, 2}, {0, 3, 4, 3, 2, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 2, 1}, {3, 0, 0, 1}, {0, 2, 0, 1}, {3, 2, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }

        case(512):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 1, 5, 2, 1}, {2, 3, 3, 2, 2, 1}, {0, 1, 3, 2, 1, 1}, {3, 4, 0, 2, 5, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 1, 3}, {2, 0, 0, 3}, {0, 1, 0, 3}, {2, 1, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(125):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 1, 3}, {5, 2, 2, 1, 1, 3}, {3, 2, 2, 1, 3, 3}, {4, 0, 3, 1, 1, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 3, 2}, {1, 0, 0, 2}, {0, 3, 0, 2}, {1, 3, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(140):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 4, 1, 3, 1}, {0, 0, 0, 3, 3, 1}, {3, 4, 0, 3, 1, 1}, {2, 5, 3, 3, 1, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 1, 0}, {3, 0, 0, 0}, {0, 1, 0, 0}, {3, 1, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});

            break;
        }

        case(401):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 0, 4, 1, 0}, {1, 3, 3, 1, 1, 0}, {2, 0, 3, 1, 0, 0}, {5, 3, 2, 1, 4, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {1, 0, 0, 3}, {0, 0, 0, 3}, {1, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});

            break;
        }
        case(14):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 0, 3}, {4, 1, 1, 0, 0, 3}, {5, 1, 1, 0, 3, 3}, {3, 2, 5, 0, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 3, 1}, {0, 0, 0, 1}, {0, 3, 0, 1}, {0, 3, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(453):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 5, 4, 0, 1}, {3, 2, 2, 0, 0, 1}, {2, 5, 2, 0, 1, 1}, {0, 1, 2, 0, 4, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 1, 2}, {0, 0, 0, 2}, {0, 1, 0, 2}, {0, 1, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(534):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 3, 5, 1, 2}, {4, 0, 0, 1, 1, 2}, {0, 3, 0, 1, 2, 2}, {1, 2, 0, 1, 5, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 2, 0}, {1, 0, 0, 0}, {0, 2, 0, 0}, {1, 2, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(345):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 6}, {4, 1, 5, 6}, {1, 2, 5, 6}, {1, 3, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 2, 2, 2}, {1, 1, 2, 2, 2, 2}, {1, 1, 1, 2, 1, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 4, 3, 2, 0}, {5, 1, 1, 2, 2, 0}, {1, 4, 1, 2, 0, 0}, {2, 0, 1, 2, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 2, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {2, 0, 0, 1}, {0, 0, 0, 1}, {2, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});

            break;
        }

        case(230):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 0, 3, 2, 2}, {0, 0, 0, 3, 3, 2}, {4, 0, 3, 2, 3, 2}, {1, 4, 5, 2, 3, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 0, 0}, {3, 0, 0, 0}, {0, 0, 2, 0}, {3, 0, 2, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(302):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 3, 2, 3, 0}, {2, 3, 3, 2, 2, 0}, {1, 3, 0, 0, 2, 0}, {5, 1, 4, 3, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {2, 0, 0, 3}, {0, 0, 0, 3}, {2, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(23):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 0, 0, 3}, {3, 2, 2, 0, 0, 3}, {5, 2, 2, 3, 0, 3}, {4, 5, 1, 0, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 0, 2}, {0, 0, 0, 2}, {0, 0, 3, 2}, {0, 0, 3, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(521):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 2, 3, 1, 5, 2}, {1, 3, 3, 1, 1, 2}, {0, 3, 2, 2, 1, 2}, {4, 0, 3, 5, 1, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 2, 0, 3}, {1, 0, 0, 3}, {0, 0, 2, 3}, {1, 0, 2, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(152):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 5, 2, 3, 1, 1}, {2, 2, 2, 3, 3, 1}, {3, 2, 5, 1, 3, 1}, {0, 3, 4, 1, 3, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 1, 0, 2}, {3, 0, 0, 2}, {0, 0, 1, 2}, {3, 0, 1, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(215):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 2, 2, 3}, {5, 1, 1, 2, 2, 3}, {4, 1, 1, 3, 2, 3}, {3, 4, 0, 2, 2, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 0, 1}, {2, 0, 0, 1}, {0, 0, 3, 1}, {2, 0, 3, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(410):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 1, 3, 0, 4, 1}, {0, 3, 3, 0, 0, 1}, {2, 3, 1, 1, 0, 1}, {3, 2, 5, 4, 0, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 1, 0, 3}, {0, 0, 0, 3}, {0, 0, 1, 3}, {0, 0, 1, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
           break;
        }
        case(41 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 4, 1, 3, 0, 0}, {1, 1, 1, 3, 3, 0}, {5, 1, 4, 0, 3, 0}, {2, 5, 3, 0, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {3, 0, 0, 1}, {0, 0, 0, 1}, {3, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(104):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 1, 1, 3}, {4, 0, 0, 1, 1, 3}, {3, 0, 0, 3, 1, 3}, {5, 3, 2, 1, 1, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 0, 0}, {1, 0, 0, 0}, {0, 0, 3, 0}, {1, 0, 3, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(543):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 4, 0, 2, 5, 1}, {3, 0, 0, 2, 2, 1}, {0, 0, 4, 1, 2, 1}, {2, 0, 1, 5, 2, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 1, 0, 0}, {2, 0, 0, 0}, {0, 0, 1, 0}, {2, 0, 1, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(354):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 5, 1, 0, 3, 2}, {4, 1, 1, 0, 0, 2}, {1, 1, 5, 2, 0, 2}, {0, 1, 2, 3, 0, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 2, 0, 1}, {0, 0, 0, 1}, {0, 0, 2, 1}, {0, 0, 2, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }
        case(435):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 6}, {1, 4, 5, 6}, {2, 1, 5, 6}, {3, 1, 2, 6}});
            mNumNewElem = 4;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 2, 2}, {1, 2, 1, 2, 2, 2}, {1, 1, 1, 1, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 3, 2, 1, 4, 0}, {5, 2, 2, 1, 1, 0}, {2, 2, 3, 0, 1, 0}, {1, 2, 0, 4, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 3, 2, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {1, 0, 0, 2}, {0, 0, 0, 2}, {1, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {std::numeric_limits<moris::size_t>::max()}});
            break;
        }

        default:
        {
            std::cout<<"WARNING INVALID PERMUTATION"<<std::endl;
            break;
        }
        }
    }

    void
    hierarchy_tet4_4na(moris::size_t const & aPermutationId)
    {
        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;
        switch(aPermutationId)
        {
         case(5420):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 5, 4, 0, 1}, {0, 0, 0, 0, 2, 1}, {0, 2, 0, 3, 2, 2}, {1, 1, 5, 2, 3, 2}, {0, 0, 1, 3, 3, 2}, {0, 1, 0, 0, 4, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 2, 1, 2}, {2, 2, 2, 0}, {3, 2, 2, 2}, {3, 2, 2, 1}, {3, 2, 2, 2}, {0, 1, 2, 2}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 5, 4, 0, 1}, {0, 0, 0, 0, 2, 1}, {0, 2, 0, 3, 2, 2}, {1, 1, 5, 2, 3, 2}, {0, 0, 1, 3, 3, 2}, {0, 1, 0, 0, 4, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 1, 2}, {0, 0, 0, 0}, {3, 2, 0, 0}, {3, 0, 2, 1}, {3, 0, 0, 0}, {0, 1, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             break;
         }
         case(5240):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 5, 2, 3, 2}, {0, 3, 3, 0, 1, 2}, {0, 1, 0, 0, 4, 1}, {3, 2, 5, 4, 0, 1}, {0, 0, 2, 0, 0, 1}, {0, 2, 0, 3, 2, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 2, 1}, {0, 0, 0, 3}, {0, 1, 0, 0}, {0, 0, 1, 2}, {0, 0, 0, 0}, {3, 2, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             break;
         }
         case(425 ):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 0, 4, 1, 0}, {5, 1, 1, 0, 3, 0}, {5, 3, 0, 2, 2, 3}, {3, 0, 0, 2, 2, 3}, {5, 0, 0, 2, 2, 3}, {5, 0, 0, 1, 4, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {0, 0, 0, 1}, {2, 3, 0, 0}, {2, 0, 3, 0}, {2, 0, 0, 0}, {1, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(3501):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 3, 5, 1, 2}, {1, 1, 1, 0, 0, 2}, {1, 0, 0, 3, 0, 0}, {2, 2, 3, 0, 3, 0}, {1, 0, 2, 3, 3, 0}, {1, 2, 0, 1, 5, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 2, 0}, {0, 0, 0, 1}, {3, 0, 0, 0}, {3, 0, 0, 2}, {3, 0, 0, 0}, {1, 2, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(1503):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 1, 5, 2, 1}, {3, 2, 2, 0, 3, 1}, {3, 3, 0, 0, 0, 3}, {4, 1, 1, 0, 0, 3}, {3, 0, 1, 0, 0, 3}, {3, 1, 0, 2, 5, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 1, 3}, {0, 0, 0, 2}, {0, 3, 0, 0}, {0, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             break;
         }
         case(3051):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 3, 0, 3, 0}, {1, 3, 3, 0, 2, 0}, {1, 2, 0, 1, 5, 2}, {4, 0, 3, 5, 1, 2}, {1, 0, 0, 1, 1, 2}, {1, 0, 0, 3, 0, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {0, 0, 0, 3}, {1, 2, 0, 0}, {1, 0, 2, 0}, {1, 0, 0, 0}, {3, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(1053):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 0, 3}, {3, 0, 0, 0, 1, 3}, {3, 1, 0, 2, 5, 1}, {2, 3, 1, 5, 2, 1}, {3, 0, 3, 2, 2, 1}, {3, 3, 0, 0, 0, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 3, 1}, {0, 0, 0, 0}, {2, 1, 0, 0}, {2, 0, 1, 3}, {2, 0, 0, 0}, {0, 3, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             break;
         }
         case(4312):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 4, 3, 2, 0}, {2, 2, 2, 0, 1, 0}, {2, 1, 0, 3, 1, 1}, {0, 0, 4, 1, 3, 1}, {2, 0, 0, 3, 3, 1}, {2, 0, 0, 2, 3, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {0, 0, 0, 2}, {3, 1, 0, 0}, {3, 0, 1, 0}, {3, 0, 0, 0}, {2, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(2314):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 2, 3, 0, 2}, {4, 0, 0, 0, 3, 2}, {4, 3, 0, 1, 1, 3}, {5, 2, 2, 1, 1, 3}, {4, 0, 2, 1, 1, 3}, {4, 2, 0, 0, 3, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 2, 3}, {0, 0, 0, 0}, {1, 3, 0, 0}, {1, 0, 3, 2}, {1, 0, 0, 0}, {0, 2, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(4132):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 4, 1, 3, 1}, {2, 3, 3, 0, 0, 1}, {2, 0, 0, 2, 3, 0}, {5, 1, 4, 3, 2, 0}, {2, 0, 1, 2, 2, 0}, {2, 1, 0, 3, 1, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 1, 0}, {0, 0, 0, 3}, {2, 0, 0, 0}, {2, 0, 0, 1}, {2, 0, 0, 0}, {3, 1, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(2134):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 1, 3}, {4, 1, 1, 0, 2, 3}, {4, 2, 0, 0, 3, 2}, {0, 3, 2, 3, 0, 2}, {4, 0, 3, 0, 0, 2}, {4, 3, 0, 1, 1, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 3, 2}, {0, 0, 0, 1}, {0, 2, 0, 0}, {0, 0, 2, 3}, {0, 0, 0, 0}, {1, 3, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});

             break;
         }
         case(245 ):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 3, 7, 6}, {4, 3, 6, 7}, {4, 3, 7, 5}, {1, 2, 7, 5}, {2, 4, 7, 5}, {4, 2, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 2}, {1, 2, 2, 3, 2, 2}, {1, 2, 3, 2, 1, 2}, {1, 2, 1, 1, 2, 2}, {1, 3, 2, 2, 2, 2}, {1, 2, 3, 2, 1, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 2, 3}, {5, 2, 2, 0, 0, 3}, {5, 0, 0, 1, 4, 0}, {1, 3, 0, 4, 1, 0}, {5, 0, 3, 1, 1, 0}, {5, 3, 0, 2, 2, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2}, {3, 3, 3, 2}, {2, 2, 3, 3}, {2, 3, 2, 2}, {2, 3, 3, 3}, {2, 2, 3, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 3, 0}, {0, 0, 0, 2}, {1, 0, 0, 0}, {1, 0, 0, 3}, {1, 0, 0, 0}, {2, 3, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {2}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}});
             break;
         }
         case(4502):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 4, 0, 2, 5, 1}, {2, 2, 2, 0, 0, 1}, {2, 0, 0, 0, 3, 0}, {1, 4, 1, 3, 0, 0}, {2, 1, 0, 3, 3, 0}, {2, 0, 1, 5, 2, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 1, 0, 0}, {0, 0, 0, 2}, {3, 0, 0, 0}, {3, 0, 0, 1}, {3, 0, 0, 0}, {2, 0, 1, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(4052):
        {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 4, 1, 3, 0, 0}, {2, 3, 3, 1, 0, 0}, {2, 0, 1, 5, 2, 1}, {3, 4, 0, 2, 5, 1}, {2, 0, 0, 2, 2, 1}, {2, 0, 0, 0, 3, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {0, 0, 0, 3}, {2, 0, 1, 0}, {2, 1, 0, 0}, {2, 0, 0, 0}, {3, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
        }
         case(1243):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 2, 2, 3}, {3, 2, 2, 1, 0, 3}, {3, 0, 1, 4, 0, 1}, {0, 1, 3, 0, 4, 1}, {3, 3, 0, 0, 0, 1}, {3, 0, 3, 2, 2, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 0, 1}, {0, 0, 0, 2}, {0, 0, 1, 0}, {0, 1, 0, 3}, {0, 0, 0, 0}, {2, 0, 3, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(2504):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 2, 3, 1, 5, 2}, {4, 1, 1, 3, 0, 2}, {4, 0, 3, 0, 0, 3}, {3, 2, 2, 0, 0, 3}, {4, 2, 0, 0, 0, 3}, {4, 0, 2, 5, 1, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 2, 0, 3}, {0, 0, 0, 1}, {0, 0, 3, 0}, {0, 3, 0, 2}, {0, 0, 0, 0}, {1, 0, 2, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(2054):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 0, 0, 3}, {4, 0, 0, 2, 0, 3}, {4, 0, 2, 5, 1, 2}, {1, 2, 3, 1, 5, 2}, {4, 3, 0, 1, 1, 2}, {4, 0, 3, 0, 0, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 0, 2}, {0, 0, 0, 0}, {1, 0, 2, 0}, {1, 2, 0, 3}, {1, 0, 0, 0}, {0, 0, 3, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});

             break;
         }
         case(5310):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 5, 1, 0, 3, 2}, {0, 0, 0, 1, 0, 2}, {0, 0, 1, 1, 3, 1}, {2, 5, 2, 3, 1, 1}, {0, 2, 0, 3, 3, 1}, {0, 0, 2, 3, 0, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 2, 0, 1}, {0, 0, 0, 0}, {3, 0, 1, 0}, {3, 1, 0, 2}, {3, 0, 0, 0}, {0, 0, 2, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(5130):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 5, 2, 3, 1, 1}, {0, 3, 3, 2, 0, 1}, {0, 0, 2, 3, 0, 2}, {4, 5, 1, 0, 3, 2}, {0, 1, 0, 0, 0, 2}, {0, 0, 1, 1, 3, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 1, 0, 2}, {0, 0, 0, 3}, {0, 0, 2, 0}, {0, 2, 0, 1}, {0, 0, 0, 0}, {3, 0, 1, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});

             break;
         }
         case(135 ):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 1, 1, 3}, {5, 1, 1, 0, 0, 3}, {5, 0, 0, 3, 2, 0}, {2, 0, 3, 2, 3, 0}, {5, 3, 0, 2, 2, 0}, {5, 0, 3, 1, 1, 3}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 0, 0}, {0, 0, 0, 1}, {2, 0, 0, 0}, {2, 0, 0, 3}, {2, 0, 0, 0}, {1, 0, 3, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});

             break;
         }
         case(315 ):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 3, 2, 3, 0}, {5, 2, 2, 3, 0, 0}, {5, 0, 3, 1, 1, 3}, {4, 0, 0, 1, 1, 3}, {5, 0, 0, 1, 1, 3}, {5, 0, 0, 3, 2, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {0, 0, 0, 2}, {1, 0, 3, 0}, {1, 3, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});

             break;
         }
         case(3421):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 3, 2, 1, 4, 0}, {1, 1, 1, 2, 0, 0}, {1, 0, 2, 2, 3, 2}, {0, 3, 0, 3, 2, 2}, {1, 0, 0, 3, 3, 2}, {1, 0, 0, 4, 1, 0}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {0, 0, 0, 1}, {3, 0, 2, 0}, {3, 2, 0, 0}, {3, 0, 0, 0}, {1, 0, 0, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(3241):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 0, 3, 2, 2}, {1, 3, 3, 0, 0, 2}, {1, 0, 0, 4, 1, 0}, {5, 3, 2, 1, 4, 0}, {1, 2, 0, 1, 1, 0}, {1, 0, 2, 2, 3, 2}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 0, 0}, {0, 0, 0, 3}, {1, 0, 0, 0}, {1, 0, 0, 2}, {1, 0, 0, 0}, {3, 0, 2, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         case(1423):
         {
             mNewElementToNode = moris::Matrix< moris::IndexMat >({{3, 0, 7, 6}, {3, 4, 6, 7}, {3, 4, 7, 5}, {2, 1, 7, 5}, {4, 2, 7, 5}, {2, 4, 7, 6}});
             mNumNewElem = 6;
             mNumElemToReplace = 1;
             mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 2, 1, 2}, {1, 2, 2, 2, 3, 2}, {1, 3, 2, 1, 2, 2}, {1, 1, 2, 2, 1, 2}, {1, 2, 3, 2, 2, 2}, {1, 3, 2, 1, 2, 2}});
             mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 1, 3, 0, 4, 1}, {3, 0, 0, 3, 0, 1}, {3, 0, 3, 2, 2, 3}, {5, 1, 1, 2, 2, 3}, {3, 1, 0, 2, 2, 3}, {3, 0, 1, 4, 0, 1}});
             mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2}, {3, 3, 3, 2}, {2, 3, 2, 3}, {2, 2, 3, 2}, {2, 3, 3, 3}, {2, 3, 2, 3}});
             mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 1, 0, 3}, {0, 0, 0, 0}, {2, 0, 3, 0}, {2, 3, 0, 1}, {2, 0, 0, 0}, {0, 0, 1, 0}});
             mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {1}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}});
             break;
         }
         default :
             std::cout<<"Template not found in the catalog"<<std::endl;
             break;

        }

    }

    void
    hierarchy_tet4_4nb(moris::size_t const & aPermutationId)
    {
        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;
        switch(aPermutationId)
        {
        case(4250):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 3, 0, 4, 0}, {1, 2, 3, 4, 1, 0}, {5, 2, 2, 1, 1, 0}, {5, 2, 2, 1, 4, 0}, {3, 2, 2, 4, 0, 0}, {0, 3, 2, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3}, {1, 0, 0, 3}, {1, 0, 0, 2}, {1, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(2450):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3, 2, 0}, {3, 4, 0, 2, 2, 0}, {5, 1, 4, 2, 2, 0}, {5, 1, 1, 2, 2, 0}, {1, 4, 1, 2, 3, 0}, {0, 0, 4, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 0}, {2, 0, 0, 0}, {2, 0, 0, 1}, {2, 0, 0, 1}, {3, 0, 0, 1}, {3, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});

            break;
        }
        case(4205):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 4, 0}, {3, 2, 2, 4, 0, 0}, {0, 3, 2, 0, 0, 0}, {0, 3, 3, 0, 4, 0}, {1, 2, 3, 4, 1, 0}, {5, 2, 2, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 3}, {0, 0, 0, 3}, {1, 0, 0, 3}, {1, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(2405):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 2, 2, 0}, {1, 4, 1, 2, 3, 0}, {0, 0, 4, 3, 3, 0}, {0, 0, 0, 3, 2, 0}, {3, 4, 0, 2, 2, 0}, {5, 1, 4, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {3, 0, 0, 1}, {3, 0, 0, 0}, {3, 0, 0, 0}, {2, 0, 0, 0}, {2, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(5031):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 3, 1, 5, 0}, {2, 0, 3, 5, 2, 0}, {3, 0, 0, 2, 2, 0}, {3, 0, 0, 2, 5, 0}, {4, 0, 0, 5, 1, 0}, {1, 3, 0, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {2, 0, 0, 3}, {2, 0, 0, 0}, {2, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(5013):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 5, 0}, {4, 0, 0, 5, 1, 0}, {1, 3, 0, 1, 1, 0}, {1, 3, 3, 1, 5, 0}, {2, 0, 3, 5, 2, 0}, {3, 0, 0, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 3}, {1, 0, 0, 3}, {2, 0, 0, 3}, {2, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(531 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 1, 3, 0, 0}, {4, 5, 1, 0, 0, 0}, {3, 2, 5, 0, 0, 0}, {3, 2, 2, 0, 0, 0}, {2, 5, 2, 0, 3, 0}, {1, 1, 5, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 0, 2}, {3, 0, 0, 2}, {3, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(513 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 0, 0, 0}, {2, 5, 2, 0, 3, 0}, {1, 1, 5, 3, 3, 0}, {1, 1, 1, 3, 0, 0}, {4, 5, 1, 0, 0, 0}, {3, 2, 5, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 2}, {3, 0, 0, 2}, {3, 0, 0, 1}, {3, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});
            break;
        }
        case(3124):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 3, 0}, {5, 1, 1, 3, 2, 0}, {2, 3, 1, 2, 2, 0}, {2, 3, 3, 2, 3, 0}, {0, 1, 3, 3, 0, 0}, {4, 1, 1, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1}, {2, 0, 0, 1}, {2, 0, 0, 3}, {2, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});

            break;
        }
        case(1342):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 2, 3, 1, 0}, {5, 3, 2, 1, 1, 0}, {4, 0, 3, 1, 1, 0}, {4, 0, 0, 1, 1, 0}, {0, 3, 0, 1, 3, 0}, {2, 2, 3, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {1, 0, 0, 2}, {1, 0, 0, 0}, {1, 0, 0, 0}, {3, 0, 0, 0}, {3, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});

            break;
        }
        case(1324):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 1, 1, 0}, {0, 3, 0, 1, 3, 0}, {2, 2, 3, 3, 3, 0}, {2, 2, 2, 3, 1, 0}, {5, 3, 2, 1, 1, 0}, {4, 0, 3, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 0}, {3, 0, 0, 0}, {3, 0, 0, 2}, {3, 0, 0, 2}, {1, 0, 0, 2}, {1, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});

            break;
        }
        case(3142):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 6, 7}, {0, 3, 6, 7}, {3, 5, 6, 7}, {5, 1, 6, 7}, {1, 2, 6, 7}, {2, 4, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 2, 1, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 3, 2, 3, 0}, {0, 1, 3, 3, 0, 0}, {4, 1, 1, 0, 0, 0}, {4, 1, 1, 0, 3, 0}, {5, 1, 1, 3, 2, 0}, {2, 3, 1, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 1}, {0, 0, 0, 1}, {2, 0, 0, 1}, {2, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}});

            break;
        }
        case(5024):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 5, 1, 0}, {3, 0, 0, 2, 5, 0}, {2, 0, 3, 2, 2, 0}, {2, 3, 3, 5, 2, 0}, {1, 3, 0, 1, 5, 0}, {4, 0, 0, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 0}, {2, 0, 0, 0}, {2, 0, 0, 3}, {2, 0, 0, 3}, {1, 0, 0, 3}, {1, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        case(524 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 0, 0}, {1, 1, 5, 3, 0, 0}, {2, 5, 2, 3, 3, 0}, {2, 2, 2, 0, 3, 0}, {3, 2, 5, 0, 0, 0}, {4, 5, 1, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1}, {3, 0, 0, 1}, {3, 0, 0, 2}, {3, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(5042):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 3, 5, 2, 0}, {1, 3, 0, 1, 5, 0}, {4, 0, 0, 1, 1, 0}, {4, 0, 0, 5, 1, 0}, {3, 0, 0, 2, 5, 0}, {2, 0, 3, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {1, 0, 0, 3}, {1, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}, {2, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        case(542 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 2, 0, 3, 0}, {3, 2, 5, 0, 0, 0}, {4, 5, 1, 0, 0, 0}, {4, 1, 1, 0, 0, 0}, {1, 1, 5, 3, 0, 0}, {2, 5, 2, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 1}, {0, 0, 0, 1}, {3, 0, 0, 1}, {3, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        case(3150):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 3, 3, 0, 0}, {2, 3, 1, 2, 3, 0}, {5, 1, 1, 2, 2, 0}, {5, 1, 1, 3, 2, 0}, {4, 1, 1, 0, 3, 0}, {0, 1, 3, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3}, {2, 0, 0, 3}, {2, 0, 0, 1}, {2, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        case(1350):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1, 3, 0}, {4, 0, 3, 1, 1, 0}, {5, 3, 2, 1, 1, 0}, {5, 2, 2, 1, 1, 0}, {2, 2, 3, 3, 1, 0}, {0, 3, 0, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 2}, {1, 0, 0, 2}, {3, 0, 0, 2}, {3, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        case(3105):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 3, 2, 0}, {4, 1, 1, 0, 3, 0}, {0, 1, 3, 0, 0, 0}, {0, 3, 3, 3, 0, 0}, {2, 3, 1, 2, 3, 0}, {5, 1, 1, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 3}, {0, 0, 0, 3}, {2, 0, 0, 3}, {2, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(1305):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 1, 0}, {2, 2, 3, 3, 1, 0}, {0, 3, 0, 3, 3, 0}, {0, 0, 0, 1, 3, 0}, {4, 0, 3, 1, 1, 0}, {5, 3, 2, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {3, 0, 0, 2}, {3, 0, 0, 0}, {3, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(4231):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 3, 4, 1, 0}, {0, 3, 2, 0, 4, 0}, {3, 2, 2, 0, 0, 0}, {3, 2, 2, 4, 0, 0}, {5, 2, 2, 1, 4, 0}, {1, 2, 3, 1, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 2}, {0, 0, 0, 2}, {1, 0, 0, 2}, {1, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(2431):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 1, 2, 3, 0}, {5, 1, 4, 2, 2, 0}, {3, 4, 0, 2, 2, 0}, {3, 0, 0, 2, 2, 0}, {0, 0, 4, 3, 2, 0}, {1, 4, 1, 3, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {2, 0, 0, 1}, {2, 0, 0, 0}, {2, 0, 0, 0}, {3, 0, 0, 0}, {3, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(4213):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 4, 0, 0}, {5, 2, 2, 1, 4, 0}, {1, 2, 3, 1, 1, 0}, {1, 3, 3, 4, 1, 0}, {0, 3, 2, 0, 4, 0}, {3, 2, 2, 0, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 2}, {1, 0, 0, 2}, {1, 0, 0, 3}, {1, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});
            break;
        }
        case(2413):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 6, 7}, {3, 0, 6, 7}, {5, 3, 6, 7}, {1, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 1, 2, 2, 2, 3}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 2, 0}, {0, 0, 4, 3, 2, 0}, {1, 4, 1, 3, 3, 0}, {1, 1, 1, 2, 3, 0}, {5, 1, 4, 2, 2, 0}, {3, 4, 0, 2, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}, {2, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 0}, {3, 0, 0, 0}, {3, 0, 0, 1}, {3, 0, 0, 1}, {2, 0, 0, 1}, {2, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}});

            break;
        }
        default :
        {
            std::cout<<"Template not found in the catalog"<<std::endl;
            break;
        }
        }
    }

    void
    hierarchy_tet4_4nc(moris::size_t const & aPermutationId)
    {
        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;
        switch(aPermutationId)
        {
        case(4520):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 3, 0, 4, 0}, {1, 2, 3, 4, 1, 0}, {2, 5, 2, 0, 1, 1}, {3, 2, 5, 4, 0, 1}, {0, 3, 2, 0, 0, 0}, {2, 2, 2, 0, 0, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3}, {1, 0, 0, 3}, {0, 1, 0, 2}, {0, 0, 1, 2}, {0, 0, 0, 3}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(2540):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3, 2, 0}, {3, 4, 0, 2, 2, 0}, {4, 5, 1, 0, 2, 2}, {1, 1, 5, 2, 3, 2}, {0, 0, 4, 3, 3, 0}, {4, 1, 1, 3, 0, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 0}, {2, 0, 0, 0}, {0, 2, 0, 1}, {3, 0, 2, 1}, {3, 0, 0, 0}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(4025):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 4, 0}, {3, 2, 2, 4, 0, 0}, {2, 0, 3, 0, 0, 0}, {1, 3, 0, 4, 1, 0}, {5, 2, 2, 1, 1, 0}, {2, 3, 3, 1, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 3}, {1, 0, 0, 3}, {1, 0, 0, 2}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(2045):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 2, 2, 0}, {1, 4, 1, 2, 3, 0}, {4, 0, 0, 0, 3, 3}, {3, 0, 0, 2, 2, 3}, {5, 1, 4, 2, 2, 0}, {4, 0, 0, 2, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {3, 0, 0, 1}, {0, 3, 0, 0}, {2, 0, 3, 0}, {2, 0, 0, 1}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(5301):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 3, 1, 5, 0}, {2, 0, 3, 5, 2, 0}, {0, 3, 0, 0, 2, 2}, {4, 0, 3, 5, 1, 2}, {1, 3, 0, 1, 1, 0}, {0, 0, 0, 1, 0, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {2, 0, 0, 3}, {0, 2, 0, 0}, {1, 0, 2, 0}, {1, 0, 0, 3}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(5103):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 5, 0}, {4, 0, 0, 5, 1, 0}, {0, 1, 3, 0, 1, 1}, {2, 3, 1, 5, 2, 1}, {3, 0, 0, 2, 2, 0}, {0, 3, 3, 2, 0, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 3}, {2, 0, 1, 3}, {2, 0, 0, 0}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(351 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 1, 3, 0, 0}, {4, 5, 1, 0, 0, 0}, {5, 3, 2, 0, 0, 0}, {2, 2, 3, 0, 3, 0}, {1, 1, 5, 3, 3, 0}, {5, 2, 2, 3, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 2}, {3, 0, 0, 2}, {3, 0, 0, 1}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(153 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 0, 0, 0}, {2, 5, 2, 0, 3, 0}, {5, 1, 1, 0, 3, 3}, {4, 1, 1, 0, 0, 3}, {3, 2, 5, 0, 0, 0}, {5, 1, 1, 0, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 2}, {3, 0, 0, 2}, {0, 3, 0, 1}, {0, 0, 3, 1}, {0, 0, 0, 2}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(3412):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 3, 2, 3, 0}, {0, 1, 3, 3, 0, 0}, {1, 4, 1, 0, 0, 0}, {5, 1, 4, 3, 2, 0}, {2, 3, 1, 2, 2, 0}, {1, 1, 1, 2, 0, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 1}, {2, 0, 0, 1}, {2, 0, 0, 3}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});

            break;
        }
        case(3214):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 3, 0}, {5, 1, 1, 3, 2, 0}, {1, 2, 3, 0, 2, 2}, {0, 3, 2, 3, 0, 2}, {4, 1, 1, 0, 0, 0}, {1, 3, 3, 0, 0, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1}, {2, 0, 0, 1}, {0, 2, 0, 3}, {0, 0, 2, 3}, {0, 0, 0, 1}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(1432):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 2, 3, 1, 0}, {5, 3, 2, 1, 1, 0}, {3, 4, 0, 0, 1, 1}, {0, 0, 4, 1, 3, 1}, {2, 2, 3, 3, 3, 0}, {3, 0, 0, 3, 0, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {1, 0, 0, 2}, {0, 1, 0, 0}, {3, 0, 1, 0}, {3, 0, 0, 2}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(1234):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{4, 0, 5, 7}, {0, 3, 5, 7}, {5, 3, 6, 7}, {1, 2, 6, 7}, {2, 4, 5, 7}, {2, 5, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 2, 1, 3}, {1, 1, 2, 1, 2, 3}, {1, 1, 2, 3, 2, 2}, {1, 2, 1, 1, 2, 2}, {1, 2, 1, 2, 2, 3}, {1, 2, 2, 2, 3, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 1, 1, 0}, {0, 3, 0, 1, 3, 0}, {3, 2, 2, 0, 3, 3}, {5, 2, 2, 1, 1, 3}, {4, 0, 3, 1, 1, 0}, {3, 2, 2, 1, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 2, 3, 2}, {2, 3, 2, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 0}, {3, 0, 0, 0}, {0, 3, 0, 2}, {1, 0, 3, 2}, {1, 0, 0, 0}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {2}, {std::numeric_limits<moris::size_t>::max()}, {1}, {1}});
            break;
        }
        case(5402):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 3, 3, 5, 2, 0}, {1, 3, 0, 1, 5, 0}, {0, 0, 4, 1, 0, 1}, {3, 4, 0, 2, 5, 1}, {2, 0, 3, 2, 2, 0}, {0, 0, 0, 0, 2, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 3}, {1, 0, 0, 3}, {0, 0, 1, 0}, {2, 1, 0, 0}, {2, 0, 0, 3}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(452 ):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{2, 2, 2, 0, 3, 0}, {3, 2, 5, 0, 0, 0}, {5, 1, 4, 0, 0, 0}, {1, 4, 1, 3, 0, 0}, {2, 5, 2, 3, 3, 0}, {5, 1, 1, 0, 3, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2}, {0, 0, 0, 2}, {0, 0, 0, 1}, {3, 0, 0, 1}, {3, 0, 0, 2}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(5204):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 0, 0, 5, 1, 0}, {3, 0, 0, 2, 5, 0}, {0, 3, 2, 2, 0, 2}, {1, 2, 3, 1, 5, 2}, {4, 0, 0, 1, 1, 0}, {0, 3, 3, 0, 1, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 0}, {2, 0, 0, 0}, {0, 0, 2, 3}, {1, 2, 0, 3}, {1, 0, 0, 0}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case( 254):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{4, 1, 1, 0, 0, 0}, {1, 1, 5, 3, 0, 0}, {5, 2, 2, 3, 0, 3}, {3, 2, 2, 0, 0, 3}, {4, 5, 1, 0, 0, 0}, {5, 2, 2, 0, 0, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1}, {3, 0, 0, 1}, {0, 0, 3, 2}, {0, 3, 0, 2}, {0, 0, 0, 1}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});
            break;
        }
        case(3510):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 3, 3, 3, 0, 0}, {2, 3, 1, 2, 3, 0}, {1, 1, 5, 2, 0, 2}, {4, 5, 1, 0, 3, 2}, {0, 1, 3, 0, 0, 0}, {1, 1, 1, 0, 0, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 3}, {2, 0, 0, 3}, {0, 0, 2, 1}, {0, 2, 0, 1}, {0, 0, 0, 3}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(1530):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 1, 3, 0}, {4, 0, 3, 1, 1, 0}, {3, 2, 5, 1, 0, 1}, {2, 5, 2, 3, 1, 1}, {0, 3, 0, 3, 3, 0}, {3, 2, 2, 0, 3, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 2}, {3, 1, 0, 2}, {3, 0, 0, 0}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(3015):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 1, 1, 3, 2, 0}, {4, 1, 1, 0, 3, 0}, {1, 3, 0, 0, 0, 0}, {2, 0, 3, 2, 3, 0}, {5, 1, 1, 2, 2, 0}, {1, 3, 3, 0, 2, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 3}, {2, 0, 0, 3}, {2, 0, 0, 1}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});
            break;
        }
        case(1035):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{5, 2, 2, 1, 1, 0}, {2, 2, 3, 3, 1, 0}, {3, 0, 0, 3, 0, 3}, {4, 0, 0, 1, 1, 3}, {5, 3, 2, 1, 1, 0}, {3, 0, 0, 0, 1, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 2}, {3, 0, 0, 2}, {0, 0, 3, 0}, {1, 3, 0, 0}, {1, 0, 0, 2}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(4321):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 3, 3, 4, 1, 0}, {0, 3, 2, 0, 4, 0}, {2, 2, 3, 0, 0, 0}, {5, 3, 2, 1, 4, 0}, {1, 2, 3, 1, 1, 0}, {2, 2, 2, 0, 1, 0}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1, 0, 0, 3}, {0, 0, 0, 3}, {0, 0, 0, 2}, {1, 0, 0, 2}, {1, 0, 0, 3}, {0, 0, 0, 2}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});
            break;
        }
        case(2341):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1, 1, 1, 2, 3, 0}, {5, 1, 4, 2, 2, 0}, {4, 0, 3, 2, 0, 2}, {0, 3, 0, 3, 2, 2}, {1, 4, 1, 3, 3, 0}, {4, 0, 0, 0, 3, 2}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 1}, {2, 0, 0, 1}, {0, 0, 2, 0}, {3, 2, 0, 0}, {3, 0, 0, 1}, {0, 0, 0, 0}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(4123):
        {

            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 2, 2, 4, 0, 0}, {5, 2, 2, 1, 4, 0}, {2, 3, 1, 1, 0, 1}, {0, 1, 3, 0, 4, 1}, {3, 2, 2, 0, 0, 0}, {2, 3, 3, 0, 0, 1}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0, 0, 0, 2}, {1, 0, 0, 2}, {0, 0, 1, 3}, {0, 1, 0, 3}, {0, 0, 0, 2}, {0, 0, 0, 3}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});

            break;
        }
        case(2143):
        {
            mNewElementToNode = moris::Matrix< moris::IndexMat >({{0, 4, 5, 7}, {3, 0, 5, 7}, {3, 5, 6, 7}, {2, 1, 6, 7}, {4, 2, 5, 7}, {5, 2, 6, 7}});
            mNumNewElem = 6;
            mNumElemToReplace = 1;
            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 2, 3}, {1, 2, 1, 2, 1, 3}, {1, 2, 1, 2, 3, 2}, {1, 1, 2, 2, 1, 2}, {1, 1, 2, 2, 2, 3}, {1, 2, 2, 3, 2, 2}});
            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3, 0, 0, 2, 2, 0}, {0, 0, 4, 3, 2, 0}, {4, 1, 1, 3, 0, 3}, {5, 1, 1, 2, 2, 3}, {3, 4, 0, 2, 2, 0}, {4, 1, 1, 0, 2, 3}});
            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2, 3, 3, 2}, {2, 3, 3, 2}, {3, 3, 2, 2}, {2, 2, 3, 2}, {2, 3, 3, 2}, {3, 3, 3, 2}});
            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{2, 0, 0, 0}, {3, 0, 0, 0}, {0, 0, 3, 1}, {2, 3, 0, 1}, {2, 0, 0, 0}, {0, 0, 0, 1}});
            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1}, {std::numeric_limits<moris::size_t>::max()}, {1}, {std::numeric_limits<moris::size_t>::max()}, {2}, {2}});
            break;
        }
        default :
        {
            std::cout<<"Template not found in the catalog"<<std::endl;
            break;
        }
        }
    }

    void
    bisected_tet(moris::size_t const & aPermutationId)
    {
        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;

        switch(aPermutationId)
        {
            case 0:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,4,2,3},
                                                                      {4,1,2,3}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 2, 1}, {1,1,2,2,1,1}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 3, 2, 3, 0, 5},{0,1,3,0,4,5}});
                mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2},{2,2,3,2}});
                mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 0, 2, 3},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{2}});
                break;
            }
            case 1:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,3},
                                                                      {0,4,2,3}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 1, 2, 1, 1, 2}, {2,1,1,1,2,1}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 3, 3, 4, 1},{3,1,2,3,1,5}});
                mNewParentFaceRanks       =
                        moris::Matrix< moris::DDSTMat >({{2, 2, 3, 2},{3,2,2,2}});
                mNewParentFaceOrdinals    =
                        moris::Matrix< moris::IndexMat>({{0, 1, 0, 3},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2},{0}});
                break;
            }
            case 2:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,3},
                                                                      {4,1,2,3}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 2, 1, 1, 1, 2}, {2,1,1,2,1,1}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 3, 2, 3, 4, 2},{3,1,2,2,4,5}});
                mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2},{3,2,2,2}});
                mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 0, 2, 3},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{0}});
                break;
            }
            case 3:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,4},
                                                                      {4,1,2,3}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 1, 1, 1, 2, 2}, {2,1,2,1,1,1}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 3, 0, 2},{0,1,2,3,4,5}});
                mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{2, 3, 2, 2},{2,2,2,3}});
                mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 0, 2, 3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{3}});
                break;
            }
            case 4:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,4,2,3},
                                                                      {0,1,2,4}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{2, 2, 1, 1, 1, 1}, {1,1,1,2,1,2}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 3, 4, 5},{0,1,2,0,4,1}});
                mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >({{2, 2, 2, 3},{2,2,3,2}});
                mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 0},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{3},{2}});
                break;
            }
            case 5:
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,3},
                                                                      {0,1,2,4}});
                mNumNewElem = 2;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >({{1, 2, 2, 1, 1, 1}, {1,1,1,2,2,1}});
                mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0, 1, 2, 3, 4, 5},{0,1,2,2,1,5}});
                mNewParentFaceRanks       =
                        moris::Matrix< moris::DDSTMat >({{2, 2, 2, 3},{3,2,2,2}});
                mNewParentFaceOrdinals    =
                        moris::Matrix< moris::IndexMat >({{0, 1, 2, 0},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{3},{0}});
                break;
            }
            default:
            {
                MORIS_ERROR(0,"Invalid edge ordinal found");
                break;
            }
        }
    }

    void
    hierarchy_tet4_2(moris::size_t const & aPermutationId)
    {
        mSpatialDimension = 3;
        mElementTopology = CellTopology::TET4;

        switch(aPermutationId)
        {
            case(14):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,5}, {0,4,2,5}, {0,5,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,1,2},{2,1,1,2,2,2},{2,2,1,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,0,4,1},{3,1,2,0,1,1},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2}, {3,2,3,2}, {2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,0,3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {0}, {std::numeric_limits<moris::size_t>::max()}});
                break;
            }

            case(41):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,4},{0,4,5,3},{0,5,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,1,2},{2,2,2,1,1,2},{2,1,1,1,2,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,0,4,1},{0,1,3,3,4,1},{3,1,2,3,1,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2}, {2,2,3,3}, {3,2,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,0,0},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {3}, {std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(15):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,5},{0,1,5,3},{0,4,2,5}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,2,2},{1,2,2,1,1,1},{2,1,1,2,2,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,2,1,1},{0,1,2,3,4,5},{3,1,2,2,1,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{3,2,3,2}, {2,2,2,3}, {3,2,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,2,0},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {0}});
                break;
            }
            case(51):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,3},{0,5,4,3},{0,5,2,4}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,1,1,2},{2,2,2,1,2,1},{2,1,1,2,2,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,3,4,1},{3,1,2,3,1,5},{3,1,2,2,1,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2}, {3,2,2,3}, {3,2,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,2,0},{0,1,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {3}, {0}});
                break;
            }
            case(45):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,4},{0,1,2,5},{0,4,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,2,2,1,2},{1,1,1,2,2,1},{2,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,0,4,1},{0,1,2,2,1,5},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,3}, {3,2,2,2}, {2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,0},{0,1,2,3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2}, {std::numeric_limits<moris::size_t>::max()}, {3}});
                break;
            }
            case(54):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,5},{0,5,2,4},{0,5,4,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,1,2,1,2},{2,2,1,2,2,1},{2,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,0,4,1},{0,1,2,2,1,5},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2}, {3,2,2,3}, {2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,2,0},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {0}, {3}});
                break;
            }
            case(53):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,5},{5,1,4,3},{1,2,5,4}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,1,1,2,2},{2,2,2,1,1,1},{1,2,2,2,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3,0,2},{0,1,2,3,4,5},{1,2,0,1,5,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2}, {2,2,2,3}, {2,2,3,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,1,2,0},{1,2,0,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()}, {3}, {2}});
                break;
            }
            case(35):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,5},{0,1,5,4},{4,1,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,1,2,2,1},{1,2,2,1,2,2},{2,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,2,1,5},{0,1,2,3,0,2},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{3,2,2,2},{2,3,2,3},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3},{0,0,2,0},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{1},{3}});
                break;
            }
            case(12):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{5,4,2,3},{1,4,5,3},{0,1,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{2,1,1,2,2,1},{1,2,2,1,2,2},{1,2,1,1,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3,1,2,2,1,5},{1,3,3,4,1,2},{0,3,2,3,4,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{3,2,2,2},{2,3,3,2},{2,3,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3,},{1,0,0,3},{0,0,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{0},{1},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(21):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{5,2,4,3},{0,5,4,3},{0,1,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,1,2},{2,2,1,1,2,2},{1,1,2,1,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{1,2,3,1,5,2},{3,3,2,3,1,2},{0,1,3,3,4,1}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2},{3,3,2,2},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{1,2,0,3},{0,0,2,3},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2},{1},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(02):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{5,1,2,3},{5,4,1,3},{0,4,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{2,1,1,2,1,1},{2,1,2,2,2,1},{1,2,1,1,2,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3,1,2,2,4,5},{3,0,3,2,0,4},{0,3,2,3,0,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{3,2,2,2},{3,2,3,2},{2,3,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3},{0,0,0,3},{0,0,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{0},{1}});
                break;
            }
            case(20):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{5,1,2,3},{4,5,2,3},{0,5,4,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,1,1},{2,2,1,2,2,1},{1,2,1,1,2,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,0,4,5},{3,3,2,2,0,5},{0,3,2,3,0,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2},{3,3,2,2},{2,3,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,0,2,3},{0,0,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{0},{1}});
                break;
            }
            case(01):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,5,2,3},{0,4,5,3},{4,1,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{2,1,1,1,2,1},{1,2,2,1,2,2},{1,1,2,2,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{3,1,2,3,1,5},{0,3,3,3,0,1},{0,1,3,0,4,1}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{3,2,2,2},{2,3,3,2},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3},{0,0,0,3},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{1},{2}});
                break;
            }
            case(10):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,5,2,3},{5,4,2,3},{5,1,4,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,1,2,1},{2,1,2,2,2,1},{1,1,2,2,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,0,5},{3,1,3,0,1,5},{0,1,3,0,4,1}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{3,2,3,2},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,1,0,3},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{0},{2}});
                break;
            }
            case(03):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,4,2,5},{4,1,2,5},{5,1,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,1,2,2},{1,1,2,2,2,2},{2,1,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,0,2},{0,1,3,0,0,2},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{2,3,3,2},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,0,0,3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{2},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(30):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,5,2,4},{4,5,2,3},{5,1,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,1,2,2},{2,2,2,1,2,1},{1,1,2,2,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,0,2},{0,3,2,3,0,5},{0,1,3,0,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{2,3,2,3},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,0,2,0},{0,1,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{2},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(04):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,4,2,5},{4,1,2,5},{0,5,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,2,2,2},{1,1,2,2,1,2},{2,2,1,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,0,0,1},{0,1,3,0,4,1},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,3,2},{2,2,3,2},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,0,3},{0,1,0,3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{2},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(40):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{5,1,2,4},{5,4,2,3},{0,5,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,2,2,1,2},{2,2,2,2,1,1},{1,2,1,1,2,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,3,0,4,1},{0,1,3,0,4,5},{0,3,2,3,0,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,3,2},{2,2,3,3},{2,3,2,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,0,3},{0,1,0,0},{0,0,2,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{2},{3},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(34):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,5,2,4},{0,1,2,5},{4,5,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{2,2,1,1,2,2},{1,1,1,2,1,2},{2,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,2,2,3,0,2},{0,1,2,0,4,1},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,3},{2,2,3,2},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,0},{0,1,0,3},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{std::numeric_limits<moris::size_t>::max()},{3}});
                break;
            }
            case(43):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,2,5},{5,1,2,4},{5,4,2,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,1,1,2,2},{2,1,2,2,1,2},{2,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3,0,2},{0,1,2,0,4,1},{0,1,2,3,4,5}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{2,2,3,3},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,1,0,0},{0,1,2,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{2},{3}});
                break;
            }
            case(25):
            {
                // done
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,3},{0,1,4,5},{1,2,4,5}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,2,1,1,1},{1,2,1,2,2,2},{1,1,2,2,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,3,4,5},{0,3,2,2,1,2},{1,2,3,1,5,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,2,2,3},{3,3,2,2},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,1,2,0},{0,0,2,3},{1,2,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{1},{2}});
                break;
            }
            case(52):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,3},{1,4,5,3},{1,2,5,4}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,1,1,2},{2,2,2,1,1,2},{1,1,2,2,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,4,2},{1,2,3,4,5,2},{1,2,3,1,5,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{2,2,3,3},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{1,2,0,0},{1,2,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{std::numeric_limits<moris::size_t>::max()},{3},{2}});
                break;
            }
            case(23):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,4,5},{4,1,2,5},{1,2,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,1,1,1,2,2},{2,1,1,2,2,2},{1,2,2,1,1,1}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,0,2},{3,1,2,2,0,2},{1,2,0,4,5,3}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{3,3,2,2},{2,2,2,3}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,0,2,3},{1,2,0,0}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{0},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            case(32):
            {
                mNewElementToNode = moris::Matrix< moris::IndexMat >({{0,1,5,4},{4,1,5,3},{1,2,5,3}});
                mNumNewElem = 3;
                mNumElemToReplace = 1;
                mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{1,2,1,1,2,2},{2,2,2,1,1,2},{1,1,2,1,1,2}});
                mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{0,3,2,3,0,2},{0,3,2,3,4,2},{1,2,3,4,5,2}});
                mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{2,3,2,2},{2,3,2,3},{2,2,3,2}});
                mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{0,0,2,3},{0,0,2,0},{1,2,0,3}});
                mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},{3},{std::numeric_limits<moris::size_t>::max()}});
                break;
            }
            //            mNewElementToNode = moris::Matrix< moris::IndexMat >({{},{},{}});
            //            mNumNewElem = 3;
            //            mNumElemToReplace = 1;
            //            mNewParentEdgeRanks = moris::Matrix< moris::DDSTMat >({{},{},{}});
            //            mNewParentEdgeOrdinals = moris::Matrix< moris::IndexMat >({{},{},{}});
            //            mNewParentFaceRanks = moris::Matrix< moris::DDSTMat >({{},{},{}});
            //            mNewParentFaceOrdinals = moris::Matrix< moris::IndexMat >({{},{},{}});
            //            mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{},{},{}});
            default :
            {
                std::cout<<"Template not found in the catalog"<<std::endl;
                break;
            }
        }
    }

};

}
#endif /* SRC_XTK_CL_XTK_CHILD_MESH_MODIFICATION_TEMPLATE_HPP_ */

