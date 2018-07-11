/*
 * cl_MeshXTK.cpp
 *
 *  Created on: Nov 11, 2016
 *      Author: doble
 */

#include "cl_Mesh_XTK.hpp"

//-----------------------------------------------------
moris::MeshXTK::MeshXTK(moris::uint D)
{
    mD = D;        // Set mesh dimensionality for indexing purposes
    mNumElems = 1; // Initialize the number of elements to one
}
//-----------------------------------------------------
moris::MeshXTK::~MeshXTK()
{

}
//-----------------------------------------------------
void
moris::MeshXTK::generate_templated_mesh(enum TemplateType    aTemplate)
{
    switch (aTemplate)
    {
        case(TemplateType::REGULAR_SUBDIVISION_HEX8):
                {

            MORIS_ASSERT(mNodesInd.n_cols()==15,"For a Hex8 regular subdivision template, there must be 15 node ids stored in MeshXTK Object");

            // Load template data from .inc files

#include"templates/reg_sub/reg_sub_0d1d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_0d2d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_0d3d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_1d0d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_1d2d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_1d3d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_2d0d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_2d1d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_2d3d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_3d0d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_3d1d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_3d2d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_3d3d.inc" // snippets XTK
            mNumElems = 24;
            mIndices.resize(16,moris::Mat<moris::uint>(1,180,0)); // Resize indices member variable
            mOffsets.resize(16,moris::Mat<moris::uint>(1,61,0));  // Resize offsets member variable

            //Insert indices
            //        mIndices.insert(0,NULL);          // 0d->0d indices data (Not implemented)
            mIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            mIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            mIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            mIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            //        mIndices.insert(5,NULL);          // 1d->1d indices data (Not implemented)
            mIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            mIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            mIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            mIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            //        mIndices.insert(10, NULL);        // 2d->2d indices data (Not Implemented)
            mIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            mIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            mIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            mIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            mIndices.insert(15, tIndices3d3d);// 3d->3d indices data

            // Insert offsets data
            //        mOffsets.insert(0, NULL);         // 0d->0d offsets data (Not Implemented)
            mOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            mOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            mOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            mOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            //        mOffsets.insert(5, NULL);         // 1d->1d offsets data (Not Implemented)
            mOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            mOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            mOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            mOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            //        mOffsets.insert(10, NULL);         // 2d->2d offsets data (Not Implemented)
            mOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            mOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            mOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            mOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            mOffsets.insert(15, tOffsets3d3d);// 3d->3d offsets data

            break;
                }
        case(TemplateType::QUAD_4):
            {
#include"templates/quad_4/quad_4_0d1d.inc" // snippets XTK
#include"templates/quad_4/quad_4_0d2d.inc" // snippets XTK
#include"templates/quad_4/quad_4_1d0d.inc" // snippets XTK
#include"templates/quad_4/quad_4_1d2d.inc" // snippets XTK
#include"templates/quad_4/quad_4_2d0d.inc" // snippets XTK
#include"templates/quad_4/quad_4_2d1d.inc" // snippets XTK

            MORIS_ASSERT(mNodesInd.n_cols()==4," For a quad 4 element there needs to be 4 nodes in XTK Mesh");
            MORIS_ASSERT(mD == 2, "Model dimension for a Quad 4 needs to be 2");
            // Resize indices and offsets
            mIndices.resize(16,moris::Mat<moris::uint>(1,8,0));
            mOffsets.resize(16,moris::Mat<moris::uint>(1,5,0));

            //Insert indices data
            mIndices.insert(0,tIndices0d1d); // 0d->1d indices data
            mIndices.insert(1,tIndices0d2d); // 0d->2d indices data
            mIndices.insert(2,tIndices1d0d); // 1d->0d indices data
            //        mIndices.insert(3,NULL);         // 1d->1d indices data (Not implemented)
            mIndices.insert(4, tIndices1d2d);// 1d->2d indices data
            mIndices.insert(5, tIndices2d0d);// 2d->0d indices data
            mIndices.insert(6, tIndices2d1d);// 2d->1d indices data
            //        mIndices.insert(7, NULL);        // 2d->2d indices data (Not Implemented)



            // Insert offsets data
            mOffsets.insert(0, tOffsets0d1d);// 0d->1d offsets data
            mOffsets.insert(1, tOffsets0d2d);// 0d->2d offsets data
            mOffsets.insert(2, tOffsets1d0d);// 1d->1d offsets data
            //        mOffsets.insert(3, NULL);        // 1d->1d offsets data (Not Implemented)
            mOffsets.insert(4, tOffsets1d2d);// 1d->2d offsets data
            mOffsets.insert(5, tOffsets2d0d);// 2d->0d offsets data
            mOffsets.insert(6, tOffsets2d1d);// 2d->1d offsets data
            //        mOffsets.insert(7, NULL);        // 2d->2d offsets data (Not Implemented)
            break;

            }

        case(TemplateType::TET_4):
                {
            // Used for template generate
#include "templates/tet_4/tet_4.inc" // snippets XTK
            MORIS_ASSERT(mNodesInd.n_cols()==4," For a tet 4 element there needs to be 4 nodes in XTK Mesh");
            MORIS_ASSERT(mD ==3,"Tet 4 is a 3 dimensional element -model dimension should be 3");
            mNumElems = 1;

            mIndices.resize(16,moris::Mat<moris::uint>(1,15,0));
            mOffsets.resize(16,moris::Mat<moris::uint>(1,15,0));
            moris::Mat<moris::uint> temp(0,0);

            mIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            mIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            mIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            mIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            mIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            mIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            mIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            mIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            mIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            mIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            mIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            mIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            mIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            mIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            mIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            mIndices.insert(15, temp);// 3d->3d indices data

            // Insert offsets data
            mOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            mOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            mOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            mOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            mOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            mOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            mOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            mOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            mOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            mOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            mOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            mOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            mOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            mOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            mOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            mOffsets.insert(15, temp);// 3d->3d offsets data
            break;

            break;
                }

        case(TemplateType:: HIERARCHY_TET4_3N):
            {
#include"templates/node_hier_tet4/node_hier_3_node_inds.inc"
#include"templates/node_hier_tet4/node_hier_3_node_offs.inc"

            mIndices.resize(8,moris::Mat<moris::uint>(1,8,0));
            mOffsets.resize(8,moris::Mat<moris::uint>(1,5,0));
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            mIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            mIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            mIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            mIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            mIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            mIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            mIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            mIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            mIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            mIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            mIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            mIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            mIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            mIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            mIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            mIndices.insert(15, temp);// 3d->3d indices data

            // Insert offsets data
            mOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            mOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            mOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            mOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            mOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            mOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            mOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            mOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            mOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            mOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            mOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            mOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            mOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            mOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            mOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            mOffsets.insert(15, temp);// 3d->3d offsets data
            break;
            }

        case(TemplateType:: HIERARCHY_TET4_4Na):
            {
#include"templates/node_hier_tet4/node_hier_4_node_a_inds.inc"
#include"templates/node_hier_tet4/node_hier_4_node_a_offs.inc"

            mIndices.resize(8,moris::Mat<moris::uint>(1,8,0));
            mOffsets.resize(8,moris::Mat<moris::uint>(1,5,0));
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            mIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            mIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            mIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            mIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            mIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            mIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            mIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            mIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            mIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            mIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            mIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            mIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            mIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            mIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            mIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            mIndices.insert(15, temp);// 3d->3d indices data

            // Insert offsets data
            mOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            mOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            mOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            mOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            mOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            mOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            mOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            mOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            mOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            mOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            mOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            mOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            mOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            mOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            mOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            mOffsets.insert(15, temp);// 3d->3d offsets data
            break;
            }
        default:
            mNumElems = 0;
            MORIS_ASSERT(mNumElems!=0,"Template not implemented");
            break;
    }
}

//-----------------------------------------------------

void
moris::MeshXTK::set_entity_ancestry(enum TemplateType                        aTemplate,
                                    moris::Cell<moris::Mat<moris::uint>>   & aParentEntities)
{
    MORIS_ASSERT(aParentEntities.size()==mD+1,"There needs to be a ancestry space allocated for each type of entity even if they are not used (node,edge,face,element)");
    switch (aTemplate)
    {
        case(TemplateType::REGULAR_SUBDIVISION_HEX8):
                {
#include"templates/reg_sub/reg_sub_anc_1d.inc" // snippets XTK
#include"templates/reg_sub/reg_sub_anc_2d.inc" // snippets XTK

            mAncestryInds.resize(4,moris::Mat<moris::uint>(1,1,UINT_MAX));
            mAncestryRank.resize(4,moris::Mat<moris::uint>(1,1,UINT_MAX));
            mAncestryInds.insert(1, tAnc1DInds);
            mAncestryRank.insert(1, tAnc1DRank);
            mAncestryInds.insert(2, tAnc2D);

            break;
                }
        default:
            MORIS_LOG_ERROR<<"Ancestry for template not implemented";
            break;
    }
}

void
moris::MeshXTK::modify_templated_mesh(enum TemplateType    aTemplate)
{
    moris::uint tModCount=0;
    switch(aTemplate)
    {
        case(TemplateType::HIERARCHY_TET4):
                {
            MORIS_ASSERT(mAuxConnectivity.n_rows()==get_num_entities(EntityRank::ELEMENT),"Auxiliary connectivity should have a row for each element, Has it been setup?");
            for(moris::uint e=0; e<mAuxConnectivity.n_rows();e++)
            {
                MORIS_ASSERT(mAuxConnectivity(e,0)==3||mAuxConnectivity(e,0)==4 ||mAuxConnectivity(e,0)==0,"Invalid element intersection pattern for nodal hierarchy");
                if(mAuxConnectivity(e,0)==3)
                {
                    enum TemplateType tTemplate= TemplateType::HIERARCHY_TET4;
                    moris::Mat<moris::uint> tSortedNodes = sort_nodes(tTemplate, mAuxConnectivity.row(e));

                    tTemplate= TemplateType::HIERARCHY_TET4_3N;
                    // case number embedded in the sorted nodes(at end).
                    MORIS_ASSERT(tSortedNodes(0,7) == UINT_MAX,"Invalid template type returned from sorting nodes");

                    // Size out embedded case number
                    tSortedNodes.resize(1,7);

                    // Access inds and offsets (do not want the member connectivities here)
                    // These connectivities will be inserted into the member connectivities
                    moris::Cell<moris::Mat<moris::uint>> tInds = access_template_indices(tTemplate);
                    moris::Cell<moris::Mat<moris::uint>> tOffs = access_template_offsets(tTemplate);


                    // Replace and insert elements only for right now
                    // Replace parent element with first child element
                    moris::uint tNumChildrenElem = tOffs(12).n_cols()-1;
                    MORIS_ASSERT(tNumChildrenElem == 4,"There should only be 4 children elements for the 3 node case");
                    moris::Mat<moris::uint> tNodeElemConn(tNumChildrenElem,4,UINT_MAX);
                    for(moris::uint i = 0; i<tNumChildrenElem; i++)
                    {
                        tNodeElemConn.row(i) = access_connectivity(tInds(12),tOffs(12),i).row(0);
                    }
                    replace_parent_entity(EntityRank::ELEMENT,e,tNodeElemConn,tSortedNodes);

                    tModCount++;
                }
                else if(mAuxConnectivity(e,0) == 4)
                {
                    //MORIS_LOG_INFO<<mAuxConnectivity.row(e);
                    enum TemplateType tTemplate= TemplateType::HIERARCHY_TET4;
                    moris::Mat<moris::uint> tSortedNodes = sort_nodes(tTemplate, mAuxConnectivity.row(e));

                    if(tSortedNodes(0,8) == 0 )      tTemplate = TemplateType::HIERARCHY_TET4_4Na;
                    else if(tSortedNodes(0,8) == 1 ) tTemplate = TemplateType::HIERARCHY_TET4_4Nb;
                    else if(tSortedNodes(0,8) == 2 ) tTemplate = TemplateType::HIERARCHY_TET4_4Nc;
                    else MORIS_LOG_ERROR<<"Invalid case number output by sort_nodes";
                    // Size out embedded case number
                    tSortedNodes.resize(1,8);

                    // Access inds and offsets (do not want the member connectivities here)
                    // These connectivities will be inserted into the member connectivities
                    moris::Cell<moris::Mat<moris::uint>> tInds = access_template_indices(tTemplate);
                    moris::Cell<moris::Mat<moris::uint>> tOffs = access_template_offsets(tTemplate);

                    moris::uint tNumChildrenElem = tOffs(12).n_cols()-1;
                    MORIS_ASSERT(tNumChildrenElem == 6,"There should only be 6 children elements for the 4 node case");
                    moris::Mat<moris::uint> tNodeElemConn(tNumChildrenElem,4,UINT_MAX);
                    for(moris::uint i = 0; i<tNumChildrenElem; i++)
                    {
                        tNodeElemConn.row(i) = access_connectivity(tInds(12),tOffs(12),i).row(0);
                    }
                    //MORIS_LOG_INFO<<"Node to Element Connectivity\n"<<tNodeElemConn;
                    replace_parent_entity(EntityRank::ELEMENT,e,tNodeElemConn,tSortedNodes);
                    tModCount++;
                }
                else if(mAuxConnectivity(e,0) == 0) continue;
                else MORIS_LOG_ERROR<<"Invalid connectivity for nodal hierarchy template, (should be 3 or 4 nodes)";
            }
            break;
                }
        default:
        {
            MORIS_LOG_ERROR<<"Invalid template specified";
            break;
        }
    }
    //MORIS_LOG_INFO<<"Elements modified: "<< tModCount;
}

void
moris::MeshXTK::init_aux_connectivity(enum EntityRank aD,
                                      enum EntityRank aDPrime1,
                                      enum EntityRank aDPrime2)
{
    moris::uint tNumD         = get_num_entities(aD);
    mAuxConnectivityNum       = get_num_entities_connected_to_entity(aD,aDPrime1,0);
    mAuxConnectivity = moris::Mat<moris::uint>(tNumD,mAuxConnectivityNum*2+1,0); // Needs to be zero for use column
    mAuxD = aD;
    mAuxDPrime1 = aDPrime1;
    mAuxDPrime2 = aDPrime2;
}

void
moris::MeshXTK::add_entity_to_aux_connectivity(moris::uint aDPrime1Ind,
                                               moris::uint aDPrime2Ind,
                                               moris::uint aFlag)
{
    if(aFlag==0)
    {
        aDPrime1Ind = aDPrime1Ind+get_num_entities(mAuxDPrime1);
    }

    // Get entities of dimension mAuxDPrime2 connected to entities of mAuxD (i.e. elements connected to given edge aDPrime2Ind)
    moris::Mat<moris::uint> tDInd = get_entities_connected_to_entity(mAuxDPrime2,mAuxD,aDPrime2Ind,0);

    for(moris::uint i = 0; i<tDInd.n_cols(); i++)
    {
        MORIS_ASSERT(mAuxConnectivity(tDInd(0,i),0)<mAuxConnectivityNum, "Entity corresponding to provided aDInd has exceeded allocated space");
        MORIS_ASSERT(tDInd(0,i)<mAuxConnectivity.n_rows(),"aDInd is outside of bounds. Has auxiliary connectivity been initialized?");

        mAuxConnectivity(tDInd(0,i),mAuxConnectivity(tDInd(0,i),0)+1)                     = aDPrime1Ind;
        mAuxConnectivity(tDInd(0,i),mAuxConnectivity(tDInd(0,i),0)+mAuxConnectivityNum+1) = aDPrime2Ind;
        mAuxConnectivity(tDInd(0,i),0)++;
    }

}

moris::Mat<moris::uint>
moris::MeshXTK::get_aux_connectivity()
{
    return mAuxConnectivity;
}


void
moris::MeshXTK::set_aux_connectivity(moris::Mat<moris::uint>  aAuxConn)
{
    mAuxConnectivity = aAuxConn;
}

moris::Mat< moris::uint >
moris::MeshXTK::sort_nodes(enum TemplateType       aTemplate,
                           moris::Mat<moris::uint> aAuxConnectivity)
{
    ////MORIS_LOG_INFO<<aAuxConnectivity;
    ////MORIS_LOG_INFO<<mNodesInd;
    //Locate highest node in auxilliary connectivity
    switch(aTemplate)
    {
        case(TemplateType::HIERARCHY_TET4):
                    {
            bool tdbgflg = 0;
            if(aAuxConnectivity(0,0) == 3)
            {
                moris::uint tIntersectionCase = UINT_MAX;
                moris::Mat<moris::uint> tHigh = {{aAuxConnectivity(0,1), aAuxConnectivity(0,5)}};
                moris::Mat<moris::uint> tMid  = {{aAuxConnectivity(0,1), aAuxConnectivity(0,5)}};
                moris::Mat<moris::uint> tLow  = {{aAuxConnectivity(0,1), aAuxConnectivity(0,5)}};

                for(moris::uint i = 0; i<2; i++)
                {
                    if(mNodesId(0,aAuxConnectivity(0,i+2))>mNodesId(0,tHigh(0,0)))
                    {
                        tMid.row(0) = tHigh.row(0);
                        tHigh(0,0) = aAuxConnectivity(0,i+2);
                        tHigh(0,1) = aAuxConnectivity(0,i+6);
                    }
                    else if(mNodesId(0,aAuxConnectivity(0,i+2))<mNodesId(0,tLow(0,0)))
                    {
                        tMid.row(0) = tLow.row(0);
                        tLow(0,0) = aAuxConnectivity(0,i+2);
                        tLow(0,1) = aAuxConnectivity(0,i+6);
                    }
                    else
                    {
                        tMid(0,0) = aAuxConnectivity(0,i+2);
                        tMid(0,1) = aAuxConnectivity(0,i+6);
                    }
                }

                if(tdbgflg)
                {
                    //MORIS_LOG_INFO<<"High: " <<tHigh;
                    //MORIS_LOG_INFO<<"Mid:  " <<tMid;
                    //MORIS_LOG_INFO<<"Low:  " <<tLow;
                }
                moris::Mat<moris::uint>  tNodes14 = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tHigh(0,1),0);
                moris::Mat<moris::uint>  tNodes13 = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tMid(0,1),0);
                moris::Mat<moris::uint>  tNodes12 = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tLow(0,1),0);

                if(tdbgflg)
                {
                    //MORIS_LOG_INFO<<"Nodes 1,4:  "<<tNodes14;
                    //MORIS_LOG_INFO<<"Nodes 1,3:  "<<tNodes13;
                    //MORIS_LOG_INFO<<"Nodes 1,2:  "<<tNodes12;
                }
                // Find the shared node in all the above lists

                // Decide which node is which (using intersections)
                moris::uint tN1 = UINT_MAX;
                moris::uint tN2 = UINT_MAX;
                moris::uint tN3 = UINT_MAX;
                moris::uint tN4 = UINT_MAX;

                // Find nodes 1,4,3
                if(tNodes14(0,0) == tNodes13(0,0) )
                {
                    tN1 = tNodes14(0,0);
                    tN4 = tNodes14(0,1);
                    tN3 = tNodes13(0,1);
                }
                else if(tNodes14(0,0) == tNodes13(0,1))
                {
                    tN1 = tNodes14(0,0);
                    tN4 = tNodes14(0,1);
                    tN3 = tNodes13(0,0);
                }
                else if(tNodes14(0,1) == tNodes13(0,0) )
                {
                    tN1 = tNodes14(0,1);
                    tN4 = tNodes14(0,0);
                    tN3 = tNodes13(0,1);
                }
                else if(tNodes14(0,1) == tNodes13(0,1))
                {
                    tN1 = tNodes14(0,1);
                    tN4 = tNodes14(0,0);
                    tN3 = tNodes13(0,0);
                }
                else MORIS_LOG_ERROR<<"Duplicate node not found, invalid edge intersection configuration";


                // Find node 2
                if(tN1 == tNodes12(0,0))
                {
                    tN2 = tNodes12(0,1);
                }
                else if(tN1 == tNodes12(0,1))
                {
                    tN2 = tNodes12(0,0);
                }
                else MORIS_LOG_ERROR<<"Node 2 not found, invalid edge intersection configuration";

                moris::Mat<moris::uint> tSortedNodes = {{tN1,tN2,tN3,tN4,tLow(0,0),tMid(0,0),tHigh(0,0),tIntersectionCase}};
                //if(tdbgflg)MORIS_LOG_INFO<<"Sorted Nodes: "<<tSortedNodes;

#ifdef MORIS_USE_TESTS
                moris::Mat<moris::uint> tDuplNodes = moris::Debug::duplicate_col_check(tSortedNodes);
                MORIS_ASSERT(tDuplNodes.n_rows()==0,"Duplicate node appeared in sorting.");
#endif
                return tSortedNodes;
            }

            else if(aAuxConnectivity(0,0) == 4)
            {
                // Sort Auxiliary nodes from highest to lowest using bubble sort (but based on first column then swap row

                moris::Mat<moris::uint> tSortedNodes(1,9);
                moris::Mat<moris::uint> tNodes     = {{aAuxConnectivity(0,1), aAuxConnectivity(0,5)},
                                                      {aAuxConnectivity(0,2), aAuxConnectivity(0,6)},
                                                      {aAuxConnectivity(0,3), aAuxConnectivity(0,7)},
                                                      {aAuxConnectivity(0,4), aAuxConnectivity(0,8)}};

                bool swapped = true;
                moris::uint j = 0;
                moris::uint n = 4; // 4 numbers to sort
                moris::Mat<moris::uint> tmp(1,2);
                while (swapped)
                {
                    swapped = false;
                    j++;
                    for (moris::uint i = 0; i < n - j; i++)
                    {
                        if (mNodesId(0,tNodes(i,0)) > mNodesId(0,tNodes(i + 1,0)))
                        {
                            tmp.row(0)        = tNodes.row(i);
                            tNodes.row(i)     = tNodes.row(i + 1);
                            tNodes.row(i + 1) = tmp.row(0);
                            swapped = true;
                        }
                    }
                }
                if(tdbgflg)
                {
                    //MORIS_LOG_INFO<<"Sorted Nodes: \n"<<tNodes;
                }
                // Determine the relationship between high and low

                moris::Mat<moris::uint>  tNodesL  = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tNodes(0,1),0);
                moris::Mat<moris::uint>  tNodesML = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tNodes(1,1),0);
                moris::Mat<moris::uint>  tNodesMH = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tNodes(2,1),0);
                moris::Mat<moris::uint>  tNodesH  = get_entities_connected_to_entity(EntityRank::EDGE,EntityRank::NODE,tNodes(3,1),0);


                if(tdbgflg)
                {
                //MORIS_LOG_INFO<<"Nodes connected to low Id: \n"<<tNodesL;
                //MORIS_LOG_INFO<<"Nodes connected to mid-low Id:\n"<<tNodesML;
                //MORIS_LOG_INFO<<"Nodes connected to mid-high Id:\n"<<tNodesMH;
                //MORIS_LOG_INFO<<"Nodes connected to high Id:\n"<<tNodesH;
                }
                moris::uint tHLOppFlag  = 1;
                moris::uint tHMHOppFlag = 1;
                moris::uint tHMLOppFlag = 1;

                // If L and H share a node then it is not case a
                for(moris::uint i = 0; i<2; i++)
                {
                    if(tNodesL(0,i) == tNodesH(0,0))
                    {
                        tHLOppFlag=0;
                    }

                    else if(tNodesL(0,i) == tNodesH(0,1))
                    {
                        tHLOppFlag=0;
                    }
                }

                // If MH and H share a node then its not case b
                for(moris::uint i = 0; i<2; i++)
                {
                    if(tNodesMH(0,i) == tNodesH(0,0))
                    {
                        tHMHOppFlag=0;
                    }

                    else if(tNodesMH(0,i) == tNodesH(0,1))
                    {
                        tHMHOppFlag=0;
                    }
                }

                // If ML and H share a node then its not case c
                for(moris::uint i = 0; i<2; i++)
                {
                    if(tNodesML(0,i) == tNodesH(0,0))
                    {
                        tHMLOppFlag=0;
                    }

                    else if(tNodesML(0,i) == tNodesH(0,1))
                    {
                        tHMLOppFlag=0;
                    }
                }
                if(tdbgflg)
                {
                    MORIS_LOG_INFO<<"High low Opposite Flag: "<<tHLOppFlag;
                    MORIS_LOG_INFO<<"High MH Opposite Flag: "<<tHMHOppFlag;
                    MORIS_LOG_INFO<<"High ML Opposite Flag: "<<tHMLOppFlag;
                }
                if(tHLOppFlag)
                {
                    tSortedNodes(0,8)=0; // Indicating that you need to use template node_hier_4_node_a.inc
                    tSortedNodes(0,4)= tNodes(0,0);
                    tSortedNodes(0,5)= tNodes(1,0);
                    tSortedNodes(0,6)= tNodes(2,0);
                    tSortedNodes(0,7)= tNodes(3,0);

                    // Get node shared by MH and H
                    moris::uint j=1;
                    moris::uint tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesMH(0,i) == tNodesH(0,0))
                        {
                            tSortedNodes(0,0) = tNodesMH(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,1);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesMH(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesMH(0,i) == tNodesH(0,1))
                        {
                            tSortedNodes(0,0) = tNodesMH(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,0);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesMH(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find nodes 1,2,3 unsuccessful");

                    // Midlow and Low to find last node
                    j=1;
                    tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesML(0,i) == tNodesL(0,0))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0,i) == tNodesL(0,1))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find node 4 unsuccessful");

                    if(tdbgflg)
                    {
                        //MORIS_LOG_INFO<<"Sorted Nodes: "<<tSortedNodes;
                    }
                }
                else if(tHMHOppFlag)
                {
                    tSortedNodes(0,8)=1;// Indicating that you need to use template node_hier_4_node_b.inc
                    tSortedNodes(0,4)= tNodes(0,0); // Low
                    tSortedNodes(0,5)= tNodes(1,0); // Mid low
                    tSortedNodes(0,6)= tNodes(2,0); // Mid high
                    tSortedNodes(0,7)= tNodes(3,0); // High

                    // Get node shared by H and L
                    moris::uint j=1;
                    moris::uint tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesL(0,i) == tNodesH(0,0))
                        {
                            tSortedNodes(0,0) = tNodesL(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,1);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesL(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesL(0,i) == tNodesH(0,1))
                        {
                            tSortedNodes(0,0) = tNodesL(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,0);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesL(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find nodes 1,2,3 unsuccessful");

                    // Midlow and MidHigh to find last node
                    j=1;
                    tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesML(0,i) == tNodesMH(0,0))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0,i) == tNodesMH(0,1))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find node 4 unsuccessful");
                    if(tdbgflg)
                    {
                        //MORIS_LOG_INFO<<"Sorted Nodes: "<<tSortedNodes;
                    }
                }
                else if (tHMLOppFlag)
                {
                    tSortedNodes(0,8)=2;// Indicating that you need to use template node_hier_4_node_b.inc
                    tSortedNodes(0,4)= tNodes(0,0); // Low
                    tSortedNodes(0,5)= tNodes(1,0); // Mid low
                    tSortedNodes(0,6)= tNodes(2,0); // Mid high
                    tSortedNodes(0,7)= tNodes(3,0); // High

                    // Get node shared by H and L
                    moris::uint j=1;
                    moris::uint tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesL(0,i) == tNodesH(0,0))
                        {
                            tSortedNodes(0,0) = tNodesL(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,1);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesL(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else if(tNodesL(0,i) == tNodesH(0,1))
                        {
                            tSortedNodes(0,0) = tNodesL(0,i); // Shared Node
                            tSortedNodes(0,1) = tNodesH(0,0);  // High nodes independent node
                            tSortedNodes(0,2) = tNodesL(0,j); // Mid highs ind node
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find nodes 1,2,3 unsuccessful");

                    // Midlow and MidHigh to find last node
                    j=1;
                    tSuccess = 0;
                    for(moris::uint i = 0; i<2; i++)
                    {
                        if(tNodesML(0,i) == tNodesMH(0,0))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else if(tNodesML(0,i) == tNodesMH(0,1))
                        {
                            tSortedNodes(0,3) = tNodesML(0,i);
                            tSuccess = 1;
                        }

                        else j = i;
                    }
                    MORIS_ASSERT(tSuccess=1,"Sorting to find node 4 unsuccessful");
                }
                else MORIS_LOG_ERROR << "Sorting Failed (invalid flagging). Did a node appear twice?";
                return tSortedNodes;
            }

            else
            {
                MORIS_LOG_ERROR << "SORTING NOT COMPLETED! Check to see if this function is called for a non-intersected element";
                moris::Mat< moris::uint > dummy;
                return dummy;
            }

            break;
                    }

        default:
        {
            MORIS_LOG_ERROR << "Sorting for specified template type not implemented";

            moris::Mat< moris::uint > dummy;
            return dummy;

            break;
        }

}

}

//-----------------------------------------------------

void
moris::MeshXTK::set_node_index(moris::Mat<moris::uint>    &aNodes)
{
    mNodesInd = aNodes;
}
//-----------------------------------------------------
void
moris::MeshXTK::set_node_ids(moris::Mat<moris::uint>    & aNodeIds)
{

    MORIS_ASSERT(aNodeIds.n_cols()>=aNodeIds.n_rows(),"Node ids should be in a column vector");
    mNodesId = aNodeIds;
}
//-----------------------------------------------------
void
moris::MeshXTK::set_pending_node_index_pointers(moris::Cell<moris::uint*>  aNodeIndPtr)
{
    mPtrPendingNodeIndex = aNodeIndPtr;
}
//-----------------------------------------------------
void
moris::MeshXTK::get_pending_node_inds()
{
    moris::Mat<moris::uint> tempNodeMat(1,mPtrPendingNodeIndex.size(),0);
    for(moris::uint i = 0; i<mPtrPendingNodeIndex.size(); i++)
    {
        tempNodeMat(0,i) = *mPtrPendingNodeIndex(i);
    }

    add_node_ind(tempNodeMat);

    mPtrPendingNodeIndex.resize(0,NULL);
}
//-----------------------------------------------------
void
moris::MeshXTK::add_node_ind(moris::Mat<moris::uint>    &aNodeInd)
{

    moris::uint tNumNewNodes = aNodeInd.n_cols();       // Get num of nodes being added
    moris::uint tNumExistingNodes = mNodesInd.n_cols();  // Get num of existing nodes
    moris::uint tNumEndNodes = tNumNewNodes+tNumExistingNodes; // Get the total number of nodes after the function

    //Resize member node variable
    mNodesInd.resize(1,tNumEndNodes);

    //Loop over new nodes add to existing node list
    for(moris::uint i = tNumExistingNodes; i < tNumEndNodes; i++)
    {
        mNodesInd(0,i) = aNodeInd(0,i-tNumExistingNodes);
    }
}
//-----------------------------------------------------
moris::uint
moris::MeshXTK::get_node_ind(moris::uint  aIndex)
{
    MORIS_ASSERT(aIndex<=mNodesInd.n_cols(),"Requested node is outside of bounds.");
    moris::uint tNID = mNodesInd(0,aIndex);
    return tNID;
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::get_all_node_inds()
{
    return mNodesInd;
}
//-----------------------------------------------------
moris::uint
moris::MeshXTK::get_num_entities(enum EntityRank aEntity)
{
    moris::uint tNum = 0;
    if(aEntity==EntityRank::NODE)
    {
        tNum = mNodesInd.n_cols();
    }
    else
    {
        moris::uint i = index_map((moris::uint)aEntity,0);
        tNum = mOffsets(i).n_cols()-1;
    }

    return tNum;
}
//-----------------------------------------------------
void
moris::MeshXTK::set_parent_element_index(moris::uint    aEID)
{
    mParentEID = aEID;
}
//-----------------------------------------------------
moris::uint
moris::MeshXTK::get_parent_element_index()
{
    return mParentEID;
}
//-----------------------------------------------------
moris::uint
moris::MeshXTK::get_num_entities_connected_to_entity(enum EntityRank aEntity,
                                                     enum EntityRank aEntityPrime,
                                                     moris::uint     aEntityIndex)
{
    moris::uint i = index_map((moris::uint)aEntity,(moris::uint)aEntityPrime);
    moris::uint tOff        = mOffsets(i)(0,aEntityIndex);   // first offset
    moris::uint tOffEnd     = mOffsets(i)(0,aEntityIndex+1); // second offset
    moris::uint tOffDiff    = tOffEnd-tOff;                 // Number of edges connected to element
    return tOffDiff;
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::get_entities_connected_to_entity(enum EntityRank aEntity,
                                                 enum EntityRank aEntityPrime,
                                                 moris::uint     aEntityIndex,
                                                 moris::uint     aType)
{
    return access_connectivity(aEntity,aEntityPrime,aEntityIndex,aType);
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::get_parent_entity(enum EntityRank aEntityRank,
                                  moris::uint     aEntityIndex)
{
    return {{mAncestryRank(moris::uint(aEntityRank))(0,aEntityIndex),mAncestryInds(moris::uint(aEntityRank))(0,aEntityIndex)}};
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::access_connectivity(enum EntityRank             aEntity,
                                    enum EntityRank             aEntityPrime,
                                    moris::uint                 aEntityIndex,
                                    moris::uint                 aType) const
{

    moris::uint ind         = index_map((moris::uint)aEntity,(moris::uint)aEntityPrime);
    moris::uint tOff        = mOffsets(ind)(0,aEntityIndex);          // first offset
    moris::uint tOffEnd     = mOffsets(ind)(0,aEntityIndex+1);        // second offset
    moris::uint tOffDiff    = tOffEnd-tOff;                           // Expected size of aMat
    MORIS_ASSERT(tOff<=tOffEnd,"Beginning offset is larger than ending offset");

    moris::Mat<moris::uint> tMat(1,tOffDiff,UINT_MAX);

    for (moris::uint m = 0; m<tOffDiff; m++)
    {
        //moris::cout<<mIndices(aIndex).data()(0,tOff);
        if(aType==1 && aEntityPrime==EntityRank::NODE)
        {
            tMat(0,m) = mNodesInd(0,mIndices(ind)(0,tOff));
        }
        else
        {
            tMat(0,m) = mIndices(ind)(0,tOff);
        }
        tOff++;
    }
    return tMat;
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::access_connectivity(moris::Mat<moris::uint>  & aInds,
                                    moris::Mat<moris::uint>  & aOffs,
                                    moris::uint                aEntityIndex) const
{
    moris::uint tOff        = aOffs(0,aEntityIndex);          // first offset
    moris::uint tOffEnd     = aOffs(0,aEntityIndex+1);        // second offset
    moris::uint tOffDiff    = tOffEnd-tOff;                   // Expected size of aMat
    MORIS_ASSERT(tOff<=tOffEnd,"Beginning offset is larger than ending offset");

    moris::Mat<moris::uint> tMat(1,tOffDiff,UINT_MAX);

    for (moris::uint m = 0; m<tOffDiff; m++)
    {

        tMat(0,m) = aInds(0,tOff);
        tOff++;
    }
    return tMat;
}
//-----------------------------------------------------
void
moris::MeshXTK::set_connectivity(moris::uint                d,
                                 moris::uint                dprime,
                                 moris::Mat<moris::uint>    aIndices,
                                 moris::Mat<moris::uint>    aOffsets)
{
    // Convert d and dprime to appropriate index for mOffsets,mIndices.
    moris::uint i = index_map(d,dprime);

    //Resize
    moris::uint SizeInd = aIndices.n_cols();
    moris::uint SizeOff = aOffsets.n_cols();

    mIndices(i).resize(1,SizeInd);
    mOffsets(i).resize(1,SizeOff);

    mIndices(i) = aIndices;
    mOffsets(i) = aOffsets;

    //    std::cout<<"\n"<<d<<" to "<<dprime<<" Offsets\n"<<aOffsets;
    //    std::cout<<"\n"<<d<<" to "<<dprime<<" Indices\n"<<aIndices;
}
//-----------------------------------------------------
moris::Mat<moris::uint>
moris::MeshXTK::get_full_connectivity(enum EntityRank    d,
                                      enum EntityRank    dprime,
                                      moris::uint        aType)
{
    moris::uint ind       = index_map((moris::uint)d,(moris::uint)dprime);
    moris::uint tNumD     = get_num_entities(d);
    moris::uint temp      = get_num_entities_connected_to_entity(d,dprime,0);
    moris::uint tOffStart = 0;
    moris::uint tOffEnd   = 0;
    moris::uint l         = 0;
    moris::Mat<moris::uint> tConnectivity(tNumD,temp,UINT_MAX);
    for(moris::uint j = 0; j<tNumD;j++)
    {
        tOffEnd   = mOffsets(ind)(0,j+1);
        for(moris::uint k = tOffStart; k<tOffEnd;k++)
        {
            if(aType == 0)
            {
                tConnectivity(j,l) = mIndices(ind)(0,k);
            }
            else if(aType == 1)
            {
                if(dprime == EntityRank::NODE)
                    tConnectivity(j,l) = mNodesInd(0,mIndices(ind)(0,k));
                else
                    tConnectivity(j,l) = mIndices(ind)(0,k);
            }

            else if(aType ==2)
            {
                if(dprime == EntityRank::NODE)
                    tConnectivity(j,l) = mNodesId(0,mIndices(ind)(0,k));
                else
                    tConnectivity(j,l) = mIndices(ind)(0,k);
            }
            else
            {
                moris::uint breaker = 0;
                MORIS_ASSERT(breaker!=0,"Incorrect type specified, 0 for local xtk index and 1 for processor local index");
            }

            l++;
        }
        l         = 0;
        tOffStart = tOffEnd;
    }
    return tConnectivity;
}
//-----------------------------------------------------
moris::uint
moris::MeshXTK::index_map(moris::uint    d,
                          moris::uint    dprime) const
{
    moris::uint i = (mD+1) * d + dprime;
    MORIS_ASSERT(i<=mIndices.size(),"The requested index is out of bounds");
    return i;
}
//-----------------------------------------------------
void
moris::MeshXTK::replace_parent_entity(enum EntityRank  aParentRank,
                                      moris::uint      aParentInd,
                                      moris::Mat<moris::uint> tChildrenConn,
                                      moris::Mat<moris::uint> tSortedNodes)
{
    moris::uint tNumChildren = tChildrenConn.n_rows();
    //MORIS_LOG_INFO<< "Replacing parent "<<aParentInd<< " with rank " <<(moris::uint)aParentRank<<" with "<<tNumChildren<<" children";

    // Replace parent with child number one
    // Later check if the parent entity has the same number of nodes as children entity
    moris::uint ind = index_map((moris::uint)aParentRank,0);
    moris::uint tOff        = mOffsets(ind)(0,aParentInd);          // first offset
    moris::uint tOffEnd     = mOffsets(ind)(0,aParentInd+1);        // second offset
    moris::uint tOffDiff    = tOffEnd-tOff;                           // Expected size of aMat

    MORIS_ASSERT(tOffDiff==tChildrenConn.n_cols(),"Size of child does not match size of parent (see comment later check)");

    //    //MORIS_LOG_INFO<<"Before first child element insertion";
    //    //MORIS_LOG_INFO<<"Indices: \n"<<mIndices(ind);
    //    //MORIS_LOG_INFO<<"Offsets: \n"<<mOffsets(ind);

    moris::uint j=0;
    for(moris::uint i = tOff; i<tOffEnd; i++)
    {
        mIndices(ind)(0,i) = tSortedNodes(0,tChildrenConn(0,j));
        j++;
    }

    // Add remaining children to end of offsets and indices
    // Resize indice and offset
    moris::uint tIndInd = mIndices(ind).n_cols(); // first available ind
    moris::uint tOffInd = mOffsets(ind).n_cols(); // first available offset

    //    //MORIS_LOG_INFO<<"After first child element insertion";
    //    //MORIS_LOG_INFO<<"Indices: \n"<<mIndices(ind);
    //    //MORIS_LOG_INFO<<"Offsets: \n"<<mOffsets(ind);
    mIndices(ind).resize(1,tIndInd+(tNumChildren-1)*tOffDiff);
    mOffsets(ind).resize(1,tOffInd+tNumChildren-1);

    for(moris::uint cr = 1; cr<tNumChildren;cr++)// child rows
    {
        for(moris::uint cc = 0; cc<tOffDiff; cc++)// child columns
        {
            mIndices(ind)(0,tIndInd) = tSortedNodes(0,tChildrenConn(cr,cc));
            tIndInd++;
        }
        mOffsets(ind)(0,tOffInd) = tIndInd;
        tOffInd++;
    }

    //    //MORIS_LOG_INFO<<"After all child element insertions";
    //    //MORIS_LOG_INFO<<"Indices: \n"<<mIndices(ind);
    //    //MORIS_LOG_INFO<<"Offsets: \n"<<mOffsets(ind);

}
//-----------------------------------------------------
moris::Cell<moris::Mat<moris::uint>>
moris::MeshXTK::access_template_indices(enum TemplateType  aTemplate)
{
    switch(aTemplate)
    {
        case(TemplateType:: HIERARCHY_TET4_3N):
                {
#include"templates/node_hier_tet4/node_hier_3_node_inds.inc"

            moris::Cell<moris::Mat<moris::uint>> tIndices(16);
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            //            tIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            tIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            tIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            tIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            tIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            //            tIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            tIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            tIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            tIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            tIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            //            tIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            tIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            tIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            tIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            tIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            //            tIndices.insert(15, temp);// 3d->3d indices data

            return tIndices;
            break;
                }
        case(TemplateType::HIERARCHY_TET4_4Na):
               {
#include"templates/node_hier_tet4/node_hier_4_node_a_inds.inc"
            moris::Cell<moris::Mat<moris::uint>> tIndices(16);
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            //            tIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            tIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            tIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            tIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            tIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            //            tIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            tIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            tIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            tIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            tIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            //            tIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            tIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            tIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            tIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            tIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            //            tIndices.insert(15, temp);// 3d->3d indices data

            return tIndices;
            break;
               }
        case(TemplateType::HIERARCHY_TET4_4Nb):
                            {
#include"templates/node_hier_tet4/node_hier_4_node_b_inds.inc"

            moris::Cell<moris::Mat<moris::uint>> tIndices(16);
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            //            tIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            tIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            tIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            tIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            tIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            //            tIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            tIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            tIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            tIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            tIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            //            tIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            tIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            tIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            tIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            tIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            //            tIndices.insert(15, temp);// 3d->3d indices data

            return tIndices;
            break;
                            }
        case(TemplateType::HIERARCHY_TET4_4Nc):
                                        {
#include"templates/node_hier_tet4/node_hier_4_node_c_inds.inc"

            moris::Cell<moris::Mat<moris::uint>> tIndices(16);
            moris::Mat<moris::uint> temp(0,0);

            //Insert indices
            //            tIndices.insert(0,temp);          // 0d->0d indices data (Not implemented)
            tIndices.insert(1,tIndices0d1d);  // 0d->1d indices data
            tIndices.insert(2,tIndices0d2d);  // 0d->2d indices data
            tIndices.insert(3,tIndices0d3d);  // 0d->3d indices data
            tIndices.insert(4,tIndices1d0d);  // 1d->0d indices data
            //            tIndices.insert(5,temp);          // 1d->1d indices data (Not implemented)
            tIndices.insert(6, tIndices1d2d); // 1d->2d indices data
            tIndices.insert(7, tIndices1d3d); // 1d->3d indices data
            tIndices.insert(8, tIndices2d0d); // 2d->0d indices data
            tIndices.insert(9, tIndices2d1d); // 2d->1d indices data
            //            tIndices.insert(10, temp);        // 2d->2d indices data (Not Implemented)
            tIndices.insert(11, tIndices2d3d);// 2d->3d indices data
            tIndices.insert(12, tIndices3d0d);// 3d->0d indices data
            tIndices.insert(13, tIndices3d1d);// 3d->1d indices data
            tIndices.insert(14, tIndices3d2d);// 3d->2d indices data
            //            tIndices.insert(15, temp);// 3d->3d indices data

            return tIndices;
            break;
                                        }

        default:
        {
            MORIS_LOG_ERROR<<"Invalid Template. Did access the wrong one by accident? Or has it not been implemented?";
            return 0;
            break;
        }
    }
}
//-----------------------------------------------------
moris::Cell<moris::Mat<moris::uint>>
moris::MeshXTK::access_template_offsets(enum TemplateType aTemplate)
{
    switch(aTemplate)
    {
        case(TemplateType:: HIERARCHY_TET4_3N):
            {
#include"templates/node_hier_tet4/node_hier_3_node_offs.inc"
            moris::Cell<moris::Mat<moris::uint>> tOffsets(16);
            moris::Mat<moris::uint> temp(0,0);


            //Insert offsets data
            //            tOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            tOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            tOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            tOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            tOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            //            tOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            tOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            tOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            tOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            tOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            //            tOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            tOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            tOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            tOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            tOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            //            tOffsets.insert(15, temp);// 3d->3d offsets data
            return tOffsets;
            break;
            }
        case(TemplateType::HIERARCHY_TET4_4Na):
                    {
#include"templates/node_hier_tet4/node_hier_4_node_a_offs.inc"

            moris::Cell<moris::Mat<moris::uint>> tOffsets(16);
            moris::Mat<moris::uint> temp(0,0);


            //Insert offsets data
            //            tOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            tOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            tOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            tOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            tOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            //            tOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            tOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            tOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            tOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            tOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            //            tOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            tOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            tOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            tOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            tOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            //            tOffsets.insert(15, temp);// 3d->3d offsets data
            return tOffsets;
            break;
                    }
        case(TemplateType::HIERARCHY_TET4_4Nb):
                    {
#include"templates/node_hier_tet4/node_hier_4_node_b_offs.inc"

            moris::Cell<moris::Mat<moris::uint>> tOffsets(16);
            moris::Mat<moris::uint> temp(0,0);


            //Insert offsets data
            //            tOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            tOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            tOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            tOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            tOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            //            tOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            tOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            tOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            tOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            tOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            //            tOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            tOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            tOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            tOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            tOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            //            tOffsets.insert(15, temp);// 3d->3d offsets data
            return tOffsets;
            break;
                    }
        case(TemplateType::HIERARCHY_TET4_4Nc):
                {
#include"templates/node_hier_tet4/node_hier_4_node_c_offs.inc"

            moris::Cell<moris::Mat<moris::uint>> tOffsets(16);
            moris::Mat<moris::uint> temp(0,0);


            //Insert offsets data
            //            tOffsets.insert(0, temp);         // 0d->0d offsets data (Not Implemented)
            tOffsets.insert(1, tOffsets0d1d); // 0d->1d offsets data
            tOffsets.insert(2, tOffsets0d2d); // 0d->2d offsets data
            tOffsets.insert(3, tOffsets0d3d); // 0d->3d offsets data
            tOffsets.insert(4, tOffsets1d0d); // 1d->1d offsets data
            //            tOffsets.insert(5, temp);         // 1d->1d offsets data (Not Implemented)
            tOffsets.insert(6, tOffsets1d2d); // 1d->2d offsets data
            tOffsets.insert(7, tOffsets1d3d); // 1d->3d offsets data
            tOffsets.insert(8, tOffsets2d0d); // 2d->0d offsets data
            tOffsets.insert(9, tOffsets2d1d); // 2d->1d offsets data
            //            tOffsets.insert(10, temp);         // 2d->2d offsets data (Not Implemented)
            tOffsets.insert(11, tOffsets2d3d);// 2d->3d offsets data
            tOffsets.insert(12, tOffsets3d0d);// 3d->0d offsets data
            tOffsets.insert(13, tOffsets3d1d);// 3d->1d offsets data
            tOffsets.insert(14, tOffsets3d2d);// 3d->2d offsets data
            //            tOffsets.insert(15, temp);// 3d->3d offsets data
            return tOffsets;
            break;
                }
        default:
        {
            MORIS_LOG_ERROR<<"Invalid Template";
            return 0;
            break;
        }
    }
}
//-----------------------------------------------------


