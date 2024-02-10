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
#include "cl_MTK_Enums.hpp"
namespace moris::xtk
{
    class Mesh_Modification_Template
    {
      public:
        Mesh_Modification_Template()
        {
        }

        // Constructor when the template does not have node inheritance information
        Mesh_Modification_Template( moris::moris_index aParentElemInd,
                moris::size_t                          aElemToReplace,
                Matrix< IndexMat > const              &aNodeInds,
                Matrix< IndexMat > const              &aParentEdgeInds,
                Matrix< moris::DDSTMat > const        &aParentEdgeRanks,
                Matrix< IndexMat > const              &aParentFaceInds,
                Matrix< moris::DDSTMat > const        &aParentFaceRanks,
                enum TemplateType                      aTemplateType,
                moris::size_t                          aPermutationId = 0 )
                : mParentNodeInds( 0, 0 )
                , mParentNodeRanks( 0, 0 )
                , mSpatialDimension( 0 )
                , mElementTopology( moris::mtk::CellTopology::UNDEFINED )
        {
            mElemIndToReplace = aElemToReplace;
            mParentElemInd    = aParentElemInd;
            mNodeInds         = aNodeInds.copy();
            mParentEdgeInds   = aParentEdgeInds.copy();
            mParentEdgeRanks  = aParentEdgeRanks.copy();
            mParentFaceInds   = aParentFaceInds.copy();
            mParentFaceRanks  = aParentFaceRanks.copy();
            template_catalog( aPermutationId, aTemplateType );
        }

        // Constructor for when the template has node inheritance information
        Mesh_Modification_Template( moris::moris_index aParentElemInd,
                moris::size_t                          aElemToReplace,
                Matrix< IndexMat > const              &aNodeInds,
                Matrix< IndexMat > const              &aParentNodeInds,
                Matrix< moris::DDSTMat > const        &aParentNodeRanks,
                Matrix< IndexMat > const              &aParentEdgeInds,
                Matrix< moris::DDSTMat > const        &aParentEdgeRanks,
                Matrix< IndexMat > const              &aParentFaceInds,
                Matrix< moris::DDSTMat > const        &aParentFaceRanks,
                enum TemplateType                      aTemplateType,
                moris::size_t                          aPermutationId = 0 )
        {
            MORIS_ASSERT( aTemplateType == TemplateType::REGULAR_SUBDIVISION_HEX8 || aTemplateType == TemplateType::REGULAR_SUBDIVISION_QUAD4 || aTemplateType == TemplateType::TET_4 || aTemplateType == TemplateType::TRI_3,
                    "Entered node inheritance constructor for template without node inheritance. Currently, only the regular subdivision hex8 and tet templates have node inheritance" );

            mTemplateType     = aTemplateType;
            mElemIndToReplace = aElemToReplace;
            mParentElemInd    = aParentElemInd;
            mNodeInds         = aNodeInds.copy();
            mParentNodeInds   = aParentNodeInds.copy();
            mParentNodeRanks  = aParentNodeRanks.copy();
            mParentEdgeInds   = aParentEdgeInds.copy();
            mParentEdgeRanks  = aParentEdgeRanks.copy();
            mParentFaceInds   = aParentFaceInds.copy();
            mParentFaceRanks  = aParentFaceRanks.copy();
            template_catalog( aPermutationId, aTemplateType );
        }

        // Constructor for when the 2D template has node inheritance information
        // Separate constructor could not be implemented because non-inheritance constructor has same parameters
        /*Mesh_Modification_Template( moris::moris_index                       aParentElemInd,
                                    moris::size_t                            aElemToReplace,
                                    Matrix< IndexMat > const & aNodeInds,
                                    Matrix< IndexMat > const & aParentNodeInds,
                                    Matrix< moris::DDSTMat >  const & aParentNodeRanks,
                                    Matrix< IndexMat > const & aParentEdgeInds,
                                    Matrix< moris::DDSTMat >  const & aParentEdgeRanks,
                                    Matrix< IndexMat > const & aParentFaceInds,
                                    Matrix< moris::DDSTMat >  const & aParentFaceRanks,
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
        moris::size_t mNumNewElem;
        moris::size_t mNumElemToReplace;
        moris::size_t mElemIndToReplace;

        // Parent Entity's parent information
        // This is the information relative to the parent this template is created from
        moris::moris_index       mParentElemInd;
        Matrix< IndexMat >       mParentNodeInds;
        Matrix< moris::DDSTMat > mParentNodeRanks;
        Matrix< IndexMat >       mParentEdgeInds;
        Matrix< moris::DDSTMat > mParentEdgeRanks;
        Matrix< IndexMat >       mParentFaceInds;
        Matrix< moris::DDSTMat > mParentFaceRanks;

        // spatial dimension of template
        moris::uint mSpatialDimension;

        // topology of children elements
        moris::mtk::CellTopology mElementTopology;

        // Node indices in the template
        Matrix< IndexMat > mNodeInds;

        // Element to Node Connectivity
        Matrix< IndexMat > mNewElementToNode;

        // Parent's of an entity ordered by ordinal relative to
        // (these are the inheritance of the new elements created by this template)
        Matrix< moris::DDSTMat > mNewParentEdgeRanks;
        Matrix< IndexMat >       mNewParentEdgeOrdinals;
        Matrix< moris::DDSTMat > mNewParentFaceRanks;
        Matrix< IndexMat >       mNewParentFaceOrdinals;
        Matrix< moris::DDSTMat > mNewElementInterfaceSides;

        // Node inheritance information
        bool                     mHasNodeInheritance = false;
        Matrix< moris::DDSTMat > mNewNodeParentRanks;
        Matrix< IndexMat >       mNewNodeParentOrdinals;

        // Reindex flag (meaning this template has been reindexed by the Child Mesh)
        bool mIsReindexed = false;

      private:
        // Note everything below this line is just template data and selecting the correct template
        void
        template_catalog( moris::size_t const &aPermutationId,
                enum TemplateType              aTemplateType )
        {
            switch ( aTemplateType )
            {
                case ( TemplateType::REGULAR_SUBDIVISION_HEX8 ):
                {
                    hex_8_reg_sub_template();
                    break;
                }
                case ( TemplateType::TET_4 ):
                {
                    tet4_template();
                    break;
                }
                case ( TemplateType::REGULAR_SUBDIVISION_QUAD4 ):
                {
                    quad_4_reg_sub_template();
                    break;
                }
                case ( TemplateType::TRI_3 ):
                {
                    tri3_template();
                    break;
                }
                case ( TemplateType::CONFORMAL_TRI3 ):
                {
                    conformal_tri_template( aPermutationId );
                    break;
                }

                case ( TemplateType::HIERARCHY_TET4_3N ):
                {
                    hierarchy_tet4_3N( aPermutationId );
                    break;
                }
                case ( TemplateType::HIERARCHY_TET4_4Na ):
                {
                    hierarchy_tet4_4na( aPermutationId );
                    break;
                }
                case ( TemplateType::HIERARCHY_TET4_4Nb ):
                {
                    hierarchy_tet4_4nb( aPermutationId );
                    break;
                }
                case ( TemplateType::HIERARCHY_TET4_4Nc ):
                {
                    hierarchy_tet4_4nc( aPermutationId );
                    break;
                }
                case ( TemplateType::BISECTED_TET4 ):
                {
                    bisected_tet( aPermutationId );
                    break;
                }
                case ( TemplateType::HIERARCHY_TET4_2 ):
                {
                    hierarchy_tet4_2( aPermutationId );
                    break;
                }
                default:
                    std::cout << "Template not found in the catalog" << std::endl;
                    break;
            }
        }

        void
        hex_8_reg_sub_template()
        {
            MORIS_ASSERT( mNodeInds.n_cols() == 15, "For a Hex8 regular subdivision template, there must be 15 node inds." );

            mNewElementToNode = Matrix< IndexMat >( { { 0, 8, 1, 14 },
                    { 1, 8, 5, 14 },
                    { 4, 5, 8, 14 },
                    { 0, 4, 8, 14 },
                    { 1, 9, 2, 14 },
                    { 2, 9, 6, 14 },
                    { 5, 6, 9, 14 },
                    { 1, 5, 9, 14 },
                    { 2, 10, 3, 14 },
                    { 2, 6, 10, 14 },
                    { 6, 7, 10, 14 },
                    { 3, 10, 7, 14 },
                    { 0, 3, 11, 14 },
                    { 3, 7, 11, 14 },
                    { 4, 11, 7, 14 },
                    { 0, 11, 4, 14 },
                    { 0, 1, 12, 14 },
                    { 1, 2, 12, 14 },
                    { 2, 3, 12, 14 },
                    { 0, 12, 3, 14 },
                    { 4, 13, 5, 14 },
                    { 5, 13, 6, 14 },
                    { 6, 13, 7, 14 },
                    { 4, 7, 13, 14 } } );

            mSpatialDimension         = 3;
            mElementTopology          = moris::mtk::CellTopology::TET4;
            mNumNewElem               = 24;
            mNumElemToReplace         = 1;
            mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 2, 2, 1, 3, 3, 3 }, { 1, 2, 2, 3, 3, 3 } } );
            mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 0, 0, 0 }, { 0, 0, 9, 0, 0, 0 }, { 4, 0, 0, 0, 0, 0 }, { 8, 0, 0, 0, 0, 0 }, { 1, 1, 1, 0, 0, 0 }, { 1, 1, 10, 0, 0, 0 }, { 5, 1, 1, 0, 0, 0 }, { 9, 1, 1, 0, 0, 0 }, { 2, 2, 2, 0, 0, 0 }, { 10, 2, 2, 0, 0, 0 }, { 6, 2, 2, 0, 0, 0 }, { 2, 2, 11, 0, 0, 0 }, { 3, 3, 3, 0, 0, 0 }, { 11, 3, 3, 0, 0, 0 }, { 3, 3, 7, 0, 0, 0 }, { 3, 3, 8, 0, 0, 0 }, { 0, 4, 4, 0, 0, 0 }, { 1, 4, 4, 0, 0, 0 }, { 2, 4, 4, 0, 0, 0 }, { 4, 4, 3, 0, 0, 0 }, { 5, 5, 4, 0, 0, 0 }, { 5, 5, 5, 0, 0, 0 }, { 5, 5, 6, 0, 0, 0 }, { 7, 5, 5, 0, 0, 0 } } );
            mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 }, { 3, 3, 3, 2 } } );
            mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 4 }, { 0, 0, 0, 4 }, { 0, 0, 0, 4 }, { 0, 0, 0, 4 }, { 0, 0, 0, 5 }, { 0, 0, 0, 5 }, { 0, 0, 0, 5 }, { 0, 0, 0, 5 } } );
            mNewElementInterfaceSides = Matrix< moris::DDSTMat >( 24, 1, std::numeric_limits< moris::size_t >::max() );

            // Mark this template as having node inheritance
            mHasNodeInheritance    = true;
            mNewNodeParentRanks    = Matrix< moris::DDSTMat >( { { 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 3 } } );
            mNewNodeParentOrdinals = Matrix< IndexMat >( { { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 0 } } );
        }

        void
        tet4_template()
        {
            mSpatialDimension         = 3;
            mElementTopology          = moris::mtk::CellTopology::TET4;
            mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 3 } } );
            mNumNewElem               = 1;
            mNumElemToReplace         = 0;
            mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 1, 1, 1 } } );
            mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 4, 5 } } );
            mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 2, 2 } } );
            mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3 } } );
            mNewElementInterfaceSides = Matrix< moris::DDSTMat >( 1, 1, std::numeric_limits< moris::size_t >::max() );
        }

        void
        quad_4_reg_sub_template()
        {
            MORIS_ASSERT( mNodeInds.n_cols() == 5, "For a Quad4 regular subdivision template, there must be 5 node inds." );
            mSpatialDimension = 2;
            mElementTopology  = moris::mtk::CellTopology::TRI3;
            mNewElementToNode = { { 0, 1, 4 },
                { 1, 2, 4 },
                { 2, 3, 4 },
                { 3, 0, 4 } };

            mNumNewElem               = 4;
            mNumElemToReplace         = 1;
            mNewParentEdgeRanks       = { { 1, 3, 3 },
                      { 1, 3, 3 },
                      { 1, 3, 3 },
                      { 1, 3, 3 } };
            mNewParentEdgeOrdinals    = { { 0, 0, 0 },
                   { 1, 0, 0 },
                   { 2, 0, 0 },
                   { 3, 0, 0 } };
            mNewElementInterfaceSides = Matrix< moris::DDSTMat >( 4, 1, std::numeric_limits< moris::size_t >::max() );

            // Mark template as having node inheritance
            mHasNodeInheritance    = true;
            mNewNodeParentRanks    = { { 0, 0, 0, 0, 3 } };
            mNewNodeParentOrdinals = { { 0, 1, 2, 3, 0 } };

            mElementTopology = moris::mtk::CellTopology::TRI3;
        }

        void
        tri3_template()
        {
            mSpatialDimension         = 2;
            mElementTopology          = moris::mtk::CellTopology::TRI3;
            mNewElementToNode         = { { 0, 1, 2 } };
            mNumNewElem               = 1;
            mNumElemToReplace         = 0;
            mNewParentEdgeRanks       = { { 1, 1, 1 } };
            mNewParentEdgeOrdinals    = { { 0, 1, 2 } };
            mNewElementInterfaceSides = Matrix< moris::DDSTMat >( 1, 1, std::numeric_limits< moris::size_t >::max() );
        }

        void
        conformal_tri_template( moris::size_t const &aPermutation )
        {

            mSpatialDimension  = 2;
            moris::size_t tMax = std::numeric_limits< moris::size_t >::max();

            switch ( aPermutation )
            {
                case 1:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 3, 2 }, { 3, 1, 4 }, { 3, 4, 2 } };
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 3, 1 }, { 1, 1, 3 }, { 3, 1, 3 } };
                    mNewParentEdgeOrdinals    = { { 0, 0, 2 }, { 0, 1, 0 }, { 0, 1, 0 } };
                    mNewElementInterfaceSides = { { tMax }, { 2 }, { 0 } };
                    break;
                }
                case 2:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 3, 4 }, { 3, 2, 4 }, { 3, 1, 2 } };
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 3, 1 }, { 3, 1, 3 }, { 1, 1, 3 } };
                    mNewParentEdgeOrdinals    = { { 0, 0, 2 }, { 0, 2, 0 }, { 0, 1, 0 } };
                    mNewElementInterfaceSides = { { 1 }, { 2 }, { tMax } };
                    break;
                }
                case 3:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 1, 3 }, { 0, 3, 4 }, { 3, 2, 4 } };
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 1, 3 }, { 3, 3, 1 }, { 1, 1, 3 } };
                    mNewParentEdgeOrdinals    = { { 0, 1, 0 }, { 0, 0, 2 }, { 1, 2, 0 } };
                    mNewElementInterfaceSides = { { tMax }, { 1 }, { 2 } };
                    break;
                }
                case 10:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 3, 2 }, { 3, 1, 2 } };
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 3, 1 }, { 1, 1, 3 } };
                    mNewParentEdgeOrdinals    = { { 0, 0, 2 }, { 0, 1, 0 } };
                    mNewElementInterfaceSides = { { 1 }, { 2 } };
                    break;
                }
                case 11:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 1, 3 }, { 0, 3, 2 } };
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 1, 3 }, { 3, 1, 1 } };
                    mNewParentEdgeOrdinals    = { { 0, 1, 0 }, { 0, 1, 2 } };
                    mNewElementInterfaceSides = { { 2 }, { 0 } };
                    break;
                }
                case 12:
                {
                    mElementTopology          = moris::mtk::CellTopology::TRI3;
                    mNewElementToNode         = { { 0, 1, 3 }, { 1, 2, 3 } };
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = { { 1, 3, 1 }, { 1, 1, 3 } };
                    mNewParentEdgeOrdinals    = { { 0, 0, 2 }, { 1, 2, 0 } };
                    mNewElementInterfaceSides = { { 1 }, { 2 } };
                    break;
                }
                default:
                    MORIS_ERROR( 0, "Invalid conformal tri 3 permutation" );
                    break;
            }
        }

        void
        hierarchy_tet4_3N( moris::size_t const &aPermutation )
        {

            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;
            switch ( aPermutation )
            {
                case ( 320 ):
                {
                    // Permutation 320
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 0, 3, 3, 0, 0, 2 }, { 1, 2, 3, 0, 2, 2 }, { 4, 5, 1, 0, 3, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 0, 3 }, { 0, 2, 0, 3 }, { 0, 2, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }

                case ( 32 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 3, 0, 3, 0 }, { 2, 2, 2, 3, 3, 0 }, { 5, 3, 2, 3, 0, 0 }, { 1, 4, 5, 3, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 3, 0, 0, 2 }, { 0, 0, 0, 2 }, { 3, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );

                    break;
                }

                case ( 203 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 2, 3 }, { 3, 0, 0, 2, 2, 3 }, { 4, 0, 0, 2, 3, 3 }, { 5, 1, 4, 2, 2, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 3, 0 }, { 2, 0, 0, 0 }, { 0, 3, 0, 0 }, { 2, 3, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 251 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 5, 2, 3, 2 }, { 1, 1, 1, 3, 3, 2 }, { 4, 5, 1, 3, 2, 2 }, { 0, 3, 4, 3, 2, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 2, 1 }, { 3, 0, 0, 1 }, { 0, 2, 0, 1 }, { 3, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }

                case ( 512 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 1, 5, 2, 1 }, { 2, 3, 3, 2, 2, 1 }, { 0, 1, 3, 2, 1, 1 }, { 3, 4, 0, 2, 5, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 1, 3 }, { 2, 0, 0, 3 }, { 0, 1, 0, 3 }, { 2, 1, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 125 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 1, 3 }, { 5, 2, 2, 1, 1, 3 }, { 3, 2, 2, 1, 3, 3 }, { 4, 0, 3, 1, 1, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 3, 2 }, { 1, 0, 0, 2 }, { 0, 3, 0, 2 }, { 1, 3, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 140 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 4, 1, 3, 1 }, { 0, 0, 0, 3, 3, 1 }, { 3, 4, 0, 3, 1, 1 }, { 2, 5, 3, 3, 1, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 1, 0 }, { 3, 0, 0, 0 }, { 0, 1, 0, 0 }, { 3, 1, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );

                    break;
                }

                case ( 401 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 0, 4, 1, 0 }, { 1, 3, 3, 1, 1, 0 }, { 2, 0, 3, 1, 0, 0 }, { 5, 3, 2, 1, 4, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 1, 0, 0, 3 }, { 0, 0, 0, 3 }, { 1, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );

                    break;
                }
                case ( 14 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 0, 3 }, { 4, 1, 1, 0, 0, 3 }, { 5, 1, 1, 0, 3, 3 }, { 3, 2, 5, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 3, 1 }, { 0, 0, 0, 1 }, { 0, 3, 0, 1 }, { 0, 3, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 453 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 5, 4, 0, 1 }, { 3, 2, 2, 0, 0, 1 }, { 2, 5, 2, 0, 1, 1 }, { 0, 1, 2, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 1, 2 }, { 0, 0, 0, 2 }, { 0, 1, 0, 2 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 534 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 3, 5, 1, 2 }, { 4, 0, 0, 1, 1, 2 }, { 0, 3, 0, 1, 2, 2 }, { 1, 2, 0, 1, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 2, 0 }, { 1, 0, 0, 0 }, { 0, 2, 0, 0 }, { 1, 2, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 345 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 6 }, { 4, 1, 5, 6 }, { 1, 2, 5, 6 }, { 1, 3, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 4, 3, 2, 0 }, { 5, 1, 1, 2, 2, 0 }, { 1, 4, 1, 2, 0, 0 }, { 2, 0, 1, 2, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 2, 0, 0, 1 }, { 0, 0, 0, 1 }, { 2, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );

                    break;
                }

                case ( 230 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 0, 3, 2, 2 }, { 0, 0, 0, 3, 3, 2 }, { 4, 0, 3, 2, 3, 2 }, { 1, 4, 5, 2, 3, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 2, 0, 0 }, { 3, 0, 0, 0 }, { 0, 0, 2, 0 }, { 3, 0, 2, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 302 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 0, 3, 2, 3, 0 }, { 2, 3, 3, 2, 2, 0 }, { 1, 3, 0, 0, 2, 0 }, { 5, 1, 4, 3, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 2, 0, 0, 3 }, { 0, 0, 0, 3 }, { 2, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 23 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 0, 0, 3 }, { 3, 2, 2, 0, 0, 3 }, { 5, 2, 2, 3, 0, 3 }, { 4, 5, 1, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 3, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 3, 2 }, { 0, 0, 3, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 521 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 2, 3, 1, 5, 2 }, { 1, 3, 3, 1, 1, 2 }, { 0, 3, 2, 2, 1, 2 }, { 4, 0, 3, 5, 1, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 2, 0, 3 }, { 1, 0, 0, 3 }, { 0, 0, 2, 3 }, { 1, 0, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 152 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 5, 2, 3, 1, 1 }, { 2, 2, 2, 3, 3, 1 }, { 3, 2, 5, 1, 3, 1 }, { 0, 3, 4, 1, 3, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 1, 0, 2 }, { 3, 0, 0, 2 }, { 0, 0, 1, 2 }, { 3, 0, 1, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 215 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 2, 2, 3 }, { 5, 1, 1, 2, 2, 3 }, { 4, 1, 1, 3, 2, 3 }, { 3, 4, 0, 2, 2, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 3, 0, 1 }, { 2, 0, 0, 1 }, { 0, 0, 3, 1 }, { 2, 0, 3, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 410 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 1 }, { 0, 3, 3, 0, 0, 1 }, { 2, 3, 1, 1, 0, 1 }, { 3, 2, 5, 4, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 1, 3 }, { 0, 0, 1, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 41 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 4, 1, 3, 0, 0 }, { 1, 1, 1, 3, 3, 0 }, { 5, 1, 4, 0, 3, 0 }, { 2, 5, 3, 0, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 3, 0, 0, 1 }, { 0, 0, 0, 1 }, { 3, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 104 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 1, 1, 3 }, { 4, 0, 0, 1, 1, 3 }, { 3, 0, 0, 3, 1, 3 }, { 5, 3, 2, 1, 1, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 3, 0, 0 }, { 1, 0, 0, 0 }, { 0, 0, 3, 0 }, { 1, 0, 3, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 543 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 4, 0, 2, 5, 1 }, { 3, 0, 0, 2, 2, 1 }, { 0, 0, 4, 1, 2, 1 }, { 2, 0, 1, 5, 2, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 1, 0, 0 }, { 2, 0, 0, 0 }, { 0, 0, 1, 0 }, { 2, 0, 1, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 354 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 5, 1, 0, 3, 2 }, { 4, 1, 1, 0, 0, 2 }, { 1, 1, 5, 2, 0, 2 }, { 0, 1, 2, 3, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 2, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 2, 1 }, { 0, 0, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 435 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 6 }, { 1, 4, 5, 6 }, { 2, 1, 5, 6 }, { 3, 1, 2, 6 } } );
                    mNumNewElem               = 4;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 2, 2 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 3, 2, 1, 4, 0 }, { 5, 2, 2, 1, 1, 0 }, { 2, 2, 3, 0, 1, 0 }, { 1, 2, 0, 4, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 1, 0, 0, 2 }, { 0, 0, 0, 2 }, { 1, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }

                default:
                {
                    std::cout << "WARNING UNDEFINED PERMUTATION" << std::endl;
                    break;
                }
            }
        }

        void
        hierarchy_tet4_4na( moris::size_t const &aPermutationId )
        {
            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;
            switch ( aPermutationId )
            {
                case ( 5420 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 5, 4, 0, 1 }, { 0, 0, 0, 0, 2, 1 }, { 0, 2, 0, 3, 2, 2 }, { 1, 1, 5, 2, 3, 2 }, { 0, 0, 1, 3, 3, 2 }, { 0, 1, 0, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 2, 1, 2 }, { 2, 2, 2, 0 }, { 3, 2, 2, 2 }, { 3, 2, 2, 1 }, { 3, 2, 2, 2 }, { 0, 1, 2, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 5, 4, 0, 1 }, { 0, 0, 0, 0, 2, 1 }, { 0, 2, 0, 3, 2, 2 }, { 1, 1, 5, 2, 3, 2 }, { 0, 0, 1, 3, 3, 2 }, { 0, 1, 0, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 1, 2 }, { 0, 0, 0, 0 }, { 3, 2, 0, 0 }, { 3, 0, 2, 1 }, { 3, 0, 0, 0 }, { 0, 1, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    break;
                }
                case ( 5240 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 5, 2, 3, 2 }, { 0, 3, 3, 0, 1, 2 }, { 0, 1, 0, 0, 4, 1 }, { 3, 2, 5, 4, 0, 1 }, { 0, 0, 2, 0, 0, 1 }, { 0, 2, 0, 3, 2, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 2, 1 }, { 0, 0, 0, 3 }, { 0, 1, 0, 0 }, { 0, 0, 1, 2 }, { 0, 0, 0, 0 }, { 3, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    break;
                }
                case ( 425 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 0, 4, 1, 0 }, { 5, 1, 1, 0, 3, 0 }, { 5, 3, 0, 2, 2, 3 }, { 3, 0, 0, 2, 2, 3 }, { 5, 0, 0, 2, 2, 3 }, { 5, 0, 0, 1, 4, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 0, 0, 0, 1 }, { 2, 3, 0, 0 }, { 2, 0, 3, 0 }, { 2, 0, 0, 0 }, { 1, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 3501 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 3, 5, 1, 2 }, { 1, 1, 1, 0, 0, 2 }, { 1, 0, 0, 3, 0, 0 }, { 2, 2, 3, 0, 3, 0 }, { 1, 0, 2, 3, 3, 0 }, { 1, 2, 0, 1, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 2, 0 }, { 0, 0, 0, 1 }, { 3, 0, 0, 0 }, { 3, 0, 0, 2 }, { 3, 0, 0, 0 }, { 1, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 1503 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 1, 5, 2, 1 }, { 3, 2, 2, 0, 3, 1 }, { 3, 3, 0, 0, 0, 3 }, { 4, 1, 1, 0, 0, 3 }, { 3, 0, 1, 0, 0, 3 }, { 3, 1, 0, 2, 5, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 1, 3 }, { 0, 0, 0, 2 }, { 0, 3, 0, 0 }, { 0, 0, 3, 1 }, { 0, 0, 0, 0 }, { 2, 1, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    break;
                }
                case ( 3051 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 3, 0, 3, 0 }, { 1, 3, 3, 0, 2, 0 }, { 1, 2, 0, 1, 5, 2 }, { 4, 0, 3, 5, 1, 2 }, { 1, 0, 0, 1, 1, 2 }, { 1, 0, 0, 3, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 0, 0, 0, 3 }, { 1, 2, 0, 0 }, { 1, 0, 2, 0 }, { 1, 0, 0, 0 }, { 3, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 1053 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 0, 3 }, { 3, 0, 0, 0, 1, 3 }, { 3, 1, 0, 2, 5, 1 }, { 2, 3, 1, 5, 2, 1 }, { 3, 0, 3, 2, 2, 1 }, { 3, 3, 0, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 3, 1 }, { 0, 0, 0, 0 }, { 2, 1, 0, 0 }, { 2, 0, 1, 3 }, { 2, 0, 0, 0 }, { 0, 3, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    break;
                }
                case ( 4312 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 4, 3, 2, 0 }, { 2, 2, 2, 0, 1, 0 }, { 2, 1, 0, 3, 1, 1 }, { 0, 0, 4, 1, 3, 1 }, { 2, 0, 0, 3, 3, 1 }, { 2, 0, 0, 2, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 0, 0, 0, 2 }, { 3, 1, 0, 0 }, { 3, 0, 1, 0 }, { 3, 0, 0, 0 }, { 2, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 2314 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 4, 0, 0, 0, 3, 2 }, { 4, 3, 0, 1, 1, 3 }, { 5, 2, 2, 1, 1, 3 }, { 4, 0, 2, 1, 1, 3 }, { 4, 2, 0, 0, 3, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 0, 0 }, { 1, 3, 0, 0 }, { 1, 0, 3, 2 }, { 1, 0, 0, 0 }, { 0, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 4132 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 4, 1, 3, 1 }, { 2, 3, 3, 0, 0, 1 }, { 2, 0, 0, 2, 3, 0 }, { 5, 1, 4, 3, 2, 0 }, { 2, 0, 1, 2, 2, 0 }, { 2, 1, 0, 3, 1, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 1, 0 }, { 0, 0, 0, 3 }, { 2, 0, 0, 0 }, { 2, 0, 0, 1 }, { 2, 0, 0, 0 }, { 3, 1, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 2134 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 1, 3 }, { 4, 1, 1, 0, 2, 3 }, { 4, 2, 0, 0, 3, 2 }, { 0, 3, 2, 3, 0, 2 }, { 4, 0, 3, 0, 0, 2 }, { 4, 3, 0, 1, 1, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 3, 2 }, { 0, 0, 0, 1 }, { 0, 2, 0, 0 }, { 0, 0, 2, 3 }, { 0, 0, 0, 0 }, { 1, 3, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );

                    break;
                }
                case ( 245 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 3, 7, 6 }, { 4, 3, 6, 7 }, { 4, 3, 7, 5 }, { 1, 2, 7, 5 }, { 2, 4, 7, 5 }, { 4, 2, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 2, 2, 3, 2, 2 }, { 1, 2, 3, 2, 1, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 3, 2, 2, 2, 2 }, { 1, 2, 3, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 2, 3 }, { 5, 2, 2, 0, 0, 3 }, { 5, 0, 0, 1, 4, 0 }, { 1, 3, 0, 4, 1, 0 }, { 5, 0, 3, 1, 1, 0 }, { 5, 3, 0, 2, 2, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 }, { 2, 3, 3, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 3, 0 }, { 0, 0, 0, 2 }, { 1, 0, 0, 0 }, { 1, 0, 0, 3 }, { 1, 0, 0, 0 }, { 2, 3, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 4502 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 4, 0, 2, 5, 1 }, { 2, 2, 2, 0, 0, 1 }, { 2, 0, 0, 0, 3, 0 }, { 1, 4, 1, 3, 0, 0 }, { 2, 1, 0, 3, 3, 0 }, { 2, 0, 1, 5, 2, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 1, 0, 0 }, { 0, 0, 0, 2 }, { 3, 0, 0, 0 }, { 3, 0, 0, 1 }, { 3, 0, 0, 0 }, { 2, 0, 1, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 4052 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 4, 1, 3, 0, 0 }, { 2, 3, 3, 1, 0, 0 }, { 2, 0, 1, 5, 2, 1 }, { 3, 4, 0, 2, 5, 1 }, { 2, 0, 0, 2, 2, 1 }, { 2, 0, 0, 0, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 0, 0, 0, 3 }, { 2, 0, 1, 0 }, { 2, 1, 0, 0 }, { 2, 0, 0, 0 }, { 3, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 1243 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 2, 2, 3 }, { 3, 2, 2, 1, 0, 3 }, { 3, 0, 1, 4, 0, 1 }, { 0, 1, 3, 0, 4, 1 }, { 3, 3, 0, 0, 0, 1 }, { 3, 0, 3, 2, 2, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 3, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 1, 0 }, { 0, 1, 0, 3 }, { 0, 0, 0, 0 }, { 2, 0, 3, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 2504 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 2, 3, 1, 5, 2 }, { 4, 1, 1, 3, 0, 2 }, { 4, 0, 3, 0, 0, 3 }, { 3, 2, 2, 0, 0, 3 }, { 4, 2, 0, 0, 0, 3 }, { 4, 0, 2, 5, 1, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 2, 0, 3 }, { 0, 0, 0, 1 }, { 0, 0, 3, 0 }, { 0, 3, 0, 2 }, { 0, 0, 0, 0 }, { 1, 0, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 2054 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 0, 0, 3 }, { 4, 0, 0, 2, 0, 3 }, { 4, 0, 2, 5, 1, 2 }, { 1, 2, 3, 1, 5, 2 }, { 4, 3, 0, 1, 1, 2 }, { 4, 0, 3, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 3, 0, 2 }, { 0, 0, 0, 0 }, { 1, 0, 2, 0 }, { 1, 2, 0, 3 }, { 1, 0, 0, 0 }, { 0, 0, 3, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );

                    break;
                }
                case ( 5310 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 5, 1, 0, 3, 2 }, { 0, 0, 0, 1, 0, 2 }, { 0, 0, 1, 1, 3, 1 }, { 2, 5, 2, 3, 1, 1 }, { 0, 2, 0, 3, 3, 1 }, { 0, 0, 2, 3, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 2, 0, 1 }, { 0, 0, 0, 0 }, { 3, 0, 1, 0 }, { 3, 1, 0, 2 }, { 3, 0, 0, 0 }, { 0, 0, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 5130 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 5, 2, 3, 1, 1 }, { 0, 3, 3, 2, 0, 1 }, { 0, 0, 2, 3, 0, 2 }, { 4, 5, 1, 0, 3, 2 }, { 0, 1, 0, 0, 0, 2 }, { 0, 0, 1, 1, 3, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 1, 0, 2 }, { 0, 0, 0, 3 }, { 0, 0, 2, 0 }, { 0, 2, 0, 1 }, { 0, 0, 0, 0 }, { 3, 0, 1, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );

                    break;
                }
                case ( 135 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 1, 1, 3 }, { 5, 1, 1, 0, 0, 3 }, { 5, 0, 0, 3, 2, 0 }, { 2, 0, 3, 2, 3, 0 }, { 5, 3, 0, 2, 2, 0 }, { 5, 0, 3, 1, 1, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 3, 0, 0 }, { 0, 0, 0, 1 }, { 2, 0, 0, 0 }, { 2, 0, 0, 3 }, { 2, 0, 0, 0 }, { 1, 0, 3, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );

                    break;
                }
                case ( 315 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 0, 3, 2, 3, 0 }, { 5, 2, 2, 3, 0, 0 }, { 5, 0, 3, 1, 1, 3 }, { 4, 0, 0, 1, 1, 3 }, { 5, 0, 0, 1, 1, 3 }, { 5, 0, 0, 3, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 0, 0, 0, 2 }, { 1, 0, 3, 0 }, { 1, 3, 0, 0 }, { 1, 0, 0, 0 }, { 2, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );

                    break;
                }
                case ( 3421 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 3, 2, 1, 4, 0 }, { 1, 1, 1, 2, 0, 0 }, { 1, 0, 2, 2, 3, 2 }, { 0, 3, 0, 3, 2, 2 }, { 1, 0, 0, 3, 3, 2 }, { 1, 0, 0, 4, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 0, 0, 0, 1 }, { 3, 0, 2, 0 }, { 3, 2, 0, 0 }, { 3, 0, 0, 0 }, { 1, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 3241 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 0, 3, 2, 2 }, { 1, 3, 3, 0, 0, 2 }, { 1, 0, 0, 4, 1, 0 }, { 5, 3, 2, 1, 4, 0 }, { 1, 2, 0, 1, 1, 0 }, { 1, 0, 2, 2, 3, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 2, 0, 0 }, { 0, 0, 0, 3 }, { 1, 0, 0, 0 }, { 1, 0, 0, 2 }, { 1, 0, 0, 0 }, { 3, 0, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                case ( 1423 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 3, 0, 7, 6 }, { 3, 4, 6, 7 }, { 3, 4, 7, 5 }, { 2, 1, 7, 5 }, { 4, 2, 7, 5 }, { 2, 4, 7, 6 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 1, 2, 2, 2, 3, 2 }, { 1, 3, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 2, 3, 2, 2, 2 }, { 1, 3, 2, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 1 }, { 3, 0, 0, 3, 0, 1 }, { 3, 0, 3, 2, 2, 3 }, { 5, 1, 1, 2, 2, 3 }, { 3, 1, 0, 2, 2, 3 }, { 3, 0, 1, 4, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 3, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 3, 3, 3 }, { 2, 3, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 0, 0, 0 }, { 2, 0, 3, 0 }, { 2, 3, 0, 1 }, { 2, 0, 0, 0 }, { 0, 0, 1, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 } } );
                    break;
                }
                default:
                    std::cout << "Template not found in the catalog" << std::endl;
                    break;
            }
        }

        void
        hierarchy_tet4_4nb( moris::size_t const &aPermutationId )
        {
            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;
            switch ( aPermutationId )
            {
                case ( 4250 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 3, 0, 4, 0 }, { 1, 2, 3, 4, 1, 0 }, { 5, 2, 2, 1, 1, 0 }, { 5, 2, 2, 1, 4, 0 }, { 3, 2, 2, 4, 0, 0 }, { 0, 3, 2, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3 }, { 1, 0, 0, 3 }, { 1, 0, 0, 2 }, { 1, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 2450 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3, 2, 0 }, { 3, 4, 0, 2, 2, 0 }, { 5, 1, 4, 2, 2, 0 }, { 5, 1, 1, 2, 2, 0 }, { 1, 4, 1, 2, 3, 0 }, { 0, 0, 4, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 0 }, { 2, 0, 0, 0 }, { 2, 0, 0, 1 }, { 2, 0, 0, 1 }, { 3, 0, 0, 1 }, { 3, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );

                    break;
                }
                case ( 4205 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 4, 0 }, { 3, 2, 2, 4, 0, 0 }, { 0, 3, 2, 0, 0, 0 }, { 0, 3, 3, 0, 4, 0 }, { 1, 2, 3, 4, 1, 0 }, { 5, 2, 2, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 3 }, { 0, 0, 0, 3 }, { 1, 0, 0, 3 }, { 1, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 2405 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 2, 2, 0 }, { 1, 4, 1, 2, 3, 0 }, { 0, 0, 4, 3, 3, 0 }, { 0, 0, 0, 3, 2, 0 }, { 3, 4, 0, 2, 2, 0 }, { 5, 1, 4, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 3, 0, 0, 1 }, { 3, 0, 0, 0 }, { 3, 0, 0, 0 }, { 2, 0, 0, 0 }, { 2, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 5031 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 3, 1, 5, 0 }, { 2, 0, 3, 5, 2, 0 }, { 3, 0, 0, 2, 2, 0 }, { 3, 0, 0, 2, 5, 0 }, { 4, 0, 0, 5, 1, 0 }, { 1, 3, 0, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 2, 0, 0, 3 }, { 2, 0, 0, 0 }, { 2, 0, 0, 0 }, { 1, 0, 0, 0 }, { 1, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 5013 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 5, 0 }, { 4, 0, 0, 5, 1, 0 }, { 1, 3, 0, 1, 1, 0 }, { 1, 3, 3, 1, 5, 0 }, { 2, 0, 3, 5, 2, 0 }, { 3, 0, 0, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 0 }, { 1, 0, 0, 0 }, { 1, 0, 0, 3 }, { 1, 0, 0, 3 }, { 2, 0, 0, 3 }, { 2, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 531 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 1, 3, 0, 0 }, { 4, 5, 1, 0, 0, 0 }, { 3, 2, 5, 0, 0, 0 }, { 3, 2, 2, 0, 0, 0 }, { 2, 5, 2, 0, 3, 0 }, { 1, 1, 5, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 0, 2 }, { 3, 0, 0, 2 }, { 3, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 513 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 0, 0, 0 }, { 2, 5, 2, 0, 3, 0 }, { 1, 1, 5, 3, 3, 0 }, { 1, 1, 1, 3, 0, 0 }, { 4, 5, 1, 0, 0, 0 }, { 3, 2, 5, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 2 }, { 3, 0, 0, 2 }, { 3, 0, 0, 1 }, { 3, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );
                    break;
                }
                case ( 3124 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 3, 0 }, { 5, 1, 1, 3, 2, 0 }, { 2, 3, 1, 2, 2, 0 }, { 2, 3, 3, 2, 3, 0 }, { 0, 1, 3, 3, 0, 0 }, { 4, 1, 1, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1 }, { 2, 0, 0, 1 }, { 2, 0, 0, 3 }, { 2, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );

                    break;
                }
                case ( 1342 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 2, 3, 1, 0 }, { 5, 3, 2, 1, 1, 0 }, { 4, 0, 3, 1, 1, 0 }, { 4, 0, 0, 1, 1, 0 }, { 0, 3, 0, 1, 3, 0 }, { 2, 2, 3, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 0, 0, 0 }, { 1, 0, 0, 0 }, { 3, 0, 0, 0 }, { 3, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );

                    break;
                }
                case ( 1324 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 1, 1, 0 }, { 0, 3, 0, 1, 3, 0 }, { 2, 2, 3, 3, 3, 0 }, { 2, 2, 2, 3, 1, 0 }, { 5, 3, 2, 1, 1, 0 }, { 4, 0, 3, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 0 }, { 3, 0, 0, 0 }, { 3, 0, 0, 2 }, { 3, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );

                    break;
                }
                case ( 3142 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 6, 7 }, { 0, 3, 6, 7 }, { 3, 5, 6, 7 }, { 5, 1, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 2, 1, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 3, 2, 3, 0 }, { 0, 1, 3, 3, 0, 0 }, { 4, 1, 1, 0, 0, 0 }, { 4, 1, 1, 0, 3, 0 }, { 5, 1, 1, 3, 2, 0 }, { 2, 3, 1, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 1 }, { 0, 0, 0, 1 }, { 2, 0, 0, 1 }, { 2, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 } } );

                    break;
                }
                case ( 5024 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 5, 1, 0 }, { 3, 0, 0, 2, 5, 0 }, { 2, 0, 3, 2, 2, 0 }, { 2, 3, 3, 5, 2, 0 }, { 1, 3, 0, 1, 5, 0 }, { 4, 0, 0, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 0 }, { 2, 0, 0, 0 }, { 2, 0, 0, 3 }, { 2, 0, 0, 3 }, { 1, 0, 0, 3 }, { 1, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                case ( 524 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 0, 0 }, { 1, 1, 5, 3, 0, 0 }, { 2, 5, 2, 3, 3, 0 }, { 2, 2, 2, 0, 3, 0 }, { 3, 2, 5, 0, 0, 0 }, { 4, 5, 1, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1 }, { 3, 0, 0, 1 }, { 3, 0, 0, 2 }, { 3, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 5042 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 3, 5, 2, 0 }, { 1, 3, 0, 1, 5, 0 }, { 4, 0, 0, 1, 1, 0 }, { 4, 0, 0, 5, 1, 0 }, { 3, 0, 0, 2, 5, 0 }, { 2, 0, 3, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 1, 0, 0, 3 }, { 1, 0, 0, 0 }, { 1, 0, 0, 0 }, { 2, 0, 0, 0 }, { 2, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                case ( 542 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 2, 0, 3, 0 }, { 3, 2, 5, 0, 0, 0 }, { 4, 5, 1, 0, 0, 0 }, { 4, 1, 1, 0, 0, 0 }, { 1, 1, 5, 3, 0, 0 }, { 2, 5, 2, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 1 }, { 0, 0, 0, 1 }, { 3, 0, 0, 1 }, { 3, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                case ( 3150 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 3, 3, 0, 0 }, { 2, 3, 1, 2, 3, 0 }, { 5, 1, 1, 2, 2, 0 }, { 5, 1, 1, 3, 2, 0 }, { 4, 1, 1, 0, 3, 0 }, { 0, 1, 3, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3 }, { 2, 0, 0, 3 }, { 2, 0, 0, 1 }, { 2, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                case ( 1350 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1, 3, 0 }, { 4, 0, 3, 1, 1, 0 }, { 5, 3, 2, 1, 1, 0 }, { 5, 2, 2, 1, 1, 0 }, { 2, 2, 3, 3, 1, 0 }, { 0, 3, 0, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 0 }, { 1, 0, 0, 0 }, { 1, 0, 0, 2 }, { 1, 0, 0, 2 }, { 3, 0, 0, 2 }, { 3, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                case ( 3105 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 3, 2, 0 }, { 4, 1, 1, 0, 3, 0 }, { 0, 1, 3, 0, 0, 0 }, { 0, 3, 3, 3, 0, 0 }, { 2, 3, 1, 2, 3, 0 }, { 5, 1, 1, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 3 }, { 0, 0, 0, 3 }, { 2, 0, 0, 3 }, { 2, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 1305 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 1, 0 }, { 2, 2, 3, 3, 1, 0 }, { 0, 3, 0, 3, 3, 0 }, { 0, 0, 0, 1, 3, 0 }, { 4, 0, 3, 1, 1, 0 }, { 5, 3, 2, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 3, 0, 0, 2 }, { 3, 0, 0, 0 }, { 3, 0, 0, 0 }, { 1, 0, 0, 0 }, { 1, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 4231 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 3, 4, 1, 0 }, { 0, 3, 2, 0, 4, 0 }, { 3, 2, 2, 0, 0, 0 }, { 3, 2, 2, 4, 0, 0 }, { 5, 2, 2, 1, 4, 0 }, { 1, 2, 3, 1, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 2 }, { 0, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 2431 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 1, 2, 3, 0 }, { 5, 1, 4, 2, 2, 0 }, { 3, 4, 0, 2, 2, 0 }, { 3, 0, 0, 2, 2, 0 }, { 0, 0, 4, 3, 2, 0 }, { 1, 4, 1, 3, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 2, 0, 0, 1 }, { 2, 0, 0, 0 }, { 2, 0, 0, 0 }, { 3, 0, 0, 0 }, { 3, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 4213 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 4, 0, 0 }, { 5, 2, 2, 1, 4, 0 }, { 1, 2, 3, 1, 1, 0 }, { 1, 3, 3, 4, 1, 0 }, { 0, 3, 2, 0, 4, 0 }, { 3, 2, 2, 0, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 0, 0, 3 }, { 1, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );
                    break;
                }
                case ( 2413 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 6, 7 }, { 3, 0, 6, 7 }, { 5, 3, 6, 7 }, { 1, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 1, 2, 2, 2, 3 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 2, 0 }, { 0, 0, 4, 3, 2, 0 }, { 1, 4, 1, 3, 3, 0 }, { 1, 1, 1, 2, 3, 0 }, { 5, 1, 4, 2, 2, 0 }, { 3, 4, 0, 2, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 0 }, { 3, 0, 0, 0 }, { 3, 0, 0, 1 }, { 3, 0, 0, 1 }, { 2, 0, 0, 1 }, { 2, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 } } );

                    break;
                }
                default:
                {
                    std::cout << "Template not found in the catalog" << std::endl;
                    break;
                }
            }
        }

        void
        hierarchy_tet4_4nc( moris::size_t const &aPermutationId )
        {
            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;
            switch ( aPermutationId )
            {
                case ( 4520 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 3, 0, 4, 0 }, { 1, 2, 3, 4, 1, 0 }, { 2, 5, 2, 0, 1, 1 }, { 3, 2, 5, 4, 0, 1 }, { 0, 3, 2, 0, 0, 0 }, { 2, 2, 2, 0, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3 }, { 1, 0, 0, 3 }, { 0, 1, 0, 2 }, { 0, 0, 1, 2 }, { 0, 0, 0, 3 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 2540 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3, 2, 0 }, { 3, 4, 0, 2, 2, 0 }, { 4, 5, 1, 0, 2, 2 }, { 1, 1, 5, 2, 3, 2 }, { 0, 0, 4, 3, 3, 0 }, { 4, 1, 1, 3, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 0 }, { 2, 0, 0, 0 }, { 0, 2, 0, 1 }, { 3, 0, 2, 1 }, { 3, 0, 0, 0 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 4025 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 4, 0 }, { 3, 2, 2, 4, 0, 0 }, { 2, 0, 3, 0, 0, 0 }, { 1, 3, 0, 4, 1, 0 }, { 5, 2, 2, 1, 1, 0 }, { 2, 3, 3, 1, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 3 }, { 1, 0, 0, 3 }, { 1, 0, 0, 2 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 2045 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 2, 2, 0 }, { 1, 4, 1, 2, 3, 0 }, { 4, 0, 0, 0, 3, 3 }, { 3, 0, 0, 2, 2, 3 }, { 5, 1, 4, 2, 2, 0 }, { 4, 0, 0, 2, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 3, 0, 0, 1 }, { 0, 3, 0, 0 }, { 2, 0, 3, 0 }, { 2, 0, 0, 1 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 5301 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 3, 1, 5, 0 }, { 2, 0, 3, 5, 2, 0 }, { 0, 3, 0, 0, 2, 2 }, { 4, 0, 3, 5, 1, 2 }, { 1, 3, 0, 1, 1, 0 }, { 0, 0, 0, 1, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 2, 0, 0, 3 }, { 0, 2, 0, 0 }, { 1, 0, 2, 0 }, { 1, 0, 0, 3 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 5103 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 5, 0 }, { 4, 0, 0, 5, 1, 0 }, { 0, 1, 3, 0, 1, 1 }, { 2, 3, 1, 5, 2, 1 }, { 3, 0, 0, 2, 2, 0 }, { 0, 3, 3, 2, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 0 }, { 1, 0, 0, 0 }, { 0, 1, 0, 3 }, { 2, 0, 1, 3 }, { 2, 0, 0, 0 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 351 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 1, 3, 0, 0 }, { 4, 5, 1, 0, 0, 0 }, { 5, 3, 2, 0, 0, 0 }, { 2, 2, 3, 0, 3, 0 }, { 1, 1, 5, 3, 3, 0 }, { 5, 2, 2, 3, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 3, 0, 0, 2 }, { 3, 0, 0, 1 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 153 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 0, 0, 0 }, { 2, 5, 2, 0, 3, 0 }, { 5, 1, 1, 0, 3, 3 }, { 4, 1, 1, 0, 0, 3 }, { 3, 2, 5, 0, 0, 0 }, { 5, 1, 1, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 2 }, { 3, 0, 0, 2 }, { 0, 3, 0, 1 }, { 0, 0, 3, 1 }, { 0, 0, 0, 2 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 3412 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 3, 2, 3, 0 }, { 0, 1, 3, 3, 0, 0 }, { 1, 4, 1, 0, 0, 0 }, { 5, 1, 4, 3, 2, 0 }, { 2, 3, 1, 2, 2, 0 }, { 1, 1, 1, 2, 0, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 1 }, { 2, 0, 0, 1 }, { 2, 0, 0, 3 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );

                    break;
                }
                case ( 3214 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 3, 0 }, { 5, 1, 1, 3, 2, 0 }, { 1, 2, 3, 0, 2, 2 }, { 0, 3, 2, 3, 0, 2 }, { 4, 1, 1, 0, 0, 0 }, { 1, 3, 3, 0, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1 }, { 2, 0, 0, 1 }, { 0, 2, 0, 3 }, { 0, 0, 2, 3 }, { 0, 0, 0, 1 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 1432 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 2, 3, 1, 0 }, { 5, 3, 2, 1, 1, 0 }, { 3, 4, 0, 0, 1, 1 }, { 0, 0, 4, 1, 3, 1 }, { 2, 2, 3, 3, 3, 0 }, { 3, 0, 0, 3, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 1, 0, 0, 2 }, { 0, 1, 0, 0 }, { 3, 0, 1, 0 }, { 3, 0, 0, 2 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 1234 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 4, 0, 5, 7 }, { 0, 3, 5, 7 }, { 5, 3, 6, 7 }, { 1, 2, 6, 7 }, { 2, 4, 5, 7 }, { 2, 5, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 3 }, { 1, 1, 2, 1, 2, 3 }, { 1, 1, 2, 3, 2, 2 }, { 1, 2, 1, 1, 2, 2 }, { 1, 2, 1, 2, 2, 3 }, { 1, 2, 2, 2, 3, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 1, 1, 0 }, { 0, 3, 0, 1, 3, 0 }, { 3, 2, 2, 0, 3, 3 }, { 5, 2, 2, 1, 1, 3 }, { 4, 0, 3, 1, 1, 0 }, { 3, 2, 2, 1, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 0 }, { 3, 0, 0, 0 }, { 0, 3, 0, 2 }, { 1, 0, 3, 2 }, { 1, 0, 0, 0 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 1 } } );
                    break;
                }
                case ( 5402 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 3, 3, 5, 2, 0 }, { 1, 3, 0, 1, 5, 0 }, { 0, 0, 4, 1, 0, 1 }, { 3, 4, 0, 2, 5, 1 }, { 2, 0, 3, 2, 2, 0 }, { 0, 0, 0, 0, 2, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 3 }, { 1, 0, 0, 3 }, { 0, 0, 1, 0 }, { 2, 1, 0, 0 }, { 2, 0, 0, 3 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 452 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 2, 2, 2, 0, 3, 0 }, { 3, 2, 5, 0, 0, 0 }, { 5, 1, 4, 0, 0, 0 }, { 1, 4, 1, 3, 0, 0 }, { 2, 5, 2, 3, 3, 0 }, { 5, 1, 1, 0, 3, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2 }, { 0, 0, 0, 2 }, { 0, 0, 0, 1 }, { 3, 0, 0, 1 }, { 3, 0, 0, 2 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 5204 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 0, 0, 5, 1, 0 }, { 3, 0, 0, 2, 5, 0 }, { 0, 3, 2, 2, 0, 2 }, { 1, 2, 3, 1, 5, 2 }, { 4, 0, 0, 1, 1, 0 }, { 0, 3, 3, 0, 1, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 0 }, { 2, 0, 0, 0 }, { 0, 0, 2, 3 }, { 1, 2, 0, 3 }, { 1, 0, 0, 0 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 254 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 4, 1, 1, 0, 0, 0 }, { 1, 1, 5, 3, 0, 0 }, { 5, 2, 2, 3, 0, 3 }, { 3, 2, 2, 0, 0, 3 }, { 4, 5, 1, 0, 0, 0 }, { 5, 2, 2, 0, 0, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1 }, { 3, 0, 0, 1 }, { 0, 0, 3, 2 }, { 0, 3, 0, 2 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );
                    break;
                }
                case ( 3510 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 3, 3, 0, 0 }, { 2, 3, 1, 2, 3, 0 }, { 1, 1, 5, 2, 0, 2 }, { 4, 5, 1, 0, 3, 2 }, { 0, 1, 3, 0, 0, 0 }, { 1, 1, 1, 0, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3 }, { 2, 0, 0, 3 }, { 0, 0, 2, 1 }, { 0, 2, 0, 1 }, { 0, 0, 0, 3 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 1530 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 1, 3, 0 }, { 4, 0, 3, 1, 1, 0 }, { 3, 2, 5, 1, 0, 1 }, { 2, 5, 2, 3, 1, 1 }, { 0, 3, 0, 3, 3, 0 }, { 3, 2, 2, 0, 3, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 0 }, { 1, 0, 0, 0 }, { 0, 0, 1, 2 }, { 3, 1, 0, 2 }, { 3, 0, 0, 0 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 3015 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 1, 1, 3, 2, 0 }, { 4, 1, 1, 0, 3, 0 }, { 1, 3, 0, 0, 0, 0 }, { 2, 0, 3, 2, 3, 0 }, { 5, 1, 1, 2, 2, 0 }, { 1, 3, 3, 0, 2, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 1 }, { 0, 0, 0, 1 }, { 0, 0, 0, 3 }, { 2, 0, 0, 3 }, { 2, 0, 0, 1 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );
                    break;
                }
                case ( 1035 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 5, 2, 2, 1, 1, 0 }, { 2, 2, 3, 3, 1, 0 }, { 3, 0, 0, 3, 0, 3 }, { 4, 0, 0, 1, 1, 3 }, { 5, 3, 2, 1, 1, 0 }, { 3, 0, 0, 0, 1, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 2 }, { 3, 0, 0, 2 }, { 0, 0, 3, 0 }, { 1, 3, 0, 0 }, { 1, 0, 0, 2 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 4321 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 3, 3, 4, 1, 0 }, { 0, 3, 2, 0, 4, 0 }, { 2, 2, 3, 0, 0, 0 }, { 5, 3, 2, 1, 4, 0 }, { 1, 2, 3, 1, 1, 0 }, { 2, 2, 2, 0, 1, 0 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 0, 0, 3 }, { 0, 0, 0, 3 }, { 0, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 0, 0, 3 }, { 0, 0, 0, 2 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );
                    break;
                }
                case ( 2341 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 1, 1, 2, 3, 0 }, { 5, 1, 4, 2, 2, 0 }, { 4, 0, 3, 2, 0, 2 }, { 0, 3, 0, 3, 2, 2 }, { 1, 4, 1, 3, 3, 0 }, { 4, 0, 0, 0, 3, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 1 }, { 2, 0, 0, 1 }, { 0, 0, 2, 0 }, { 3, 2, 0, 0 }, { 3, 0, 0, 1 }, { 0, 0, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 4123 ):
                {

                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 2, 2, 4, 0, 0 }, { 5, 2, 2, 1, 4, 0 }, { 2, 3, 1, 1, 0, 1 }, { 0, 1, 3, 0, 4, 1 }, { 3, 2, 2, 0, 0, 0 }, { 2, 3, 3, 0, 0, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 2 }, { 1, 0, 0, 2 }, { 0, 0, 1, 3 }, { 0, 1, 0, 3 }, { 0, 0, 0, 2 }, { 0, 0, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );

                    break;
                }
                case ( 2143 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 5, 7 }, { 3, 0, 5, 7 }, { 3, 5, 6, 7 }, { 2, 1, 6, 7 }, { 4, 2, 5, 7 }, { 5, 2, 6, 7 } } );
                    mNumNewElem               = 6;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 2, 3 }, { 1, 2, 1, 2, 1, 3 }, { 1, 2, 1, 2, 3, 2 }, { 1, 1, 2, 2, 1, 2 }, { 1, 1, 2, 2, 2, 3 }, { 1, 2, 2, 3, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 0, 0, 2, 2, 0 }, { 0, 0, 4, 3, 2, 0 }, { 4, 1, 1, 3, 0, 3 }, { 5, 1, 1, 2, 2, 3 }, { 3, 4, 0, 2, 2, 0 }, { 4, 1, 1, 0, 2, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 }, { 2, 3, 3, 2 }, { 3, 3, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 2, 0, 0, 0 }, { 3, 0, 0, 0 }, { 0, 0, 3, 1 }, { 2, 3, 0, 1 }, { 2, 0, 0, 0 }, { 0, 0, 0, 1 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 2 } } );
                    break;
                }
                default:
                {
                    std::cout << "Template not found in the catalog" << std::endl;
                    break;
                }
            }
        }

        void
        bisected_tet( moris::size_t const &aPermutationId )
        {
            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;

            switch ( aPermutationId )
            {
                case 0:
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 2, 3 },
                                    { 4, 1, 2, 3 } } );
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 1 }, { 1, 1, 2, 2, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 5 }, { 0, 1, 3, 0, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 } } );
                    break;
                }
                case 1:
                {
                    mNewElementToNode      = Matrix< IndexMat >( { { 0, 1, 4, 3 },
                                 { 0, 4, 2, 3 } } );
                    mNumNewElem            = 2;
                    mNumElemToReplace      = 1;
                    mNewParentEdgeRanks    = Matrix< moris::DDSTMat >( { { 1, 1, 2, 1, 1, 2 }, { 2, 1, 1, 1, 2, 1 } } );
                    mNewParentEdgeOrdinals = Matrix< IndexMat >( { { 0, 1, 3, 3, 4, 1 }, { 3, 1, 2, 3, 1, 5 } } );
                    mNewParentFaceRanks =
                            Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals =
                            Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 0 } } );
                    break;
                }
                case 2:
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 4, 3 },
                                    { 4, 1, 2, 3 } } );
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 1, 2 }, { 2, 1, 1, 2, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 4, 2 }, { 3, 1, 2, 2, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 0 } } );
                    break;
                }
                case 3:
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 4 },
                                    { 4, 1, 2, 3 } } );
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 1, 2, 2 }, { 2, 1, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 0, 2 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 3 } } );
                    break;
                }
                case 4:
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 2, 3 },
                                    { 0, 1, 2, 4 } } );
                    mNumNewElem               = 2;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 1, 1, 1, 1 }, { 1, 1, 1, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 2, 3 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 0 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 3 }, { 2 } } );
                    break;
                }
                case 5:
                {
                    mNewElementToNode      = Matrix< IndexMat >( { { 0, 1, 4, 3 },
                                 { 0, 1, 2, 4 } } );
                    mNumNewElem            = 2;
                    mNumElemToReplace      = 1;
                    mNewParentEdgeRanks    = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 1, 1 }, { 1, 1, 1, 2, 2, 1 } } );
                    mNewParentEdgeOrdinals = Matrix< IndexMat >( { { 0, 1, 2, 3, 4, 5 }, { 0, 1, 2, 2, 1, 5 } } );
                    mNewParentFaceRanks =
                            Matrix< moris::DDSTMat >( { { 2, 2, 2, 3 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals =
                            Matrix< IndexMat >( { { 0, 1, 2, 0 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 3 }, { 0 } } );
                    break;
                }
                default:
                {
                    MORIS_ERROR( 0, "Invalid edge ordinal found" );
                    break;
                }
            }
        }

        void
        hierarchy_tet4_2( moris::size_t const &aPermutationId )
        {
            mSpatialDimension = 3;
            mElementTopology  = moris::mtk::CellTopology::TET4;

            switch ( aPermutationId )
            {
                case ( 14 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 4, 5 }, { 0, 4, 2, 5 }, { 0, 5, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 2, 1, 1, 2, 2, 2 }, { 2, 2, 1, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 1 }, { 3, 1, 2, 0, 1, 1 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 2, 3, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 0, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 0 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }

                case ( 41 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 4 }, { 0, 4, 5, 3 }, { 0, 5, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 2, 2, 2, 1, 1, 2 }, { 2, 1, 1, 1, 2, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 1 }, { 0, 1, 3, 3, 4, 1 }, { 3, 1, 2, 3, 1, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 2, 3, 3 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 0, 0 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 3 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 15 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 4, 5 }, { 0, 1, 5, 3 }, { 0, 4, 2, 5 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 2, 2 }, { 1, 2, 2, 1, 1, 1 }, { 2, 1, 1, 2, 2, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 2, 1, 1 }, { 0, 1, 2, 3, 4, 5 }, { 3, 1, 2, 2, 1, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 2, 3, 2 }, { 2, 2, 2, 3 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 2, 0 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 0 } } );
                    break;
                }
                case ( 51 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 3 }, { 0, 5, 4, 3 }, { 0, 5, 2, 4 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 1, 1, 2 }, { 2, 2, 2, 1, 2, 1 }, { 2, 1, 1, 2, 2, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 3, 4, 1 }, { 3, 1, 2, 3, 1, 5 }, { 3, 1, 2, 2, 1, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 2, 2, 3 }, { 3, 2, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 2, 0 }, { 0, 1, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 3 }, { 0 } } );
                    break;
                }
                case ( 45 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 4 }, { 0, 1, 2, 5 }, { 0, 4, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 2, 1, 2 }, { 1, 1, 1, 2, 2, 1 }, { 2, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 0, 4, 1 }, { 0, 1, 2, 2, 1, 5 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 3 }, { 3, 2, 2, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 0 }, { 0, 1, 2, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { std::numeric_limits< moris::size_t >::max() }, { 3 } } );
                    break;
                }
                case ( 54 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 5 }, { 0, 5, 2, 4 }, { 0, 5, 4, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 2, 1, 2 }, { 2, 2, 1, 2, 2, 1 }, { 2, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 0, 4, 1 }, { 0, 1, 2, 2, 1, 5 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 2, 2, 3 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 2, 0 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 0 }, { 3 } } );
                    break;
                }
                case ( 53 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 4, 3 }, { 1, 2, 5, 4 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 1, 2, 2 }, { 2, 2, 2, 1, 1, 1 }, { 1, 2, 2, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 0, 2 }, { 0, 1, 2, 3, 4, 5 }, { 1, 2, 0, 1, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 2, 2, 3 }, { 2, 2, 3, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 2, 0 }, { 1, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 3 }, { 2 } } );
                    break;
                }
                case ( 35 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 5 }, { 0, 1, 5, 4 }, { 4, 1, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 2, 2, 1 }, { 1, 2, 2, 1, 2, 2 }, { 2, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 2, 1, 5 }, { 0, 1, 2, 3, 0, 2 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 2, 2, 2 }, { 2, 3, 2, 3 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3 }, { 0, 0, 2, 0 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 3 } } );
                    break;
                }
                case ( 12 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 5, 4, 2, 3 }, { 1, 4, 5, 3 }, { 0, 1, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 1, 1, 2, 2, 1 }, { 1, 2, 2, 1, 2, 2 }, { 1, 2, 1, 1, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 1, 2, 2, 1, 5 }, { 1, 3, 3, 4, 1, 2 }, { 0, 3, 2, 3, 4, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 2, 2, 2 }, { 2, 3, 3, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { {
                                                                           0,
                                                                           1,
                                                                           2,
                                                                           3,
                                                                   },
                               { 1, 0, 0, 3 },
                               { 0, 0, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 0 }, { 1 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 21 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 5, 2, 4, 3 }, { 0, 5, 4, 3 }, { 0, 1, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 2, 2, 1, 1, 2, 2 }, { 1, 1, 2, 1, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 1, 2, 3, 1, 5, 2 }, { 3, 3, 2, 3, 1, 2 }, { 0, 1, 3, 3, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 1, 2, 0, 3 }, { 0, 0, 2, 3 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 1 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 02 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 5, 1, 2, 3 }, { 5, 4, 1, 3 }, { 0, 4, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 1, 1, 2, 1, 1 }, { 2, 1, 2, 2, 2, 1 }, { 1, 2, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 1, 2, 2, 4, 5 }, { 3, 0, 3, 2, 0, 4 }, { 0, 3, 2, 3, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 2, 2, 2 }, { 3, 2, 3, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3 }, { 0, 0, 0, 3 }, { 0, 0, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 0 }, { 1 } } );
                    break;
                }
                case ( 20 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 5, 1, 2, 3 }, { 4, 5, 2, 3 }, { 0, 5, 4, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 1 }, { 2, 2, 1, 2, 2, 1 }, { 1, 2, 1, 1, 2, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 5 }, { 3, 3, 2, 2, 0, 5 }, { 0, 3, 2, 3, 0, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 3, 3, 2, 2 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 0, 2, 3 }, { 0, 0, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 0 }, { 1 } } );
                    break;
                }
                case ( 01 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 5, 2, 3 }, { 0, 4, 5, 3 }, { 4, 1, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 1, 1, 1, 2, 1 }, { 1, 2, 2, 1, 2, 2 }, { 1, 1, 2, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 3, 1, 2, 3, 1, 5 }, { 0, 3, 3, 3, 0, 1 }, { 0, 1, 3, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 3, 2, 2, 2 }, { 2, 3, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3 }, { 0, 0, 0, 3 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 10 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 5, 2, 3 }, { 5, 4, 2, 3 }, { 5, 1, 4, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 1 }, { 2, 1, 2, 2, 2, 1 }, { 1, 1, 2, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 5 }, { 3, 1, 3, 0, 1, 5 }, { 0, 1, 3, 0, 4, 1 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 2, 3, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 0, 3 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 0 }, { 2 } } );
                    break;
                }
                case ( 03 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 5, 1, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 1, 1, 2, 2, 2, 2 }, { 2, 1, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 0, 1, 3, 0, 0, 2 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 3, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 0, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 30 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 5, 2, 4 }, { 4, 5, 2, 3 }, { 5, 1, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 2, 2, 2, 1, 2, 1 }, { 1, 1, 2, 2, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 0, 3, 2, 3, 0, 5 }, { 0, 1, 3, 0, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 2, 0 }, { 0, 1, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 04 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 4, 2, 5 }, { 4, 1, 2, 5 }, { 0, 5, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 2, 2, 2 }, { 1, 1, 2, 2, 1, 2 }, { 2, 2, 1, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 0, 0, 1 }, { 0, 1, 3, 0, 4, 1 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 3, 2 }, { 2, 2, 3, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 0, 3 }, { 0, 1, 0, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 2 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 40 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 5, 1, 2, 4 }, { 5, 4, 2, 3 }, { 0, 5, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 2, 2, 1, 2 }, { 2, 2, 2, 2, 1, 1 }, { 1, 2, 1, 1, 2, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 3, 0, 4, 1 }, { 0, 1, 3, 0, 4, 5 }, { 0, 3, 2, 3, 0, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 3, 2 }, { 2, 2, 3, 3 }, { 2, 3, 2, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 0, 3 }, { 0, 1, 0, 0 }, { 0, 0, 2, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 2 }, { 3 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 34 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 5, 2, 4 }, { 0, 1, 2, 5 }, { 4, 5, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 1, 1, 2, 2 }, { 1, 1, 1, 2, 1, 2 }, { 2, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 2, 2, 3, 0, 2 }, { 0, 1, 2, 0, 4, 1 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 3 }, { 2, 2, 3, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 0 }, { 0, 1, 0, 3 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { std::numeric_limits< moris::size_t >::max() }, { 3 } } );
                    break;
                }
                case ( 43 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 2, 5 }, { 5, 1, 2, 4 }, { 5, 4, 2, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 1, 2, 2 }, { 2, 1, 2, 2, 1, 2 }, { 2, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 0, 2 }, { 0, 1, 2, 0, 4, 1 }, { 0, 1, 2, 3, 4, 5 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 2, 3, 3 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 1, 0, 0 }, { 0, 1, 2, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 2 }, { 3 } } );
                    break;
                }
                case ( 25 ):
                {
                    // done
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 3 }, { 0, 1, 4, 5 }, { 1, 2, 4, 5 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 2, 1, 1, 1 }, { 1, 2, 1, 2, 2, 2 }, { 1, 1, 2, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 3, 4, 5 }, { 0, 3, 2, 2, 1, 2 }, { 1, 2, 3, 1, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 2, 2, 3 }, { 3, 3, 2, 2 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 1, 2, 0 }, { 0, 0, 2, 3 }, { 1, 2, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 1 }, { 2 } } );
                    break;
                }
                case ( 52 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 3 }, { 1, 4, 5, 3 }, { 1, 2, 5, 4 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 1, 2 }, { 2, 2, 2, 1, 1, 2 }, { 1, 1, 2, 2, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 4, 2 }, { 1, 2, 3, 4, 5, 2 }, { 1, 2, 3, 1, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 2, 3, 3 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 1, 2, 0, 0 }, { 1, 2, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { std::numeric_limits< moris::size_t >::max() }, { 3 }, { 2 } } );
                    break;
                }
                case ( 23 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 4, 5 }, { 4, 1, 2, 5 }, { 1, 2, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 1, 1, 1, 2, 2 }, { 2, 1, 1, 2, 2, 2 }, { 1, 2, 2, 1, 1, 1 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 3, 1, 2, 2, 0, 2 }, { 1, 2, 0, 4, 5, 3 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 3, 3, 2, 2 }, { 2, 2, 2, 3 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 2, 3 }, { 1, 2, 0, 0 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 0 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                case ( 32 ):
                {
                    mNewElementToNode         = Matrix< IndexMat >( { { 0, 1, 5, 4 }, { 4, 1, 5, 3 }, { 1, 2, 5, 3 } } );
                    mNumNewElem               = 3;
                    mNumElemToReplace         = 1;
                    mNewParentEdgeRanks       = Matrix< moris::DDSTMat >( { { 1, 2, 1, 1, 2, 2 }, { 2, 2, 2, 1, 1, 2 }, { 1, 1, 2, 1, 1, 2 } } );
                    mNewParentEdgeOrdinals    = Matrix< IndexMat >( { { 0, 3, 2, 3, 0, 2 }, { 0, 3, 2, 3, 4, 2 }, { 1, 2, 3, 4, 5, 2 } } );
                    mNewParentFaceRanks       = Matrix< moris::DDSTMat >( { { 2, 3, 2, 2 }, { 2, 3, 2, 3 }, { 2, 2, 3, 2 } } );
                    mNewParentFaceOrdinals    = Matrix< IndexMat >( { { 0, 0, 2, 3 }, { 0, 0, 2, 0 }, { 1, 2, 0, 3 } } );
                    mNewElementInterfaceSides = Matrix< moris::DDSTMat >( { { 1 }, { 3 }, { std::numeric_limits< moris::size_t >::max() } } );
                    break;
                }
                //            mNewElementToNode = Matrix< IndexMat >({{},{},{}});
                //            mNumNewElem = 3;
                //            mNumElemToReplace = 1;
                //            mNewParentEdgeRanks = Matrix< moris::DDSTMat >({{},{},{}});
                //            mNewParentEdgeOrdinals = Matrix< IndexMat >({{},{},{}});
                //            mNewParentFaceRanks = Matrix< moris::DDSTMat >({{},{},{}});
                //            mNewParentFaceOrdinals = Matrix< IndexMat >({{},{},{}});
                //            mNewElementInterfaceSides = Matrix< moris::DDSTMat >({{},{},{}});
                default:
                {
                    std::cout << "Template not found in the catalog" << std::endl;
                    break;
                }
            }
        }
    };

}    // namespace moris::xtk
#endif /* SRC_XTK_CL_XTK_CHILD_MESH_MODIFICATION_TEMPLATE_HPP_ */
