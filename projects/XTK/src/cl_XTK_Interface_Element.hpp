/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Interface_Element.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_

#include "cl_Matrix.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "op_div.hpp"
#include "fn_trans.hpp"

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "assert.hpp"

#include "cl_MTK_Cell_Info_Tet4.hpp"

// TODO:PUT THIS IN A MORE APPROPRIATE SPOT
inline Matrix< IndexMat >
get_nodes_on_face_map( moris::uint aNumVertsPerElem,
        moris::mtk::Geometry_Type  aGeometryType,
        moris::uint                aSideOrd )
{
    if ( aNumVertsPerElem == 4 && aGeometryType == moris::mtk::Geometry_Type::TET )
    {
        moris::mtk::Cell_Info_Tet4 tTet4;
        return tTet4.get_node_to_face_map( aSideOrd );
    }
    else
    {
        MORIS_ERROR( 0, "Unsupported node on face map" );
        return Matrix< IndexMat >( 0, 0 );
    }
}

namespace moris::xtk
{

    class Interface_Element
    {
      public:
        Interface_Element()
                : mElementId( MORIS_ID_MAX )
                , mElementIndex( MORIS_INDEX_MAX )
                , mElementOwner( MORIS_INDEX_MAX )
                , mElementPairs( 0 )
                , mElementSideOrdinals( 0, 0 )
        {
        }

        /*
         * Provide this interface element with the element index pair and corresponding side ordinals which these two elements
         * share
         */
        void
        set_element_pair_and_side_ordinal(
                const Vector< const moris::mtk::Cell* >& aElementIndexPair,
                const Matrix< IndexMat >&                aElementSideOrdinals )
        {
            // Verify lengths
            MORIS_ASSERT( aElementIndexPair.size() == 2, "Provided aElementIndexPair vector needs to have length 2" );
            MORIS_ASSERT( aElementSideOrdinals.numel() == 2, "Provided aElementSideOrdinals vector needs to have length 2" );

            // Verify that we haven't already set these
            MORIS_ASSERT( mElementPairs.size() == 0, "This interface element has already been given an element index pair" );
            MORIS_ASSERT( mElementSideOrdinals.numel() == 0, "This interface element has already been given an element side ordinal pair" );

            // Copy the vectors
            mElementPairs        = aElementIndexPair;
            mElementSideOrdinals = aElementSideOrdinals.copy();
        }

        /*!
         * Set the element identifiers (index and id)
         *
         */
        void
        set_element_identifiers(
                moris::moris_index aElementIndex,
                moris::moris_id    aElementId )
        {
            MORIS_ASSERT( mElementId == MORIS_ID_MAX, "Element Id has already been set" );
            MORIS_ASSERT( mElementIndex == MORIS_INDEX_MAX, "Element index has already been set" );

            mElementId    = aElementId;
            mElementIndex = aElementIndex;
        }

        void
        set_element_owner( moris::moris_index aElementOwner )
        {
            MORIS_ASSERT( mElementOwner == MORIS_INDEX_MAX, "Owner already set" );
            mElementOwner = aElementOwner;
        }

        moris::moris_id
        get_element_id()
        {
            return mElementId;
        }

        moris::moris_index
        get_element_index()
        {
            return mElementIndex;
        }

        moris::moris_index
        get_element_owner()
        {
            return mElementOwner;
        }

        /*
         * Returns the cell pointers in this interface element
         */
        Vector< const moris::mtk::Cell* > const &
        get_cell_ptrs()
        {
            MORIS_ASSERT( mElementPairs.size() == 2, "This interface element has not been given an element index pair" );
            return mElementPairs;
        }

        /*!
         * returns the side ordinals of the element pair on the interface
         */
        Matrix< IndexMat > const &
        get_element_pair_side_ordinals() const
        {
            MORIS_ASSERT( mElementSideOrdinals.numel() == 2, "This interface element has not been given an element index pair" );
            return mElementSideOrdinals;
        }

        /*!
         * return the outward facing normal relative to the provided pair index
         * @param[in] - aPairIndex
         */
        Matrix< moris::F31RMat >
        get_outward_normal( moris::uint aPairIndex ) const
        {
            MORIS_ASSERT( aPairIndex < mElementPairs.size(), "Pair Index provided is out of bounds" );

            // Vertex coordinates connected to cell
            Matrix< DDRMat > tVertexCoords = mElementPairs( aPairIndex )->get_vertex_coords();

            // Get the nodes which need to be used to compute normal
            moris::mtk::Cell_Info_Tet4 tTet4;
            Matrix< IndexMat >         tEdgeNodesForNormal = tTet4.get_node_map_outward_normal( mElementSideOrdinals( aPairIndex ) );

            // Get vector along these edges
            Matrix< moris::F31RMat > tEdge0Vector = moris::linalg_internal::trans( tVertexCoords.get_row( tEdgeNodesForNormal( 1, 0 ) ) - tVertexCoords.get_row( tEdgeNodesForNormal( 0, 0 ) ) );
            Matrix< moris::F31RMat > tEdge1Vector = moris::linalg_internal::trans( tVertexCoords.get_row( tEdgeNodesForNormal( 1, 1 ) ) - tVertexCoords.get_row( tEdgeNodesForNormal( 0, 1 ) ) );

            // Take the cross product to get the normal
            Matrix< F31RMat > tOutwardNormal = moris::cross( tEdge0Vector, tEdge1Vector );

            // Normalize
            Matrix< F31RMat > tUnitOutwardNormal = tOutwardNormal / moris::norm( tOutwardNormal );

            return tUnitOutwardNormal;
        }

        /*
         * Extract the interface element as a standard element rather than a double-sided side set type structure
         * If elements are
         */
        Matrix< IndexMat >
        extract_as_standard_element_loc_inds()
        {
            moris::uint        tNumVertsPerElem = mElementPairs( 0 )->get_number_of_vertices();
            Matrix< IndexMat > tElementToNode( 1, 6 );

            MORIS_ASSERT( tNumVertsPerElem == 4 && mElementPairs( 0 )->get_geometry_type() == moris::mtk::Geometry_Type::TET, "Only supporting tet4 geometry" );

            moris::uint tCount = 0;
            for ( moris::uint i = 0; i < mElementPairs.size(); i++ )
            {
                Matrix< IndexMat > tSingleElemToNode = mElementPairs( i )->get_vertex_inds();
                Matrix< IndexMat > tNodesOnFaceMap   = get_nodes_on_face_map( tNumVertsPerElem, mElementPairs( i )->get_geometry_type(), mElementSideOrdinals( i ) );
                for ( moris::uint j = 0; j < tNodesOnFaceMap.numel(); j++ )
                {
                    tElementToNode( tCount ) = tSingleElemToNode( tNodesOnFaceMap( j ) );
                    tCount++;
                }
            }

            return tElementToNode;
        };

        /*!
         * Return the extrected interface element topology
         */
        mtk::CellTopology
        get_extracted_cell_topology()
        {
            moris::uint               tNumVertsPerElem = mElementPairs( 0 )->get_number_of_vertices();
            moris::mtk::Geometry_Type tGeometryType    = mElementPairs( 0 )->get_geometry_type();

            if ( tNumVertsPerElem == 4 && tGeometryType == moris::mtk::Geometry_Type::TET )
            {
                return mtk::CellTopology::PRISM6;
            }
            else
            {
                return mtk::CellTopology::UNDEFINED;
            }
        }

      private:
        moris::moris_id                   mElementId;
        moris::moris_index                mElementIndex;
        moris::moris_index                mElementOwner;
        Vector< const moris::mtk::Cell* > mElementPairs;
        Matrix< IndexMat >                mElementSideOrdinals;
    };

}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_ */
