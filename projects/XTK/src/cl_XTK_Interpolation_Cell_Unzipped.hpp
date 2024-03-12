/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Interpolation_Cell_Unzipped.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_UNZIPPED_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_UNZIPPED_HPP_

#include "cl_XTK_Interpolation_Cell.hpp"
#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

using namespace moris;

namespace moris
{
    namespace mtk
    {
        class Vertex;
    }
}    // namespace moris

//------------------------------------------------------------------------------

namespace moris::xtk
{

    class Interpolation_Vertex_Unzipped;

    class Interpolation_Cell_Unzipped : public Interpolation_Cell
    {
        //------------------------------------------------------------------------------

      protected:
        // friend class
        friend class Enrichment;
        friend class Ghost_Stabilization;
        friend class Enriched_Interpolation_Mesh;

        moris::mtk::Cell*                             mBaseCell;
        moris::moris_index                            mSubPhaseIndex;
        moris::moris_index                            mBulkPhaseIndex;
        Vector< xtk::Interpolation_Vertex_Unzipped* > mVertices;

        Vector< moris::moris_index > mSpgIndices;
        Vector< moris::moris_index > mBulkPhaseIndices;

        //------------------------------------------------------------------------------

        Vector< xtk::Interpolation_Vertex_Unzipped* >&
        get_xtk_interpolation_vertices();

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        Interpolation_Cell_Unzipped(){};

        //------------------------------------------------------------------------------

        // constructor for SPG based initialization
        Interpolation_Cell_Unzipped(
                moris::mtk::Cell*                        aBaseCell,
                moris_index                              aSubphaseIndex,
                moris_index                              aBulkPhaseIndex,
                moris_id                                 aCellId,
                moris_index                              aCellIndex,
                moris_id                                 aCellOwner,
                std::shared_ptr< moris::mtk::Cell_Info > aConnectivity );

        //------------------------------------------------------------------------------

        // constructor for SPG based initialization
        // FIXME: this can go eventually, leave for now to distinquish constructors
        Interpolation_Cell_Unzipped(
                moris::mtk::Cell*                        aBaseCell,
                moris_index                              aPrimarySubPhaseIndex,
                moris_index                              aPrimaryBulkPhaseIndex,
                moris_index                              aEnrCellIndex,
                moris_id                                 aEnrCellOwner,
                moris_index                              aNumMeshIndices,
                std::shared_ptr< moris::mtk::Cell_Info > aConnectivity,
                bool                                     aIsSpgBasedConstruction );

        //------------------------------------------------------------------------------
        // MTK Interpolation Cell Implementation
        // see base class for documentation
        //------------------------------------------------------------------------------

        uint
        get_level() const;

        //------------------------------------------------------------------------------

        uint
        get_number_of_vertices() const;

        //------------------------------------------------------------------------------

        Vector< mtk::Vertex* >
        get_vertex_pointers() const;

        //------------------------------------------------------------------------------

        void
        remove_vertex_pointer( moris_index aIndex );

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_vertex_coords() const;

        //------------------------------------------------------------------------------

        void
        set_vertices( Vector< xtk::Interpolation_Vertex_Unzipped* > const & aVertexPointers );

        //------------------------------------------------------------------------------

        moris::mtk::Cell const *
        get_base_cell() const;

        //------------------------------------------------------------------------------

        moris::mtk::Cell*
        get_base_cell();

        //------------------------------------------------------------------------------
        // End Mtk Interpolation Cell Implementation
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Accessor functions of XTK specific data structures
        //------------------------------------------------------------------------------

        /**
         * @brief Return the processor local sub phase index of this interpolation cell
         */
        moris_index
        get_subphase_index() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Return the bulk phase index of the subphase
         */
        moris_index
        get_bulkphase_index() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Set the SPG and BP indices for a given DM list index
         *
         * @param aMeshListIndex
         * @param aSpgIndex
         * @param aBulkPhaseIndex
         */
        void
        set_SPG_and_BP_indices_for_DM_list_index(
                const moris_index aMeshListIndex,
                const moris_index aSpgIndex,
                const moris_index aBulkPhaseIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Retrieve the SPG index for the primary phase for a given DM list index from the Enr. IP cell
         *
         * @param aMeshListIndex
         * @return moris_index
         */
        moris_index
        get_SPG_index_for_DM_list_index( const moris_index aMeshListIndex ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Retrieve the bulk-phase index for the primary phase for a given DM list index from the Enr. IP cell
         *
         * @param aMeshListIndex
         * @return moris_index
         */
        moris_index
        get_bulkphase_index_for_DM_list_index( const moris_index aMeshListIndex ) const;

        //------------------------------------------------------------------------------

        Vector< xtk::Interpolation_Vertex_Unzipped* > const &
        get_xtk_interpolation_vertices() const;

        //------------------------------------------------------------------------------

        // memory
        size_t
        capacity();

        //------------------------------------------------------------------------------

        virtual const luint*
        get_ijk() const
        {
            return mBaseCell->get_ijk();
        }

        //------------------------------------------------------------------------------

    };    // class Interpolation_Cell_Unzipped

    //------------------------------------------------------------------------------
    // free functions
    //------------------------------------------------------------------------------

    /*
     * Outputting of this data structure
     */
    inline std::ostream&
    operator<<( std::ostream& os, const xtk::Interpolation_Cell_Unzipped& dt )
    {
        os << "Cell Id: " << std::right << std::setw( 9 ) << dt.get_id()
           << " | Cell Index: " << std::setw( 9 ) << dt.get_index()
           << " | Base Cell Id: " << std::setw( 9 ) << dt.get_base_cell()->get_id()
           << " | Subphase: " << std::setw( 9 ) << dt.get_subphase_index()
           << " | Bulkphase: " << std::setw( 9 ) << dt.get_bulkphase_index();

        // vertex interpolation
        std::cout << "\n Interpolation Cell Vertices:" << std::endl;
        Vector< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = dt.get_xtk_interpolation_vertices();

        for ( moris::uint i = 0; i < tVertices.size(); i++ )
        {
            os << "    " << *tVertices( i ) << std::endl;
        }

        return os;
    }

    //------------------------------------------------------------------------------

    inline std::ostream&
    operator<<( std::ostream& os, xtk::Interpolation_Cell_Unzipped const * const & dt )
    {
        os << *dt;
        return os;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERPOLATION_CELL_UNZIPPED_HPP_ */
