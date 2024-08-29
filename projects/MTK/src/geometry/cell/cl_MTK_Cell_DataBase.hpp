/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_DataBase.hpp
 *
 */

#ifndef SRC_cl_MTK_Cell_DataBase_HPP_
#define SRC_cl_MTK_Cell_DataBase_HPP_

#include "cl_MTK_Cell.hpp"

namespace moris::mtk
{
    class Mesh;

    class Cell_DataBase : public mtk::Cell
    {
      private:
        mtk::Cell* mBaseCell;

        moris_index mCellIndex2;
        mtk::Mesh*  mMesh;

      public:
        //------------------------------------------------------------------------------

        /**
         * @brief Construct a new Cell_DataBase object default constructor
         *
         */

        Cell_DataBase() = default;

        //------------------------------------------------------------------------------

        /**
         * @brief Construct the main constructor for the class // used for the IP cells
         *
         * @param aCell an mtk cell which we use for id/index/owner
         * @param aCellInfo a cell info
         * @param aVertices beginning of list of vertices
         */
        Cell_DataBase( mtk::Cell&                        aCell,
                std::shared_ptr< moris::mtk::Cell_Info > aCellInfo,
                moris_index                              aCellIndex2,
                mtk::Mesh*                               aMesh );

        //------------------------------------------------------------------------------

        /**
         * @brief Construct the main constructor for the class // used for the IG cells
         *
         * @param aCell an mtk cell which we use for id/index/owner
         * @param aVertices beginning of list of vertices
         */
        Cell_DataBase( mtk::Cell& aCell,
                moris_index       aCellIndex2,
                mtk::Mesh*        aMesh );

        //------------------------------------------------------------------------------

        /**
         * @brief Construct the main constructor for the class // used for the IG cells
         *
         * @param aCell an mtk cell which we use for id/index/owner
         * @param aVertices beginning of list of vertices
         */
        Cell_DataBase( moris_index aCellIndex2,
                mtk::Mesh*         aMesh );

        //------------------------------------------------------------------------------

        /**
         * @brief Destroy the Cell_DataBase object
         *
         */

        ~Cell_DataBase() override = default;

        //------------------------------------------------------------------------------

        /**
         * @return Ptrs of vertices connected to this cell
         */

        Vector< Vertex* >
        get_vertex_pointers() const override;

        //------------------------------------------------------------------------------

        void
        remove_vertex_pointer( moris_index aIndex ) override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the vertex coords of the cell
         *
         * @return Matrix< DDRMat >  the matrix containing the cell coordinates with (NumVertices, SpatialDim )
         */
        Matrix< DDRMat >
        get_vertex_coords() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief  Returns the level that this cell is on. For most meshes this returns 0. However,
         *       for HMR this is not trivial
         *
         * @return uint
         */

        uint
        get_level() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief  Returns the level that this cell is on. For most meshes this returns 0. However,
         *       for HMR this is not trivial
         *
         * @return uint
         */
        mtk::Cell const *
        get_base_cell() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the base cell object
         *
         * @return mtk::Cell*
         */

        mtk::Cell*
        get_base_cell() override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the base cell object
         *
         * @return mtk::Cell*
         */

        const luint*
        get_ijk() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief memory usage of the object
         *
         * @return size_t
         */

        size_t
        capacity() override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the cell info object
         *
         * @return mtk::Cell_Info const*
         */

        mtk::Cell_Info const *
        get_cell_info() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the cell info sp object
         *
         * @return std::shared_ptr< mtk::Cell_Info >
         */

        std::shared_ptr< mtk::Cell_Info >
        get_cell_info_sp() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the id object
         *
         * @return moris_id
         */

        moris_id
        get_id() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the index object
         *
         * @return moris_index
         */

        moris_index
        get_index() const override;

        //------------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this cell
         * ( this information is needed for STK )
         */
        moris_id
        get_owner() const override;

        //------------------------------------------------------------------------------

    };    // class Cell_DataBase

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* cl_MTK_Cell_DataBase.hpp */
