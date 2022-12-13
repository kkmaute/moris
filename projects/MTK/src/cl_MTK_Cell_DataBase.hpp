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
         * @brief Construct a new Cell_DataBase object defualt constrcutor
         *
         */

        Cell_DataBase() = default;

        //------------------------------------------------------------------------------

        /**
         * @brief Construct the main constrcuctor for the class // used for the IP cells
         *
         * @param aCell an mtk cell which we use for id/index/owner
         * @param aCellInfo a cell info
         * @param aVerteices begining of list of vertices
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
         * @param aVerteices begining of list of vertices
         */
        Cell_DataBase( mtk::Cell& aCell,
                moris_index       aCellIndex2,
                mtk::Mesh*        aMesh );

        //------------------------------------------------------------------------------

        /**
         * @brief Construct the main constructor for the class // used for the IG cells
         *
         * @param aCell an mtk cell which we use for id/index/owner
         * @param aVerteices begining of list of vertices
         */
        Cell_DataBase( moris_index aCellIndex2,
                mtk::Mesh*         aMesh );

        //------------------------------------------------------------------------------

        /**
         * @brief Destroy the Cell_DataBase object
         *
         */

        ~Cell_DataBase() = default;

        //------------------------------------------------------------------------------

        /**
         * @return Ptrs of vertices connected to this cell
         */

        virtual moris::Cell< Vertex* >
        get_vertex_pointers() const override;

        //------------------------------------------------------------------------------

        virtual void
        remove_vertex_pointer( moris_index aIndex ) override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the vertex coords of the cell
         *
         * @return Matrix< DDRMat >  the matrix containing the cell coordinates with (NumVertices, SpatialDim )
         */
        virtual Matrix< DDRMat >
        get_vertex_coords() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief  Returns the level that this cell is on. For most meshes this returns 0. However,
         *       for HMR this is not trivial
         *
         * @return uint
         */

        virtual uint
        get_level() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief  Returns the level that this cell is on. For most meshes this returns 0. However,
         *       for HMR this is not trivial
         *
         * @return uint
         */
        virtual mtk::Cell const *
        get_base_cell() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the base cell object
         *
         * @return mtk::Cell*
         */

        virtual mtk::Cell*
        get_base_cell() override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the base cell object
         *
         * @return mtk::Cell*
         */

        virtual const luint*
        get_ijk() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief memory usage of the object
         *
         * @return size_t
         */

        size_t
        capacity();

        //------------------------------------------------------------------------------

        /**
         * @brief Get the cell info object
         *
         * @return mtk::Cell_Info const*
         */

        virtual mtk::Cell_Info const *
        get_cell_info() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the cell info sp object
         *
         * @return std::shared_ptr< mtk::Cell_Info >
         */

        virtual std::shared_ptr< mtk::Cell_Info >
        get_cell_info_sp() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the id object
         *
         * @return moris_id
         */

        virtual moris_id
        get_id() const override;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the index object
         *
         * @return moris_index
         */

        virtual moris_index
        get_index() const;

        //------------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this cell
         * ( this information is needed for STK )
         */
        virtual moris_id
        get_owner() const;

        //------------------------------------------------------------------------------

    }; // class Cell_DataBase

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* cl_MTK_Cell_DataBase.hpp */
