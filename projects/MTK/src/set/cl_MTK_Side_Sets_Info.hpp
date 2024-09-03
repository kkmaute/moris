/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Sets_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::mtk
{
    /////////////////////////
    // STRUC FOR SIDE SET  //
    /////////////////////////
    /*
     * To declare a side set in the mesh the following information
     * is needed:
     *  mElemIdsAndSideOrds - Element Id and the side ordinal of this element
     *                         col 0 - Element Id
     *                         col 1 - Side ordinal
     *
     *  mSideSetName        - Name of the side set
     */

    struct MtkSideSetInfo
    {
        Matrix< IdMat >  *mElemIdsAndSideOrds;
        std::string       mSideSetName;
        bool              mParallelConsistencyReq = true;
        enum CellTopology mSideTopology;

        MtkSideSetInfo()
                : mElemIdsAndSideOrds()
                , mSideSetName()
                , mSideTopology( CellTopology::UNDEFINED )
        {
        }

        bool
        sideset_has_name()
        {
            return !mSideSetName.empty();
        }
    };

    inline std::ostream &
    operator<<( std::ostream &os, mtk::MtkSideSetInfo const *const &dt )
    {
        os << "Side Set Name: " << dt->mSideSetName << " | Number of Sides: " << dt->mElemIdsAndSideOrds->n_rows() << " | Side Topology: " << get_enum_str( dt->mSideTopology ) << "  | Parallel Consistent: " << dt->mParallelConsistencyReq;

        for ( moris::uint i = 0; i < dt->mElemIdsAndSideOrds->n_rows(); i++ )
        {
            os << "       Cell Id: " << std::setw( 6 ) << ( *dt->mElemIdsAndSideOrds )( i, 0 ) << " | Side Ord: " << std::setw( 6 ) << ( *dt->mElemIdsAndSideOrds )( i, 1 ) << "\n";
        }
        return os;
    }
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_SIDE_SETS_INFO_HPP_ */
