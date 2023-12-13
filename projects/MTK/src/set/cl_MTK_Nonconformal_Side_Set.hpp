/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Nonconformal_Side_Set.hpp
 *
 */
#ifndef MORIS_CL_MTK_NONCONFORMAL_SIDE_SET_HPP
#define MORIS_CL_MTK_NONCONFORMAL_SIDE_SET_HPP


#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Nonconformal_Side_Cluster.hpp"

namespace moris::mtk
{
    class Nonconformal_Side_Set : public Double_Side_Set
    {
      public:
        Nonconformal_Side_Set(
                std::string const                                      &aName,
                moris::Cell< Cluster const * > const &aClusters,
                Matrix< IndexMat > const                               &aColors,
                uint const                                             &aSpatialDim )
                : Double_Side_Set( aName, aClusters, aColors, aSpatialDim )
        {
            mSetType = mtk::SetType::NONCONFORMAL_SIDESET;
        };
    };
}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_NONCONFORMAL_SIDE_SET_HPP
