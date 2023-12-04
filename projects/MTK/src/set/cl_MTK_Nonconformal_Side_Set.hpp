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
#include "cl_MTK_Side_Set.hpp"
#include "MTK/src/set/mapper/cl_MTK_PointMapper.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Nonconformal_Side_Cluster.hpp"

namespace moris::mtk
{
    class Nonconformal_Side_Set : public Double_Side_Set
    {
      private:
        std::shared_ptr< Integration_Rule > mIntegrationRule;

      public:
        Nonconformal_Side_Set(
                std::string const                  &aName,
                moris::Cell< Cluster const * >      aClusters,
                Matrix< IndexMat > const           &aColors,
                uint const                         &aSpatialDim,
                std::shared_ptr< Integration_Rule > aIntegrationRule )
                : Double_Side_Set( aName, aClusters, aColors, aSpatialDim )
                , mIntegrationRule( aIntegrationRule ){};
    };
}    // namespace moris::mtk


#endif    // MORIS_CL_MTK_NONCONFORMAL_SIDE_SET_HPP
