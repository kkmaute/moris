/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mesh_Builder.hpp
 *
 */

#ifndef INCLUDE_MESH_CL_MESH_BUILDER_HPP_
#define INCLUDE_MESH_CL_MESH_BUILDER_HPP_

// Linear Algebra headers

// Mesh Headers
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Bucket.hpp"
#include "mesh/cl_Mesh_Side_Set_Input.hpp"
#include "mesh/cl_Mesh_Node_Set_Input.hpp"

// XTKL: Matrix headers
#include "linalg/cl_XTK_Matrix.hpp"
#include "containers/cl_XTK_Cell.hpp"

using namespace xtk;

namespace mesh
{
    template< typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix >
    class Mesh_Builder
    {
      public:
        virtual ~Mesh_Builder()
        {
        }

        virtual std::shared_ptr< mesh::Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > > build_mesh_from_string( std::string const &aMeshFileName,
                Vector< std::string > const                                                                                                &aScalarFieldNames,
                bool                                                                                                                        aCreateFacesAndEdges ) = 0;

        virtual std::shared_ptr< Mesh_Data< Real, Integer, Real_Matrix, Integer_Matrix > >
        build_mesh_from_data( Integer const                                                   aSpatialDimension,
                Vector< Bucket< Integer, Integer_Matrix > > const                            &aElementBuckets,
                Vector< Side_Set_Input< Integer, Integer_Matrix > > const                    &aSideSets,
                Vector< Node_Set_Input< Real, Integer, Real_Matrix, Integer_Matrix > > const &aNodeSets,
                moris::Matrix< Real_Matrix > const                                           &aNodeCoordinates,
                moris::Matrix< Integer_Matrix > const                                        &aLocaltoGlobalNodeMap,
                Vector< std::string > const                                                  &aElementPartNames,
                Vector< enum EntityTopology >                                                &aElementPartTopologys,
                Vector< std::string > const                                                  &aSideSetNames,
                Vector< std::string > const                                                  &aNodeSetNames,
                Vector< std::string > const                                                  &aInterfaceNodeSetNames,
                Vector< std::string > const                                                  &aInterfaceSideSetNames,
                Vector< std::string > const                                                  &aRealScalarFields,
                Vector< std::string > const                                                  &aRealVectorFields,
                Vector< std::string > const                                                  &aElementRealFields,
                Sensitivity< Real, Integer, Real_Matrix, Integer_Matrix >                    &aSensitivityData,
                bool const                                                                   &aSetupDataForInternalUse = false ) const = 0;
    };
}    // namespace mesh

#endif /* INCLUDE_MESH_CL_MESH_BUILDER_HPP_ */
