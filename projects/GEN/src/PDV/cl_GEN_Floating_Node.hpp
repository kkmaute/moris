/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Floating_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"

namespace moris::gen
{
    // Forward declare necessary classes
    class Geometry;
    class Background_Node;
    class Basis_Node;
    class Parent_Node;

    class Floating_Node : public Derived_Node
    {
      private:
        moris_id             mPDVStartingID;
        moris_id             mNodeID    = -1;
        moris_index          mNodeOwner = -1;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBackgroundNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry
         */
        Floating_Node(
                uint                              aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Matrix< DDRMat >&           aParametricCoordinates,
                mtk::Geometry_Type                aGeometryType,
                mtk::Interpolation_Order          aInterpolationOrder );

        /**
         * Gets if this node's position depends on ADVs. This means either the interface geometry or the parent nodes depend on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() const override;

        /**
         * Gets if this floating node can be determined that it is on a specific interface without any field evaluation.
         *
         * @param aGeometry Potential interface geometry
         * @return If this node is on the requested interface
         */
        bool is_on_interface( const Geometry& aGeometry ) const override;

        /**
         * Gets the number of PDVs on this floating node.
         *
         * @return Number of PDVs
         */
        uint get_num_pdvs() override;

        /**
         * Sets the starting index to be able to use the intersection coordinates of this node as PDVs
         *
         * @param aPDVStartingID The global index of the first PDV on the host
         */
        void set_starting_pdv_id( moris_id aPDVStartingID ) override;

        /**
         * Get the starting global index for the intersection coordinate PDVs
         *
         * @return The global index of the first PDV on the host
         */
        moris_id get_starting_pdv_id() override;

        /**
         * Set the node ID for this node.
         *
         * @param aNodeID Node ID
         */
        void set_id( moris_id aNodeID ) override;

        /**
         * Set the owning processor for this node.
         *
         * @param aNodeOwner Owning processor
         */
        void set_owner( moris_index aNodeOwner ) override;

        /**
         * Get the ID for this node.
         *
         * @return Node ID
         */
        moris_id get_id() const override;

        /**
         * Get the owning processor for this node.
         *
         * @return Owning processor
         */
        moris_index get_owner() override;

      protected:
        /**
         * Gets the geometry that this floating node was created on its interface.
         *
         * @return Geometry reference
         */
        virtual Geometry& get_interface_geometry() = 0;

        /**
         * Gets the geometry that this floating node was created on its interface (const version)
         *
         * @return Const geometry reference
         */
        virtual const Geometry& get_interface_geometry() const = 0;
    };
}    // namespace moris::gen
