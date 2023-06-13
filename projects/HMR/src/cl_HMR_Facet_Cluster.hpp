/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Facet_Cluster.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_FACET_CLUSTER_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_FACET_CLUSTER_HPP_

#include "cl_Cell.hpp"

namespace moris::hmr
{
    class Facet;

    class Facet_Cluster
    {
        bool mTrivial; /*! This indicates that the leader facet contains all the information (i.e. follower = leader)*/
        Facet const *              mLeaderFacet; /*Facet from coarser refinement level*/
        moris::Cell<Facet const *> mFollowerFacets; /*Facet from finer refinement level*/
        // ----------------------------------------------------------------------------
    public:
        // ----------------------------------------------------------------------------
        Facet_Cluster():mTrivial(false),mLeaderFacet(nullptr),mFollowerFacets(0,nullptr)
        {};
        // ----------------------------------------------------------------------------
        /*!
         * Flag this facet cluster as trivial
         */
        void
        mark_trivial(){ mTrivial = true;};
        // ----------------------------------------------------------------------------
        /*!
         * Add the leader facet to cluster
         */
        void
        add_leader_facet(Facet const * aNewLeaderFacet)
        {
            mLeaderFacet = aNewLeaderFacet;
        }
        // ----------------------------------------------------------------------------
        /*!
         * Add a new follower facet to the cluster
         */
        void
        add_follower_facet(Facet const * aNewFollowerFacet)
        {
            mFollowerFacets.push_back(aNewFollowerFacet);
        };
        // ----------------------------------------------------------------------------
        /*!
         * returns the leader facet
         */
        Facet const *
        get_hmr_leader_facet() { return mLeaderFacet; };
        // ----------------------------------------------------------------------------
        /*!
         * Retrieve the hmr follower facets
         */
        moris::Cell<Facet const *> const &
        get_hmr_follower_facets(){ return mFollowerFacets; };
        // ----------------------------------------------------------------------------
    };
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_FACET_CLUSTER_HPP_ */

