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

namespace moris
{
namespace hmr
{
class Facet;

class Facet_Cluster
{
    bool mTrivial; /*! This indicates that the master facet contains all the information (i.e. slave = master)*/
    Facet const *              mMasterFacet; /*Facet from coarser refinement level*/
    moris::Cell<Facet const *> mSlaveFacets; /*Facet from finer refinement level*/
    // ----------------------------------------------------------------------------
public:
    // ----------------------------------------------------------------------------
    Facet_Cluster():mTrivial(false),mMasterFacet(nullptr),mSlaveFacets(0,nullptr)
    {};
    // ----------------------------------------------------------------------------
    /*!
     * Flag this facet cluster as trivial
     */
    void
    mark_trivial(){ mTrivial = true;};
    // ----------------------------------------------------------------------------
    /*!
     * Add the master facet to cluster
     */
    void
    add_master_facet(Facet const * aNewMasterFacet)
    {
        mMasterFacet = aNewMasterFacet;
    }
    // ----------------------------------------------------------------------------
    /*!
     * Add a new slave facet to the cluster
     */
    void
    add_slave_facet(Facet const * aNewSlaveFacet)
    {
        mSlaveFacets.push_back(aNewSlaveFacet);
    };
    // ----------------------------------------------------------------------------
    /*!
     * returns the master facet
     */
    Facet const *
    get_hmr_master_facet() { return mMasterFacet; };
    // ----------------------------------------------------------------------------
    /*!
     * Retrieve the hmr slave facets
     */
    moris::Cell<Facet const *> const &
    get_hmr_slave_facets(){ return mSlaveFacets; };
    // ----------------------------------------------------------------------------
};
}
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_FACET_CLUSTER_HPP_ */

