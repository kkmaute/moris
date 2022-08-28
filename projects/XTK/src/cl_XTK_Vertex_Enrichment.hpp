/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Vertex_Enrichment.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_

#include "cl_Matrix.hpp"
#include "assert.h"
#include <unordered_map>
#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

using namespace moris;

namespace xtk
{
class Enriched_Interpolation_Mesh;
}

namespace xtk
{
class Vertex_Enrichment : public mtk::Vertex_Interpolation
{

protected:
    moris::moris_index               mNodeIndex;
    moris::Matrix< moris::IndexMat > mBasisIndices;
    moris::Matrix< moris::IndexMat > mBasisIds;
    moris::Matrix< moris::IndexMat > mBasisOwners;
    moris::Matrix< moris::DDRMat >   mBasisWeights;
    mtk::Vertex_Interpolation *      mBaseVertexInterp;
    std::unordered_map<moris::moris_index, moris::moris_index> mBasisMap; /*From basis to local index*/

public:

    Vertex_Enrichment();

    // Functions to interface with vertex interpolation in MTK
    //------------------------------------------------------------------------------
    /**
     * returns the IDs of the interpolation coefficients
     */
    Matrix< IdMat >
    get_ids() const;

    //------------------------------------------------------------------------------

    /**
     * returns the indices of the interpolation coefficients
     */
    Matrix< IndexMat >
    get_indices() const;

    //------------------------------------------------------------------------------

    /**
     * returns the proc owners of the IDs of this vertex
     */
    Matrix< IdMat >
    get_owners() const;

    //------------------------------------------------------------------------------
    /**
     * set the interpolation weights
     */
    void
    set_weights( const moris::Matrix< DDRMat > & aWeights )
    {
        MORIS_ERROR(0,"set_weights not implemented in xtk vertex interpolation");
    }

    //------------------------------------------------------------------------------

    /**
     * returns the interpolation weights
     */
    const Matrix< DDRMat > *
    get_weights() const;

    //------------------------------------------------------------------------------
    /**
     * set the coefficient objects
     */
    void
    set_coefficients( moris::Cell< mtk::Vertex* > & aCoefficients )
    {
        MORIS_ERROR(0,"set_coefficients not implemented in xtk vertex interpolation");
    }

    //------------------------------------------------------------------------------

    /**
     * returns the pointers to the coefficient objects
     */
    moris::Cell< mtk::Vertex* > &
    get_coefficients()
    {
        MORIS_ERROR(0,"get_coefficients not implemented in xtk vertex interpolation");
        return mCoefficients;
    }

    //------------------------------------------------------------------------------

    /**
     * returns the pointers to the coefficient objects (const version)
     */
    const moris::Cell< mtk::Vertex* > &
    get_coefficients() const
    {
        MORIS_ERROR(0,"get_coefficients not implemented in xtk vertex interpolation");
        return mCoefficients;
    }

    //------------------------------------------------------------------------------

    /**
     * returns the number of coefficients attributed to this basis
     */
    uint
    get_number_of_coefficients() const
    {
        MORIS_ERROR(0,"get_number_of_coefficients not implemented in xtk vertex interpolation");
        return 0;
    }

    void
    set_node_index(moris::moris_index aNodeIndex);

    mtk::Vertex_Interpolation const *
	get_base_vertex_interpolation() const;
    /*
     * Add the basis information which includes weights, enrichment level, and basis index.
     * There is no "smartness" in this function. Duplicates should have been removed prior to call
     * An assertion will catch duplicates in debug mode
     */
    void
    add_basis_information( moris::Matrix<moris::IndexMat> const & aBasisIndices,
                           moris::Matrix<moris::IndexMat> const & aBasisId);

    void
    add_basis_weights(moris::Matrix<moris::IndexMat> const & aBasisIndices,
                      moris::Matrix<moris::DDRMat>   const & aBasisWeight);

    void
    add_basis_owners( moris::Matrix<moris::IndexMat> const & aBasisIndices,
                      moris::Matrix<moris::IndexMat> const & aBasisOwners);

    void
    add_base_vertex_interpolation(mtk::Vertex_Interpolation * aBaseVertInterp);

    std::unordered_map<moris::moris_index, moris::moris_index> &
    get_basis_map();

    moris::uint
    local_basis_index(moris::uint aBasisIndex);

    void
    condense_out_basis_with_0_weight();

    moris::Matrix< moris::IndexMat > const &
    get_basis_indices() const;

    moris::Matrix< moris::IndexMat > const &
    get_basis_ids() const;

    moris::Matrix< moris::DDRMat > const &
    get_basis_weights() const;

    moris::Matrix< moris::DDRMat > &
    get_basis_weights();

    bool
    basis_exists_in_enrichment(moris_index aBasisIndex) const;

    friend class Enriched_Interpolation_Mesh;

    bool
    has_interpolation() const;

    size_t
    capacity();

};

/*
 * They are considered the same if they have the same basis weights and basis coefficients,
 * not necessarily in the same order
 */
inline
bool
operator==(const Vertex_Enrichment & aA,
           const Vertex_Enrichment & aB )
{
    // get basis indices of aA,aB
    moris::Matrix< moris::IndexMat > const & tBasisIndicesA = aA.get_basis_indices();
    moris::Matrix< moris::IndexMat > const & tBasisIndicesB = aB.get_basis_indices();

    // if they do not have the same number of basis they cannot be equal
    if(tBasisIndicesA.numel() != tBasisIndicesB.numel())
    {
        return false;
    }

    // iterate through basis in aA and ask aB if they exist in their basis
    for(moris::uint i = 0; i < tBasisIndicesA.numel(); i++)
    {
        if(!aB.basis_exists_in_enrichment(tBasisIndicesA(i)))
        {
            return false;
        }
    }
    // only check weights in debug mode (since this operation is performed when constructing the enriched
    // interpolation mesh only we expect the vertex enrichments being check to have the same
    // base vertex

    // get basis weights of aA,aB
//    moris::Matrix< moris::DDRMat > const & tWeightsA = aA.get_basis_weights();
//    moris::Matrix< moris::DDRMat > const & tWeightsB = aB.get_basis_weights();

    return true;
}

inline
std::ostream &
operator<<(std::ostream & os, const xtk::Vertex_Enrichment & dt)
{
    moris::Matrix< moris::IndexMat > const & tBasisIndices = dt.get_basis_indices();
    moris::Matrix< moris::IndexMat > const & tBasisOwner   = dt.get_owners();
    moris::Matrix< moris::DDRMat >   const & tBasisWeights = dt.get_basis_weights();

    // base vertex
    mtk::Vertex_Interpolation const * tBaseVertIp = dt.get_base_vertex_interpolation();
    moris::Matrix< moris::IndexMat > tBackBasisIndices = tBaseVertIp->get_indices();
    for(moris::uint i = 0; i < tBasisIndices.numel(); i++)
    {
        os<<"Basis Index: "       <<std::setw(9)<<tBasisIndices(i);
        os<<" | Basis Weight: "   <<std::setw(9)<<tBasisWeights(i);
        os<<" | Basis Owner: "    <<std::setw(9)<<tBasisOwner(i);
        os<<" | Back Basis Index:"<<std::setw(9)<<tBackBasisIndices(i)<<std::endl;
        os<<std::endl;
    }

    return os;
}

}

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_ */

