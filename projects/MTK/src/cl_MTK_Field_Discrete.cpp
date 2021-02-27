
#include "cl_MTK_Field_Discrete.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Discrete::Field_Discrete(
                mtk::Mesh_Pair * aMeshPairs,
                uint     const & aDiscretizationMeshIndex,
                uint     const & mNumberOfFields)
                : Field(aMeshPairs,mNumberOfFields),
                  mDiscretizationMeshIndex(aDiscretizationMeshIndex)
                  {
            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair->mInterpolationMesh;

            // get number of coefficients
            mNumberOfCoefficients = tIPmesh->get_max_num_coeffs_on_proc(mDiscretizationMeshIndex);

            // set coefficient vector
            mCoefficients.set_size( mNumberOfCoefficients, 1);
                  }

        // ----------------------------------------------------------------------------------------------

        uint Field_Discrete::get_discretization_order() const
        {
            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair->mInterpolationMesh;

            // return discretization order
            return tIPmesh->get_discretization_order( this->get_discretization_mesh_index() );
        }

        // ----------------------------------------------------------------------------------------------

        moris_index Field_Discrete::get_discretization_mesh_index() const
        {
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::get_discretization_mesh_index - discretization mesh index not set.\n");

            return mDiscretizationMeshIndex;
        }

        //------------------------------------------------------------------------------

        void Field_Discrete::set_coefficient_vector(const Matrix< DDRMat > & aCoefficients)
        {
            // set coefficient initialization flag to true
            mCoefficientsAreInitialized = true;

            // copy coefficients
            mCoefficients = aCoefficients;
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Discrete::compute_nodal_values()
        {
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::compute_nodal_values - discretization mesh index not set.\n");

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair->mInterpolationMesh;

            // check that number of coefficient on mesh matches size of coefficient vector
            MORIS_ASSERT( mCoefficients.n_rows() == tIPmesh->get_max_num_coeffs_on_proc(mDiscretizationMeshIndex),
                    "Field_Discrete::compute_nodal_values -  number of coefficient on mesh does not match size of coefficient vector.\n");

            // check that coefficient vector has been initialized
            MORIS_ASSERT( mCoefficientsAreInitialized,
                    "Field_Discrete::compute_nodal_values - coefficient vector has not been initialized.\n");

            // make sure that nodal value matrix is properly sized
            mNodalValues.resize( tIPmesh->get_num_nodes(), mNumberOfFields );

            //loop over all nodes
            for (uint tNodeIndex=0;tNodeIndex<mNodalValues.n_cols();++tNodeIndex)
            {
                // get t-matrix and coefficient indices
                const Matrix<IndexMat> & tCoefIndices = tIPmesh->get_coefficient_indices_of_node(tNodeIndex, mDiscretizationMeshIndex);
                const Matrix<DDRMat>   & tMatrix      = tIPmesh->get_t_matrix_of_node_loc_ind(tNodeIndex, mDiscretizationMeshIndex);

                // compute nodal value : t-matrix * coefficients
                real tValue=0.0;

                for (uint tCoef=0;tCoef<tCoefIndices.numel();++tCoef)
                {
                    tValue += tMatrix(tCoef) * mCoefficients(tCoefIndices(tCoef));
                }

                // apply nodal value to all fields
                for (uint tFieldIndex=0;tFieldIndex<mNumberOfFields;++tFieldIndex)
                {
                    mNodalValues(tNodeIndex,tFieldIndex)=tValue;
                }
            }
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Discrete::compute_derivatives_of_field_value(
                Matrix< DDRMat >       & aDerivatives,
                Matrix< IndexMat >     & aCoefIndices,
                uint             const & aNodeIndex,
                uint             const & aFieldIndex)
        {
            MORIS_ASSERT( mDiscretizationMeshIndex > -1,
                    "Field_Discrete::compute_nodal_values - discretization mesh index not set.\n");

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair->mInterpolationMesh;

            aCoefIndices = tIPmesh->get_coefficient_indices_of_node(aNodeIndex, mDiscretizationMeshIndex);
            aDerivatives = tIPmesh->get_t_matrix_of_node_loc_ind(aNodeIndex, mDiscretizationMeshIndex);
        }

        // ----------------------------------------------------------------------------------------------

    }
}

