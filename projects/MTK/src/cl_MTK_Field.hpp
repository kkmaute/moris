/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_MTK_FIELD_HPP_
#define PROJECTS_HMR_SRC_CL_MTK_FIELD_HPP_

#include <memory>

#include "typedefs.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        class Mesh;
        class Mesh_Manager;

        //------------------------------------------------------------------------------

        /**
         * Base class of interpolation mesh based nodaly discretized scalar or vector field; it is assumed that the values
         * at a node are a function of space and some coefficients; the base implementation provides access
         * to the nodal values, allows to set them for projection on to the coefficients; and access to the
         * coefficients. The coefficients have an ID and ownership.
         *
         * The base class requires implementations for how the nodal values and their derivatives with
         * respect to the coefficients are computed.
         */
        class Field
        {
            protected:

                //! Name of mtk::field
                std::string mLabel;

                //! Index of mtk::field
                uint mFieldIndex = MORIS_UINT_MAX;

                //! Mesh pair
                Mesh_Pair mMeshPair;

                //! Mesh pair label
                std::string mMeshPairLabel = "";

                //! Number of nodal fields
                uint mNumberOfFields = 1;

                //! Number of coefficients
                sint mNumberOfCoefficients = -1;

                //! Nodal field matrix: number of nodes x number of nodal fields
                Matrix< DDRMat > mValues;

                //! Coefficients vector: number of coefficients x 1
                Matrix< DDRMat > mCoefficients;

                //! Map from field indices to mesh coefficient IDs and owning processor rank
                Matrix < IdMat > mFieldIndexToMeshCoefficientIdAndOwnerMap;

                //! Flag that nodal values need to be updated
                bool mUpdateNodalValues = true;

                //! Lock flag
                bool mFieldIsLocked = true;

                //! Lock flag
                bool mFieldIsDiscrete = false;

                enum Field_Entity_Type mFieldEntityType = Field_Entity_Type::UNDEFINED;

                enum Field_Implementation mFieldImplementation = Field_Implementation::MTK;

                //------------------------------------------------------------------------------
            protected:
                //------------------------------------------------------------------------------
public:
                /**
                 *  @brief required implementation for computing and storing nodal values
                 *         for current coefficients
                 */
                virtual void compute_nodal_values()
                {
                    MORIS_ERROR(false,"mtk::Field::compute_nodal_values - function not implemented.\n");
                }
protected:
                //------------------------------------------------------------------------------

                /**
                 *  @brief required implementation  for computing derivatives nodal values for all
                 *         coefficients this value depends on
                 *
                 *  @param[in] aNodeIndex - index of node
                 *  @param[in] aFieldIndex - index of field; default 0
                 *
                 *  @param[out] aDerivatives - vector of derivatives
                 *  @param[out] aCoefIndices - vector coefficient indices
                 */
                virtual void compute_derivatives_of_field_value(
                        Matrix< DDRMat >       & aDerivatives,
                        Matrix< IndexMat >     & aCoefIndices,
                        uint             const & aNodeIndex,
                        uint             const & aFieldIndex)
                {
                    MORIS_ERROR(false,"mtk::Field::compute_derivatives_of_field_value - function not implemented.\n");
                }
                //------------------------------------------------------------------------------

                /**
                 * @brief set vector of nodal values
                 */
                virtual void set_nodal_value_vector( const Matrix< DDRMat > & aNodalValues );

                //------------------------------------------------------------------------------

                /**
                 * @brief fills coefficient vector
                 *
                 * @param[in]  vector of coefficients
                 */
                virtual void set_coefficient_vector(const Matrix< DDRMat > & aCoefficients);

                //------------------------------------------------------------------------------

                /**
                  * @brief updates coefficient vector
                  */
                virtual void get_coefficient_vector()
                {
                }

                //------------------------------------------------------------------------------

                /**
                  * @brief get vector of IDs and owner rank for all coefficient used by field
                  *
                  * @ return matrix of IDs and owners: number of coefficients x 2
                  */
                virtual const Matrix<IdMat> & get_coefficient_id_and_owner_vector()
                {
                    MORIS_ERROR(false,"mtk::Field::get_coefficient_id_vector - function not implemented.\n");

                   return mFieldIndexToMeshCoefficientIdAndOwnerMap;
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief updates internal data associated with coefficients
                 */

                virtual void update_coefficient_data()
                {
                }

                //------------------------------------------------------------------------------

                /**
                 * @ brief determines whether nodal value vector needs to be updated
                 */
                bool nodal_values_need_update();

                //------------------------------------------------------------------------------
            public :
                //------------------------------------------------------------------------------

                //------------------------------------------------------------------------------

                /**
                 *  @brief Field constructor
                 *
                 * @param[in]   aMeshPair - pointer to mesh pair
                 * @param[in]   aNumDim - dimension of nodal field; default 1
                 *
                 */
                Field(
                        Mesh_Pair        aMeshPair,
                        uint     const & aNumberOfFields = 1,
                        enum Field_Entity_Type aFieldEntityType = Field_Entity_Type::NODAL );

                Field( enum Field_Entity_Type aFieldEntityType = Field_Entity_Type::NODAL);

                //------------------------------------------------------------------------------

                virtual ~Field();

                //------------------------------------------------------------------------------

                /**
                 *  @brief returns the field entity type. can be NODAL or ELEMENTAL
                 *
                 * @return field entity type
                 */
                enum Field_Entity_Type get_field_entity_type();

                //------------------------------------------------------------------------------

                /**
                 *  @brief returns mesh pair
                 *
                 * @return mesh pair pointer
                 */
                Mesh_Pair get_mesh_pair();

                //------------------------------------------------------------------------------

                /**
                 *  @brief returns mesh pair label
                 *
                 * @return mesh pair label
                 */
                const std::string & get_mesh_pair_label();

                //------------------------------------------------------------------------------

                /**
                 *  @brief set mesh pair
                 *
                 * @param[in] aMeshPair - mesh pair pointer
                 */
                void set_mesh_pair( Mesh_Pair aMeshPair);

                //------------------------------------------------------------------------------

                /**
                 * @brief returns number of nodal fields; e.g., 1 for scalar fields
                 *
                 * @return number of nodal fields
                 */
                uint get_number_of_fields() const
                {
                    return mNumberOfFields;
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief triggers update of field values and coeffs
                 *
                 */

                void update_field();

                //------------------------------------------------------------------------------

                /**
                 * @brief access to number of nodes
                 *
                 * @return number of nodes
                 */
                uint get_number_of_nodes() const
                {
                    return mValues.n_rows();
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief access to number of nodes
                 *
                 * @return number of nodes
                 */
                uint get_number_of_coefficients() const
                {
                    MORIS_ERROR( mNumberOfCoefficients > -1,
                            "mtk::Field::get_number_of_coefficients - coefficient vector not set.\n");

                    return mNumberOfCoefficients;
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief returns all values for current coefficients
                 *
                 * @return matrix of value
                 */
                const
                Matrix< DDRMat > & get_values();

                //------------------------------------------------------------------------------

                /**
                 * @brief returns all maximal value
                 *
                 * @return maximal value
                 */
                moris::real get_max_value();

                //------------------------------------------------------------------------------

                /**
                 * @brief returns value of the field entity; if nodal and value is not updated all nodal values
                 *        will be computed first
                 *
                 *        Note: function will be removed soon as not consistent with child implementation
                 *
                 * @param[in]  aIndex - node index
                 * @param[in]  aFieldIndex - field index
                 *
                 * @return nodal value
                 */
                moris::real get_value(
                        const uint & aIndex,
                        const uint & aFieldIndex = 0);

                //------------------------------------------------------------------------------

                /**
                 * @brief returns value of entity; if nodal and value is not updated all nodal values
                 *        will be computed first
                 *
                 *        Note: function will be removed soon as not consistent with child implementation
                 *
                 * @param[in]  aIndices - vector of indices
                 * @param[in]  aFieldIndices - field indices
                 *
                 * @param[out] aValues - values
                 */
                void get_value(
                        Matrix< IndexMat > const & aIndex,
                        Matrix< DDRMat >         & aValues,
                        Matrix< IndexMat > const & aFieldIndex = 0);

                //------------------------------------------------------------------------------

                /**
                 *  @brief return derivatives of  nodal values for all
                 *         coefficients this value depends on
                 *
                 *  @param[in] aNodeIndex - index of node
                 *  @param[in] aFieldIndex - index of field; default 0
                 *
                 *  @param[out] aDerivatives - vector of derivatives
                 *  @param[out] aCoefIndices - vector coefficient indices
                 */
                void get_derivatives_of_field_value(
                        Matrix< DDRMat >       & aDerivatives,
                        Matrix< IndexMat >     & aCoefIndices,
                        uint             const & aNodeIndex,
                        uint             const & aFieldIndex)
                {
                    // call to child implementation
                    this->compute_derivatives_of_field_value(
                            aDerivatives,
                            aCoefIndices,
                            aNodeIndex,
                            aFieldIndex);
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief sets all nodal values
                 *
                 * @param[in]  aNodalField - matrix of size number of nodes x field dimension
                 *
                 */
                void set_values( const Matrix< DDRMat > & aNodalValues );

                //------------------------------------------------------------------------------

                /**
                 * @brief returns all coefficients
                 *
                 * @return  coefficient vector
                 *
                 */
                const Matrix< DDRMat > & get_coefficients();

                //------------------------------------------------------------------------------

                /**
                 * @brief sets all coefficients; input vector needs to have same dimension of coefficient
                 *        vector; needs to be a column vector
                 *
                 * @param[in]  aaCoefficients - vector of size number of coefficient x 1
                 *
                 */
                void set_coefficients( const Matrix< DDRMat > & aCoefficients );

                //------------------------------------------------------------------------------

                /**
                  * @brief get vector of IDs and owner rank for all coefficient used by field
                  *
                  * @ return matrix of IDs and owners: number of coefficients x 2
                  */
                const Matrix<IdMat> & get_coefficient_ids_and_owners();

                //------------------------------------------------------------------------------

                /**
                 * @brief returns name of field
                 */
                const std::string & get_label() const
                {
                    return mLabel;
                };

                //------------------------------------------------------------------------------

                /**
                 * @brief sets name of field
                 *
                 * @param[in] aLabel - name of field
                 */
                void set_label( const std::string & aLabel )
                {
                    mLabel = aLabel;
                };

                //------------------------------------------------------------------------------

                /**
                 * @brief indicates if this is a discrete field. This can be used to determine if a field has to be mapped
                 *
                 */
                bool get_field_is_discrete()
                {
                    return mFieldIsDiscrete;
                };

                //------------------------------------------------------------------------------

                /**
                 * @brief returns coordinate of a node in the field
                 *
                 * @param[in]  aNodeIndex - node index
                 *
                 */
                Matrix< DDRMat > get_node_coordinate( const moris_index & aNodeIndex ) const
                {
                    return mMeshPair.get_interpolation_mesh()->get_node_coordinate( aNodeIndex );
                }

                //------------------------------------------------------------------------------

                /**
                 * @brief returns the interpolation order of the Lagrange Mesh
                 */
                uint get_lagrange_order() const
                {
                    return mMeshPair.get_interpolation_mesh()->get_order();
                }

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief returns order of underlying discretization; unless child class implementation is provided
                 *        interpolation order of the Lagrange mesh, it is assumed that discretization oder is equal to
                 *        order of Lagrange mesh
                 */
                virtual uint get_discretization_order() const
                {
                    return this->get_lagrange_order();
                }

                // ----------------------------------------------------------------------------------------------

                /**
                 * @brief returns discretization mesh index of underlying discretization; unless child class implementation
                 *        is provided, it is assumed that discretization mesh index is 0
                 */
                virtual moris_index get_discretization_mesh_index() const
                {
                    return 0;
                }

                //------------------------------------------------------------------------------

                /**
                 * Expert user function. Allows for the underlying mesh as well as the nodal and coefficients values to be changed.
                 * Such a change can result in an unwanted behavior.
                 */
                void unlock_field()
                {
                    mFieldIsLocked = false;
                };

                //------------------------------------------------------------------------------

                /**
                 * @brief check whether field is locked and if yes, throw error
                 */
                void error_if_locked() const;

                //------------------------------------------------------------------------------

                /**
                 * return the mtk::field index
                 */
                uint get_field_index() const
                {
                    return mFieldIndex;
                }

                //------------------------------------------------------------------------------

                /**
                 *  @brief save nodal values to hdf5 file
                 *
                * @param[in]  aFilePath - name of hdf5 file
                * @param[in]  aCreateNewFile - flag for removing old (if exists) and creating new file
                 */
                void save_nodal_values_to_hdf5(
                        const std::string & aFilePath,
                        const bool aCreateNewFile=true );

                //------------------------------------------------------------------------------

                /**
                 *  @brief save coefficients to hdf5 file
                 *
                * @param[in]  aFilePath - name of hdf5 file
                * @param[in]  aCreateNewFile - flag for removing old (if exists) and creating new file
                 */
                void save_coefficients_to_hdf5(
                        const std::string & aFilePath,
                        const bool aCreateNewFile=true );

                //------------------------------------------------------------------------------

                /**
                 *  @brief load coefficients from hdf5 file; coefficients are read from file and nodal
                 *         values are computed
                 *
                 * @param[in]  aFilePath - name of hdf5 file
                 */
                void load_coefficients_from_hdf5(
                        const std::string & aFilePath );

                //------------------------------------------------------------------------------

                /**
                 *  @brief load nodal value from hdf5 file
                 *
                 * @param[in]  aFilePath - name of hdf5 file
                 */
                void load_nodal_values_from_hdf5(
                        const std::string & aFilePath );

                //------------------------------------------------------------------------------

                /**
                 *  @brief save coefficients to binary file; will always overwrite existing file
                 *
                * @param[in]  aFilePath - name of binary file
                 */
                void save_coefficients_to_binary( const std::string & aFilePath );

                //------------------------------------------------------------------------------

                /**
                 *  @brief save nodal values to binary file; will always overwrite existing file
                 *
                * @param[in]  aFilePath - name of binary file
                 */
                void save_nodal_values_to_binary( const std::string & aFilePath );

                //------------------------------------------------------------------------------

                /**
                 *  @brief save coefficients to binary file; will always overwrite existing file
                 *
                * @param[in]  aFilePath - name of hdf5 file
                 */
                void load_field_from_exodus(
                        const std::string    & aFileName,
                        const moris_index      aTimeIndex = 0,
                        const Matrix<DDUMat> & aFiledIndices = {{0}} );

                //------------------------------------------------------------------------------

                /**
                 *  @brief save coefficients to binary file; will always overwrite existing file
                 *
                * @param[in]  aFilePath - name of hdf5 file
                 */
                void save_field_to_exodus( const std::string & aFileName );

                enum Field_Implementation get_field_implementation()
                {
                    return mFieldImplementation;
                };

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_MRK_FIELD_HPP_ */

