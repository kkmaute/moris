/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Design_Variable_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_MSI_DESIGN_VARIABLE_INTERFACE_HPP_
#define SRC_FEM_CL_MSI_DESIGN_VARIABLE_INTERFACE_HPP_

#include <utility>

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
    }

    namespace gen
    {
        struct Design_Extraction_Operator;
    }

    namespace MSI
    {
        class Equation_Model;
        class Design_Variable_Interface
        {
            //------------------------------------------------------------------------------

          private:
            Matrix< DDRMat >                       mTime;
            std::shared_ptr< MSI::Equation_Model > mModel         = nullptr;
            bool                                   mdQIdpImported = false;
            sol::Dist_Vector*                      mdQIdp         = nullptr;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Design_Variable_Interface() {};

            //------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~Design_Variable_Interface() {};

            //------------------------------------------------------------------------------

            /**
             * set model pointer
             * @param[ in ] aModel Model pointer
             */
            void
            set_equation_model( std::shared_ptr< MSI::Equation_Model > aModel )
            {
                mModel = std::move( aModel );
            }

            //------------------------------------------------------------------------------

            /**
             * set time
             * @param[ in ] aTime Time
             */
            void
            set_time( const Matrix< DDRMat >& aTime )
            {
                mTime = aTime;
            }

            //------------------------------------------------------------------------------

            /**
             * set requested IQI type for sensitivity analysis
             * @param[ in ] aRequestedIQIType
             */
            void set_requested_IQIs( const Vector< std::string >& aRequestedIQIs );

            //------------------------------------------------------------------------------

            /**
             * get unique dv types for set
             * @param[ in ] aIntegrationMeshSetIndex
             * @param[ in ] aDvTypes
             */
            virtual void get_ip_unique_dv_types_for_set(
                    const moris::moris_index      aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) = 0;

            virtual void get_ig_unique_dv_types_for_set(
                    const moris::moris_index      aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) = 0;

            //------------------------------------------------------------------------------

            virtual Vector< std::shared_ptr< gen::Design_Extraction_Operator > >
                    get_IG_Desgin_Extraction_Operators( Matrix< moris::IndexMat > ) = 0;

            //------------------------------------------------------------------------------

            virtual Vector< sint >
            build_local_adv_indices( Vector< std::shared_ptr< gen::Design_Extraction_Operator > >& ) = 0;

            //------------------------------------------------------------------------------

            virtual void
            populate_adv_geo_weights(
                    const std::shared_ptr< gen::Design_Extraction_Operator >& aOperator,
                    Matrix< DDRMat >&                                         aAdvGeoWeights,
                    const uint                                                aNumAdvs ) = 0;

            //------------------------------------------------------------------------------

            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ] aIntegrationMeshSetIndex  integration Mesh index
             * @param[ in ] aDvTypes                  list of group of dv types
             */
            virtual void get_ip_dv_types_for_set(
                    const moris::moris_index                aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) = 0;

            virtual void get_ig_dv_types_for_set(
                    const moris::moris_index                aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) = 0;

            //------------------------------------------------------------------------------

            virtual void set_GenMeshMap( const Vector< moris_index >& aGenMeshMap )
            {
                MORIS_ERROR( false, "MSI_Design_Variable_Interface::set_GenMeshMap() - This function is not defined in this class" );
            }

            virtual void get_ig_pdv_value(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< moris::Matrix< DDRMat > >&  aDvValues,
                    Vector< Vector< bool > >&           aIsActiveDv ) = 0;

            //------------------------------------------------------------------------------

            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             */
            virtual void get_ip_pdv_value(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< moris::Matrix< DDRMat > >&  aDvValues ) = 0;

            //------------------------------------------------------------------------------

            /**
             * reshape pdv values
             * i.e. reshape a cell of matrix to a matrix
             * @param[ in ] aPdvValues         a cell of matrices with pdv values
             * @param[ in ] aReshapedPdvValues a matrix of pdv values
             */
            virtual void
            reshape_pdv_values(
                    const Vector< moris::Matrix< DDRMat > >& aPdvValues,
                    moris::Matrix< DDRMat >&                 aReshapedPdvValues )
            {
                MORIS_ASSERT( aPdvValues.size() != 0,
                        "GEN_Design_Variable_Interface::reshape_pdv_value - pdv value vector is empty." );

                // get the number of rows and columns
                uint tRows = aPdvValues( 0 ).numel();
                uint tCols = aPdvValues.size();

                // set size for the reshaped matrix
                aReshapedPdvValues.set_size( tRows, tCols );

                // fill values into the reshaped matrix
                for ( uint iCol = 0; iCol < tCols; iCol++ )
                {
                    aReshapedPdvValues( { 0, tRows - 1 }, { iCol, iCol } ) =
                            aPdvValues( iCol ).matrix_data();
                }
            }

            //------------------------------------------------------------------------------

            /**
             * return local to global dv map
             */
            virtual const moris::Matrix< DDSMat >& get_my_local_global_map() = 0;

            //------------------------------------------------------------------------------

            /**
             * return local to global dv type map
             * @param[ in ] aVertexIndex   List of vertex indices
             * @param[ in ] aDvType        List of Dv types
             * @param[ in ] aDvIds         List of Dv Ids
             */
            virtual void get_ip_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< moris::Matrix< IdMat > >&   aDvIds ) = 0;

            virtual void get_ig_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< moris::Matrix< IdMat > >&   aDvIds ) = 0;

            //------------------------------------------------------------------------------

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ip_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) = 0;

            //------------------------------------------------------------------------------

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ig_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) = 0;

            //------------------------------------------------------------------------------

            /**
             * returns the dQIdp
             * @param[ out ] dQIdp matrix filled with dQIdp
             */
            virtual sol::Dist_Vector* get_dQIdp();

            //------------------------------------------------------------------------------

            /**
             * returns the dQIdp
             * @param[ out ] dQIdp matrix filled with dQIdp
             */
            void set_dQIdp_dist_vect( sol::Dist_Vector* adQIdp );

            //------------------------------------------------------------------------------

        };    // class Design_Variable_Interface

        //------------------------------------------------------------------------------

    }    // namespace MSI
}    // namespace moris

#endif /* SRC_FEM_CL_MSI_SOLVER_INTERFACE_HPP_ */
