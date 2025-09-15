/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Design_Variable_Interface.hpp
 *
 */

#pragma once
#include <utility>

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Map.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MSI_QI.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Vector;
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

          protected:
            // QI Values that were requested to be used for optimization objectives or constraints
            Vector< std::string > mRequestedQIs = {};

            // All QI values that were specified in the input file, names mapped to QI objects
            std::map< std::string, MSI::QI > mQIs;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Design_Variable_Interface( Vector< std::string >& aRequestedQIs )
                    : mRequestedQIs( aRequestedQIs ) {};

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
             * Registers a single quantity of interest. Adds to list of mQIs and adds the name to mQINameToIndexMap.
             * @param[ in ] aQI Quantity of interest to register. See QI constructor for options with construction.
             */
            void register_QI( const std::string& aName, QI aQI );

            //------------------------------------------------------------------------------

            /**
             * Templated constructor to register a single quantity of interest. Forwards the arguments to the appropriate QI constructor.
             * See MSI_QI constructor for options with construction.
             *
             * @tparam Args Variadic template to accept any number of arguments of any type
             * @param[ in ] args Arguments to forward to the QI constructor
             */
            template< typename... Args >
            void register_QI( const std::string& aName, Args&&... args )
            {
                this->register_QI( aName, QI( std::forward< Args >( args )... ) );
            }

            //------------------------------------------------------------------------------


            template< typename T >
            void register_QIs( const Vector< std::string >& aQIName, const Vector< Module_Type >& aModule, Vector< T >& aQIValue )
            {
                MORIS_ASSERT( aQIName.size() == aQIValue.size(),
                        "Design_Variable_Interface::register_QI - Size of aQIName and aQIValue must be the same." );

                for ( uint iQI = 0; iQI < aQIName.size(); iQI++ )
                {
                    this->register_QI( aQIName( iQI ), QI( aModule( iQI ), aQIValue( iQI ) ) );
                }
            }

            //------------------------------------------------------------------------------

            template< typename T, typename U >
            void register_QIs( const Vector< std::string >& aQIName, const Vector< Module_Type >& aModule, Vector< T >& aQIValue, Vector< U >& adQIdADV = Vector< U >() )
            {
                MORIS_ASSERT( aQIName.size() == aQIValue.size() && aQIName.size() == adQIdADV.size(),
                        "Design_Variable_Interface::register_QI - Size of aQIName, aQIValue and adQIdADV must be the same." );

                for ( uint iQI = 0; iQI < aQIName.size(); iQI++ )
                {
                    this->register_QI( aQIName( iQI ), QI( aModule( iQI ), aQIValue( iQI ), adQIdADV( iQI ) ) );
                }
            }

            //------------------------------------------------------------------------------

            /**
             * set requested IQI type for sensitivity analysis
             * @param[ in ] aRequestedIQIType
             */
            void set_requested_QIs( const Vector< std::string >& aRequestedIQIs );

            //------------------------------------------------------------------------------

            const Vector< std::string >&
            get_all_QI_names() const;

            //------------------------------------------------------------------------------

            /**
             * Gets the names of all requested QIs that come from a specific module (e.g. FEM, GEN, etc.)
             *
             * @param[ in ] aModule Module type to filter requested QI names by
             * @return Vector of strings containing the names of the requested QIs from the specified module
             */
            const Vector< std::string >
            get_QI_names( Module_Type aModule ) const;

            //------------------------------------------------------------------------------

            /**
             * Gets the value of the QIs that were requested for optimization as a vector of scalars
             */
            const Vector< real >
            get_all_QI_values() const;

            //------------------------------------------------------------------------------

            /**
             * Gets the value of the requested QI values as a vector of matrices. Since QI values are scalars, this function should eventually be deprecated.
             */
            const Vector< Matrix< DDRMat > >
            get_all_QI_values_mat() const;

            //-----------------------------------------------------------------------------

            /**
             * Brendan documentation
             */
            real get_QI( const std::string& aQIName ) const;

            /**
             * Brendan documentation
             */
            const Matrix< DDRMat >& get_dQIdADV( const std::string& aQIName ) const;

            /*
             * Brendan documentation
             */
            void update_QI( const std::string& aQIName, real aValue );

            void update_QI( const std::string& aQIName, const Matrix< DDRMat >& adQIdADV );

            void update_QI( const std::string& aQIName, real aValue, const Matrix< DDRMat >& adQIdADV );

            /**
             * get unique dv types for set
             * @param[ in ] aIntegrationMeshSetIndex
             * @param[ in ] aDvTypes
             */
            virtual void
            get_ip_unique_dv_types_for_set(
                    const moris_index             aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) const = 0;

            virtual void get_ig_unique_dv_types_for_set(
                    const moris_index             aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ] aIntegrationMeshSetIndex  integration Mesh index
             * @param[ in ] aDvTypes                  list of group of dv types
             */
            virtual void get_ip_dv_types_for_set(
                    const moris_index                       aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const = 0;

            virtual void get_ig_dv_types_for_set(
                    const moris_index                       aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const = 0;

            //------------------------------------------------------------------------------

            virtual void set_GenMeshMap( const Vector< moris_index >& aGenMeshMap )
            {
                MORIS_ERROR( false, "MSI_Design_Variable_Interface::set_GenMeshMap() - This function is not defined in this class" );
            }

            virtual void get_ig_pdv_value(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< DDRMat > >&         aDvValues,
                    Vector< Vector< bool > >&           aIsActiveDv ) const = 0;

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
                    Vector< Matrix< DDRMat > >&         aDvValues ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * reshape pdv values
             * i.e. reshape a cell of matrix to a matrix
             * @param[ in ] aPdvValues         a cell of matrices with pdv values
             * @param[ in ] aReshapedPdvValues a matrix of pdv values
             */
            virtual void
            reshape_pdv_values(
                    const Vector< Matrix< DDRMat > >& aPdvValues,
                    Matrix< DDRMat >&                 aReshapedPdvValues ) const
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
            virtual const Matrix< DDSMat >& get_my_local_global_map() = 0;
            // BRENDAN to make this work in parallel

            //------------------------------------------------------------------------------

            /**
             * return local to global dv type map
             * @param[ in ] aVertexIndex   List of vertex indices
             * @param[ in ] aDvType        List of Dv types
             * @param[ in ] aDvIds         List of Dv Ids
             */
            // BRENDAN NEED THESE FOR XQIs TO GET SENSITIVITIE
            virtual void get_ip_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< IdMat > >&          aDvIds ) const = 0;

            virtual void get_ig_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< IdMat > >&          aDvIds ) const = 0;

            //------------------------------------------------------------------------------

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ip_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) const = 0;

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

            // /**
            //  * returns the dQIdp
            //  * @param[ out ] dQIdp matrix filled with dQIdp
            //  */
            // void set_dQIdp_dist_vect( sol::Dist_Vector* adQIdp );

            //------------------------------------------------------------------------------

        };    // class Design_Variable_Interface

        //------------------------------------------------------------------------------

    }    // namespace MSI
}    // namespace moris