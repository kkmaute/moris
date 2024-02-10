/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Design_Variable_Interface_Proxy.hpp
 *
 */

#ifndef SRC_MSI_CL_DESIGN_VARIABLE_INTERFACE_PROXY_HPP_
#define SRC_MSI_CL_DESIGN_VARIABLE_INTERFACE_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_MSI_Design_Variable_Interface.hpp" // COM/src

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace MSI
    {
        class Design_Variable_Interface_Proxy : public Design_Variable_Interface
        {
        private:
          Vector< Vector< enum gen::PDV_Type >>     mDvTypes;
          Vector< enum gen::PDV_Type >            mDvTypesUnique;
            moris::Matrix< DDRMat >        mDvValues;
            moris::Matrix< DDSMat >        mIsActiveDv;
            Vector< moris::Matrix< IdMat > > mDvIds;
            moris::Matrix< DDSMat >        mMap;
            moris::Matrix< DDUMat >        mConstraintDofs;
            moris::map< gen::PDV_Type, sint > mDvToIndexMap;

        public :
            Design_Variable_Interface_Proxy()
            {
                mDvTypes.resize( 2 );
                mDvTypes( 0 ).resize( 2 );     mDvTypes( 0 )( 0 ) = gen::PDV_Type::X_COORDINATE;   mDvTypes( 0 )( 1 ) = gen::PDV_Type::Y_COORDINATE;
                mDvTypes( 1 ).resize( 1 );     mDvTypes( 1 )( 0 ) = gen::PDV_Type::DENSITY;

                mDvTypesUnique.resize( 3 );
                mDvTypesUnique = { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE, gen::PDV_Type::DENSITY };

                mDvToIndexMap[ gen::PDV_Type::X_COORDINATE ]   = 0;
                mDvToIndexMap[ gen::PDV_Type::Y_COORDINATE ]   = 1;
                mDvToIndexMap[ gen::PDV_Type::DENSITY ] = 2;

                mDvValues.set_size( 6, 3 );
                mDvValues( 0, 0 ) = 0;      mDvValues( 0, 1 ) = 0;             mDvValues( 0, 2 ) = 2;
                mDvValues( 1, 0 ) = 1;      mDvValues( 1, 1 ) = 0;             mDvValues( 1, 2 ) = 3;
                mDvValues( 2, 0 ) = 1;      mDvValues( 2, 1 ) = 1;             mDvValues( 2, 2 ) = 1;
                mDvValues( 3, 0 ) = 0;      mDvValues( 3, 1 ) = 1;             mDvValues( 3, 2 ) = 5;
                mDvValues( 4, 0 ) = 0.5;    mDvValues( 4, 1 ) = 0;             mDvValues( 4, 2 ) = 0;
                mDvValues( 5, 0 ) = 0;      mDvValues( 5, 1 ) = 0.5;           mDvValues( 5, 2 ) = 0;

                mIsActiveDv.set_size( 6, 3 );
                mIsActiveDv( 0, 0 ) = 0;      mIsActiveDv( 0, 1 ) = 0;             mIsActiveDv( 0, 2 ) = 1;
                mIsActiveDv( 1, 0 ) = 0;      mIsActiveDv( 1, 1 ) = 0;             mIsActiveDv( 1, 2 ) = 1;
                mIsActiveDv( 2, 0 ) = 0;      mIsActiveDv( 2, 1 ) = 0;             mIsActiveDv( 2, 2 ) = 1;
                mIsActiveDv( 3, 0 ) = 0;      mIsActiveDv( 3, 1 ) = 0;             mIsActiveDv( 3, 2 ) = 1;
                mIsActiveDv( 4, 0 ) = 1;      mIsActiveDv( 4, 1 ) = 0;             mIsActiveDv( 4, 2 ) = 0;
                mIsActiveDv( 5, 0 ) = 0;      mIsActiveDv( 5, 1 ) = 1;             mIsActiveDv( 5, 2 ) = 0;

                mDvIds.resize( 3 );
                for ( uint Ik = 0; Ik < mDvIds.size(); Ik++ )
                {
                    mDvIds( Ik ).set_size( 6, 1 );
                }

                mDvIds( 0 )( 0 ) = 4;      mDvIds( 1 )( 0 ) = gNoID;  mDvIds( 2 )( 0 ) = gNoID;
                mDvIds( 0 )( 1 ) = 5;      mDvIds( 1 )( 1 ) = gNoID;  mDvIds( 2 )( 1 ) = gNoID;
                mDvIds( 0 )( 2 ) = 6;      mDvIds( 1 )( 2 ) = gNoID;  mDvIds( 2 )( 2 ) = gNoID;
                mDvIds( 0 )( 3 ) = 7;      mDvIds( 1 )( 3 ) = gNoID;  mDvIds( 2 )( 3 ) = gNoID;
                mDvIds( 0 )( 4 ) = gNoID;  mDvIds( 1 )( 4 ) = 0;      mDvIds( 2 )( 4 ) = 2;
                mDvIds( 0 )( 5 ) = gNoID;  mDvIds( 1 )( 5 ) = 1;      mDvIds( 2 )( 5 ) = 3;

                mMap.set_size( 8, 1 );
                mMap( 0 ) = 0;
                mMap( 1 ) = 1;
                mMap( 2 ) = 2;
                mMap( 3 ) = 3;
                mMap( 4 ) = 4;
                mMap( 5 ) = 5;
                mMap( 6 ) = 6;
                mMap( 7 ) = 7;

            }

            // ----------------------------------------------------------------------------------------------

            ~Design_Variable_Interface_Proxy(){};

            //------------------------------------------------------------------------------

            void get_ip_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type > & aDvTypes )
            {
                aDvTypes = { gen::PDV_Type::DENSITY };
            };

            //------------------------------------------------------------------------------

            void get_ig_unique_dv_types_for_set( const moris::moris_index    aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type > & aDvTypes )
            {
                aDvTypes = { gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE };
            };

            //------------------------------------------------------------------------------

            void get_ip_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type >> & aDvTypes )
            {
                aDvTypes = {{ gen::PDV_Type::DENSITY }};
            };

            //------------------------------------------------------------------------------

            void get_ig_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type >> & aDvTypes )
            {
                aDvTypes = {{ gen::PDV_Type::X_COORDINATE, gen::PDV_Type::Y_COORDINATE }};
            };

            //------------------------------------------------------------------------------

            void get_ig_pdv_value(
                    const moris::Matrix< IndexMat >& aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector<moris::Matrix< DDRMat > >&  aDvValues,
                    Vector< Vector< bool > >&            aIsActiveDv )
            {
                // Get the number of node indices requested
                uint tNumIndices = aNodeIndices.length();

                // Get the number of dv types requested
                uint tNumTypes = aDvTypes.size();

                // Set cell size
                aIsActiveDv.resize( tNumTypes );
                aDvValues.resize(tNumTypes);

                for ( uint Ik = 0; Ik < tNumTypes; Ik++ )
                {
                    // Matrix size
                    aDvValues(Ik).set_size(tNumIndices, 1, MORIS_REAL_MAX);
                    aIsActiveDv(Ik).resize(tNumIndices, true);

                    sint tIndex = mDvToIndexMap.find( aDvTypes (Ik) );

                    for ( uint Ii = 0; Ii < aNodeIndices.numel(); Ii++ )
                    {
                        // get node index
                        uint tNodeIndex = aNodeIndices( Ii );

                        if( mIsActiveDv( tNodeIndex, tIndex ) == 1 )
                        {
                            aDvValues( Ik )( Ii ) = mDvValues( tNodeIndex, tIndex );
                        }

                        aIsActiveDv( Ik )( Ii ) = mIsActiveDv( tNodeIndex, tIndex );
                    }
                }
            }

            //------------------------------------------------------------------------------

            void get_ip_pdv_value(
                    const moris::Matrix< IndexMat >  & aNodeIndices,
                    const Vector< enum gen::PDV_Type >      & aDvTypes,
                    Vector<moris::Matrix< DDRMat > >   & aDvValues )
            {
                // Get the number of node indices requested
                uint tNumIndices = aNodeIndices.length();

                // Get the number of dv types requested
                uint tNumTypes = aDvTypes.size();

                // Set cell size
                aDvValues.resize(tNumTypes);

                for ( uint Ik = 0; Ik < aDvTypes.size(); Ik++ )
                {
                    // Set matrix size
                    aDvValues(Ik).set_size( tNumIndices, 1, MORIS_REAL_MAX );

                    sint tIndex = mDvToIndexMap.find( aDvTypes (Ik));

                    for ( uint Ii = 0; Ii < aNodeIndices.numel(); Ii++ )
                    {
                        // get node index
                        uint tNodeIndex = aNodeIndices( Ii );

                        if( mIsActiveDv( tNodeIndex, tIndex ) == 1 )
                        {
                            aDvValues( Ik )( Ii ) = mDvValues( tNodeIndex, tIndex );
                        }
                    }
                }
            }

            //------------------------------------------------------------------------------

            void reshape_pdv_values(
                    const Vector< moris::Matrix< DDRMat > > & aPdvValues,
                    moris::Matrix< DDRMat >                      & aReshapedPdvValues )
            {
                MORIS_ASSERT( aPdvValues.size() != 0,
                        "GEN_Design_Variable_Interface::reshape_pdv_value - pdv value vector is empty.");

                // get the number of rows and columns
                uint tRows = aPdvValues( 0 ).numel();
                uint tCols = aPdvValues.size();

                // set size for the reshaped matrix
                aReshapedPdvValues.set_size( tRows, tCols );

                for( uint iCol = 0; iCol < tCols; iCol++ )
                {
                    aReshapedPdvValues( { 0, tRows - 1 }, { iCol, iCol } )
                = aPdvValues( iCol ).matrix_data();
                }
            }

            //------------------------------------------------------------------------------

            const moris::Matrix< DDSMat > & get_my_local_global_map()
            {
                return mMap;
            }

            //------------------------------------------------------------------------------

            void get_ip_dv_ids_for_type_and_ind( const Matrix<IndexMat> & aNodeIndices,
                    const Vector< enum gen::PDV_Type >               & aDvTypes,
                    Vector<moris::Matrix< IdMat > >     & aDvIds )
            {
                aDvIds.resize( aDvTypes.size() );

                for ( uint Ik = 0; Ik < aDvTypes.size(); Ik++ )
                {
                    aDvIds( Ik ).set_size( aNodeIndices.length(), 1, MORIS_UINT_MAX );

                    for ( uint Ii = 0; Ii < aNodeIndices.length(); Ii++ )
                    {
                        aDvIds( Ik )( Ii ) = mDvIds( Ik )( aNodeIndices( Ii ) );
                    }
                }
            }
            //------------------------------------------------------------------------------

            void get_ig_dv_ids_for_type_and_ind( const Matrix<IndexMat>               & aNodeIndices,
                    const Vector< enum gen::PDV_Type >          & aDvTypes,
                    Vector< moris::Matrix< IdMat > > & aDvIds )
            {
                aDvIds.resize( aDvTypes.size() );

                for ( uint Ik = 0; Ik < aDvTypes.size(); Ik++ )
                {
                    aDvIds( Ik ).set_size( aNodeIndices.length(), 1, MORIS_UINT_MAX );

                    for ( uint Ii = 0; Ii < aNodeIndices.length(); Ii++ )
                    {
                        aDvIds( Ik )( Ii ) = mDvIds( Ik )( aNodeIndices( Ii ) );
                    }
                }
            }

            //------------------------------------------------------------------------------

            void get_ip_requested_dv_types( Vector< enum gen::PDV_Type > & aDvTypes )
            {
                aDvTypes = {{ gen::PDV_Type::DENSITY }};
            }

            //------------------------------------------------------------------------------

            void get_ig_requested_dv_types( Vector< enum gen::PDV_Type > & aDvTypes )
            {
                MORIS_ERROR( false, "Design_Variable_Interface_Proxy::get_ig_requested_dv_types() - not implemented in the child class" );
            }

            //------------------------------------------------------------------------------

            void set_requested_IQIs( const Vector< std::string> & aRequestedIQINames );
        };
    }
}
#endif /* SRC_MSI_CL_DESIGN_VARIABLE_INTERFACE_PROXY_HPP_ */

