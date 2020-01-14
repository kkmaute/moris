/*
 * cl_Design_Variable_Interface_Proxy.hpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
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
    Cell< enum MSI::Dv_Type >      mDvTypes;
    moris::Matrix< DDRMat >        mDvValues;
    moris::Matrix< DDSMat >        mIsActiveDv;
    Cell< moris::Matrix< IdMat > > mDvIds;
    moris::Matrix< DDSMat >        mMap;

public :
    Design_Variable_Interface_Proxy()
    {
        mDvTypes.resize( 3 );     mDvTypes( 0 ) = MSI::Dv_Type::XCOORD;     mDvTypes( 1 ) = MSI::Dv_Type::YCOORD;      mDvTypes( 2 ) = MSI::Dv_Type::DENSITY;

        mDvValues.set_size( 6, 3 );
        mDvValues( 0, 0 ) = 0;      mDvValues( 0, 1 ) = 0;             mDvValues( 0, 2 ) = 0;
        mDvValues( 1, 0 ) = 1;      mDvValues( 1, 1 ) = 0;             mDvValues( 1, 2 ) = 0;
        mDvValues( 2, 0 ) = 1;      mDvValues( 2, 1 ) = 1;             mDvValues( 2, 2 ) = 0;
        mDvValues( 3, 0 ) = 0;      mDvValues( 3, 1 ) = 1;             mDvValues( 3, 2 ) = 0;
        mDvValues( 4, 0 ) = 0.5;    mDvValues( 4, 1 ) = 0;             mDvValues( 4, 2 ) = 0;
        mDvValues( 5, 0 ) = 0;      mDvValues( 5, 1 ) = 0.5;           mDvValues( 5, 2 ) = 0;

        mIsActiveDv.set_size( 6, 3 );
        mIsActiveDv( 0, 0 ) = 0;      mIsActiveDv( 0, 1 ) = 0;             mIsActiveDv( 0, 2 ) = 1;
        mIsActiveDv( 1, 0 ) = 0;      mIsActiveDv( 1, 1 ) = 0;             mIsActiveDv( 1, 2 ) = 1;
        mIsActiveDv( 2, 0 ) = 0;      mIsActiveDv( 2, 1 ) = 0;             mIsActiveDv( 2, 2 ) = 1;
        mIsActiveDv( 3, 0 ) = 0;      mIsActiveDv( 3, 1 ) = 0;             mIsActiveDv( 3, 2 ) = 1;
        mIsActiveDv( 4, 0 ) = 1;      mIsActiveDv( 4, 1 ) = 1;             mIsActiveDv( 4, 2 ) = 0;
        mIsActiveDv( 5, 0 ) = 1;      mIsActiveDv( 5, 1 ) = 1;             mIsActiveDv( 5, 2 ) = 0;

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

        // create map object
        Matrix_Vector_Factory tMatFactory( MapType::Epetra );

        mVectorMap = tMatFactory.create_map( this->get_my_local_global_map() );

        mVector = tMatFactory.create_vector( nullptr, mVectorMap, VectorType::FREE );
    }

    // ----------------------------------------------------------------------------------------------
    ~Design_Variable_Interface_Proxy(){};

    //------------------------------------------------------------------------------

    void get_dv_types_for_set( const moris::moris_index          aIntegrationMeshSetIndex,
                                     Cell< enum MSI::Dv_Type > & aDvTypes )
    {
        aDvTypes = mDvTypes;
    };

    //------------------------------------------------------------------------------

    void get_pdv_value( const moris::Cell< moris::moris_index > & aNodeIndices,
                        const Cell< enum MSI::Dv_Type >         & aDvTypes,
                              Cell<moris::Matrix< DDRMat > >          & aDvValues,
                              Cell<moris::Matrix< DDSMat > >          & aIsActiveDv )
    {
        aDvValues  .resize( aDvTypes.size() );
        aIsActiveDv.resize( aDvTypes.size() );

        for ( uint Ik = 0; Ik < aDvTypes.size(); Ik++ )
        {
            aDvValues(Ik)  .set_size( aNodeIndices.size(), 1, MORIS_REAL_MAX );
            aIsActiveDv(Ik).set_size( aNodeIndices.size(), 1, MORIS_SINT_MAX );

            for ( uint Ii = 0; Ii < aNodeIndices.size(); Ii++ )
            {
                aDvValues( Ik )( Ii ) = mDvValues( Ii, Ik );
                aIsActiveDv( Ik )( Ii ) = mIsActiveDv( Ii, Ik );
            }
        }
    }

    //------------------------------------------------------------------------------

    moris::Matrix< DDSMat > get_my_local_global_map()
    {
         return mMap;
    }

    //------------------------------------------------------------------------------

    void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                      const Cell< enum MSI::Dv_Type >         & aDvTypes,
                                            Cell<moris::Matrix< IdMat > >     & aDvIds )
    {
        aDvIds.resize( aDvTypes.size() );

        for ( uint Ik = 0; Ik < aDvTypes.size(); Ik++ )
        {
            aDvIds( Ik ).set_size( aNodeIndices.size(), 1, MORIS_UINT_MAX );

            for ( uint Ii = 0; Ii < aNodeIndices.size(); Ii++ )
            {
               aDvIds( Ik )( Ii ) = mDvIds( Ik )( aNodeIndices( Ii ) );
            }
        }
    }
};
}}
#endif /* SRC_MSI_CL_DESIGN_VARIABLE_INTERFACE_PROXY_HPP_ */
