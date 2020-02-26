/*
 * cl_GEN_Pdv_Host.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_

#include "cl_GEN_Pdv.hpp"
#include "cl_GEN_Property.hpp"

#include "cl_GEN_Dv_Enums.hpp"

// XTK includes
#include "cl_XTK_Topology.hpp"

namespace moris
{
    namespace ge
    {

        class GEN_Pdv_Host
        {
        private :
            Cell< std::shared_ptr< GEN_Pdv > > mPdvList;                // list of all dv types on host

            Cell< std::shared_ptr< GEN_Property > > mPdvProperties;     // list of properties for the dv types

            moris_index mIndex;                                         // vertex index which this host lives on

            Matrix< IdMat > mTypeToIDMap;                               // map between dv type and global ID (row vector)
            //------------------------------------------------------------------------------
        public:
            //------------------------------------------------------------------------------

            GEN_Pdv_Host( const uint aNumPdvs,
                          const moris_index aVertexIndex ) : mIndex( aVertexIndex )
            {
                mPdvList.resize( aNumPdvs, nullptr );

                mTypeToIDMap.resize( aNumPdvs, 1 );
            };

            //------------------------------------------------------------------------------
            ~GEN_Pdv_Host(){};
            //------------------------------------------------------------------------------
            void update_pdv_list( const uint aNumNewPdvs )
            {
                uint tNumCurrPdvs = mPdvList.size();

                mPdvList.resize( tNumCurrPdvs + aNumNewPdvs, nullptr );

                mTypeToIDMap.resize( tNumCurrPdvs + aNumNewPdvs, 1 );
            }
            //------------------------------------------------------------------------------
            moris_index get_index()
            {
                return mIndex;
            }
            //------------------------------------------------------------------------------
            void create_pdv( std::shared_ptr< GEN_Property > aPropertyPointer,
                             enum GEN_DV                     aPdvType,
                             const Matrix< IndexMat >      & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                if( mPdvList( tPos ) == nullptr )
                {
                    mPdvList( tPos ) = std::make_shared< GEN_Pdv >( aPropertyPointer );
                }
            }
            //------------------------------------------------------------------------------
            void create_pdv( moris::real              aPdvVal,
                             enum GEN_DV              aPdvType,
                             const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                if( mPdvList( tPos ) == nullptr )
                {
                    mPdvList( tPos ) = std::make_shared< GEN_Pdv >( aPdvVal );
                }
            }
            //------------------------------------------------------------------------------
            std::shared_ptr< GEN_Pdv > get_pdv_by_type(       enum GEN_DV          aPdvType,
                                                        const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                return mPdvList( tPos );
            }
            //------------------------------------------------------------------------------
            sint is_active_type(       enum GEN_DV          aPdvType,
                                 const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                if( mPdvList( tPos ) == nullptr )
                {
                    return 0;
                }
                else
                {
                    return 1;
                }
            }
            //------------------------------------------------------------------------------
            void assign_id_to_type( uint aId, enum GEN_DV aPdvType, const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                if( mPdvList( tPos ) != nullptr )
                {
                    mTypeToIDMap( tPos ) = aId;
                }
                else
                {
                    mTypeToIDMap( tPos ) = gNoID;
                }
            }
            //------------------------------------------------------------------------------
            uint get_global_index_for_dv_type( enum GEN_DV aPdvType, const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                return mTypeToIDMap( tPos );
            }

        };

    }   // end ge namepsace
}       // end moris namespace

#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_ */
