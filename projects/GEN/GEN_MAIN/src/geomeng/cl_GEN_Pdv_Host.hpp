/*
 * cl_GEN_Pdv_Host.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_

// GEN_MAIN
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Pdv.hpp"
#include "cl_GEN_Property.hpp"

// GEN_CORE
#include "cl_GEN_Dv_Enums.hpp"

// XTK includes
#include "cl_XTK_Topology.hpp"

namespace moris
{
    namespace ge
    {

        class GEN_Pdv_Host
        {
//------------------------------------------------------------------------------
        private :
            // list of all dv types on host
            Cell< std::shared_ptr< GEN_Pdv > > mPdvList;

            // list of properties for the dv types
            Cell< std::shared_ptr< GEN_Property > > mPdvProperties;

            // vertex index which this host lives on
            moris_index mIndex;

            // map between pdv and global ID (row vector)
            Matrix< IdMat > mTypeToIDMap;

        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aNumPdvs     number of pdv type
             * @param[ in ] aVertexIndex vertex index
             */
            GEN_Pdv_Host( const uint aNumPdvs,
                          const moris_index aVertexIndex )
            : mIndex( aVertexIndex )
            {
                // set size for the pdv list
                mPdvList.resize( aNumPdvs, nullptr );

                // set size for the id map
                mTypeToIDMap.resize( aNumPdvs, 1 );
            };

//------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~GEN_Pdv_Host(){};

//------------------------------------------------------------------------------
            Matrix< IdMat > get_type_to_id_map(  )
            {
                return mTypeToIDMap;
            }

//------------------------------------------------------------------------------
            void update_pdv_list( const uint aNumNewPdvs )
            {
                uint tNumCurrPdvs = mPdvList.size();

                mPdvList.resize( tNumCurrPdvs + aNumNewPdvs, nullptr );

                mTypeToIDMap.resize( tNumCurrPdvs + aNumNewPdvs, 1 );
            }

//------------------------------------------------------------------------------
            /**
             * get index of vertex associated to pdv host
             */
            moris_index get_index()
            {
                return mIndex;
            }
//------------------------------------------------------------------------------
            /**
             * create pdv
             * @param[ in ] aFieldPointer     a field pointer for dv type
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            void create_pdv( std::shared_ptr< GEN_Field > aFieldPointer,
                             enum GEN_DV                  aPdvType,
                             const Matrix< IndexMat >   & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // if pdv not created
                if( mPdvList( tPos ) == nullptr )
                {
                    // create a pdv with field pointer
                    mPdvList( tPos ) = std::make_shared< GEN_Pdv >( aFieldPointer, mIndex );
                }
            }
//------------------------------------------------------------------------------
            /**
             * create pdv
             * @param[ in ] aPropertyPointer  a property pointer for dv type
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            void create_pdv( std::shared_ptr< GEN_Property > aPropertyPointer,
                             enum GEN_DV                     aPdvType,
                             const Matrix< IndexMat >      & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // if pdv not created
                if( mPdvList( tPos ) == nullptr )
                {
                    // create a pdv with property pointer
                    mPdvList( tPos ) = std::make_shared< GEN_Pdv >( aPropertyPointer );
                }
            }

//------------------------------------------------------------------------------
            /**
             * create pdv
             * @param[ in ] aPdvVal           a pdv value
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            void create_pdv( moris::real                aPdvVal,
                             enum GEN_DV                aPdvType,
                             const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // if pdv not created
                if( mPdvList( tPos ) == nullptr )
                {
                    // create a pdv with pdv value
                    mPdvList( tPos ) = std::make_shared< GEN_Pdv >( aPdvVal );
                }
            }

//------------------------------------------------------------------------------
            /**
             * FIXME add phase
             * get pdv by type
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            std::shared_ptr< GEN_Pdv > get_pdv_by_type(       enum GEN_DV          aPdvType,
                                                        const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // return pdf for index
                return mPdvList( tPos );
            }

//------------------------------------------------------------------------------
            /**
             * check if active per phase and type
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            sint is_active_type(       enum GEN_DV          aPdvType,
                                 const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // if pdv not created
                if( mPdvList( tPos ) == nullptr )
                {
                    // return false
                    return 0;
                }
                // if pdv created
                else
                {
                    // return true
                    return 1;
                }
            }

//------------------------------------------------------------------------------
            /**
             * assign an id to a pdv with phase and type
             * @param[ in ] aId               an id
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            void assign_id_to_type(      uint                  aId,
                                    enum GEN_DV                aPdvType,
                                    const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // if pdv created
                if( mPdvList( tPos ) != nullptr )
                {
                    // assign id
                    mTypeToIDMap( tPos ) = aId;
                }
                // if pdv not created
                else
                {
                    // assign default -1
                    mTypeToIDMap( tPos ) = gNoID;
                }
            }

//------------------------------------------------------------------------------
            /**
             * get global index for pdv with phase and type
             * @param[ in ] aPdvType          a dv type
             * @param[ in ] aGlobalPdvTypeMap a map from dv type enum to index
             */
            uint get_global_index_for_dv_type(       enum GEN_DV          aPdvType,
                                               const Matrix< IndexMat > & aGlobalPdvTypeMap )
            {
                // get index for dv type
                moris_index tPos = aGlobalPdvTypeMap( static_cast<sint>(aPdvType) );

                // return id from map
                return mTypeToIDMap( tPos );
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    }   // end ge namepsace
}       // end moris namespace

#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_ */
