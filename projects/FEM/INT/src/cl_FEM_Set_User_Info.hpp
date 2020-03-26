/*
 * cl_FEM_Set_User_Info.hpp
 *
 *  Created on: Oct 23, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_SET_USER_INFO_HPP_
#define SRC_FEM_CL_FEM_SET_USER_INFO_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include "cl_FEM_IWG.hpp"                 //FEM/MSI/src
#include "cl_FEM_IQI.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
    class IWG;
    class IQI;

//------------------------------------------------------------------------------
        /**
         * Set_User_Info
         */
        class Set_User_Info
        {
        protected :

        // cell of IWG pointers
        moris::Cell< std::shared_ptr< IWG > > mIWGs;

        // cell of IQIs pointers
        moris::Cell< std::shared_ptr< IQI > > mIQIs;

        // set mesh index
        uint mMeshIndex;

        // mesh set name
        std::string mMeshSetName;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Set_User_Info(){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
             ~Set_User_Info(){};

//------------------------------------------------------------------------------
             /**
              * print names
              */
             void print_names()
             {
                 // print the mesh set name
                 std::cout<<"Mesh set name: "<<mMeshSetName<<std::endl;

                 // print IWG names
                 for ( uint iIWG = 0; iIWG < mIWGs.size(); iIWG++ )
                 {
                     std::cout<<"IWG name: "<<mIWGs( iIWG )->get_name()<<std::endl;
                 }

                 // print IQI names
                 for ( uint iIQI = 0; iIQI < mIQIs.size(); iIQI++ )
                 {
                     std::cout<<"IQI name: "<<mIQIs( iIQI )->get_name()<<std::endl;
                 }
             }

//------------------------------------------------------------------------------
            /**
             * set the mesh index
             * @param[ in ] aMeshIndex mesh index for set
             */
            void set_mesh_index( uint aMeshIndex )
            {
                mMeshIndex = aMeshIndex;
            };

//------------------------------------------------------------------------------
            /**
             * get the mesh index
             * @param[ out ] mMeshIndex mesh index for set
             */
            uint get_mesh_index()
            {
                return mMeshIndex;
            };

//------------------------------------------------------------------------------
            /**
             * set the mesh set name
             * @param[ in ] aMeshSetName mesh set name
             */
            void set_mesh_set_name( std::string aMeshSetName )
            {
                mMeshSetName = aMeshSetName;
            };

//------------------------------------------------------------------------------
            /**
             * get the mesh index
             * @param[ out ] mMeshSetName mesh set name
             */
            std::string get_mesh_set_name()
            {
                return mMeshSetName;
            };

//------------------------------------------------------------------------------
            /**
             * set IWGs
             * @param[ in ] aIWGs list of IWG pointers
             */
            void set_IWGs( const moris::Cell< std::shared_ptr< fem::IWG > > & aIWGs )
            {
                mIWGs = aIWGs;
            }

//------------------------------------------------------------------------------
            /**
             * set IWG
             * @param[ in ] aIWG IWG pointer
             */
            void set_IWG( std::shared_ptr< fem::IWG > aIWG )
            {
                mIWGs.push_back( aIWG );
            }

//------------------------------------------------------------------------------
            /**
             * get IWGs
             * @param[ out ] mIWGs list of IWG pointers
             */
            const moris::Cell< std::shared_ptr< fem::IWG > > & get_IWGs() const
            {
                return mIWGs;
            }

//------------------------------------------------------------------------------
            /**
             * set IQIs
             * @param[ in ] aIQIs list of IQI pointers
             */
            void set_IQIs( const moris::Cell< std::shared_ptr< fem::IQI > > & aIQIs )
            {
                mIQIs = aIQIs;
            }

//------------------------------------------------------------------------------
            /**
             * set IQI
             * @param[ in ] aIQI IQI pointer
             */
            void set_IQI( std::shared_ptr< fem::IQI > aIQI )
            {
                mIQIs.push_back( aIQI );
            }

//------------------------------------------------------------------------------
            /**
             * get IQIs
             * @param[ out ] mIQIs list of IQI pointers
             */
            const moris::Cell< std::shared_ptr< fem::IQI > > & get_IQIs() const
            {
                return mIQIs;
              }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_SET_USER_INFO_HPP_ */
