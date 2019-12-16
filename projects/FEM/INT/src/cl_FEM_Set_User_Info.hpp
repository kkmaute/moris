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

        // set type
        fem::Element_Type mSetType;

        // cell of IWG pointers
        moris::Cell< std::shared_ptr< IWG > > mIWGs;

        // cell of IQIs pointers
        moris::Cell< std::shared_ptr< IQI > > mIQIs;

        // set mesh index
        uint mMeshIndex;

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
             * set set type
             * @param[ in ] aSetType an enum for the set type
             */
            void set_set_type( fem::Element_Type aSetType )
            {
                mSetType = aSetType;
            };

//------------------------------------------------------------------------------
            /**
             * return set type
             * @param[ out ] mSetType n enum for the set type
             */
            fem::Element_Type get_set_type()
            {
                return mSetType;
            };

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
             * set IWGs
             * @param[ in ] aIWGs list of IWG pointers
             */
            void set_IWGs( const moris::Cell< std::shared_ptr< fem::IWG > > & aIWGs )
            {
                mIWGs = aIWGs;
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
