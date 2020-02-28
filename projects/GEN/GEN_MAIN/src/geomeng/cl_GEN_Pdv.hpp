/*
 * cl_GEN_Pdv.hpp
 *
 *  Created on: Jan 14, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HPP_

#include "cl_GEN_Property.hpp"

namespace moris
{
namespace ge
{
    class GEN_Pdv
    {
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
        // pdv index
        moris_index      mIndex;

        // pdv value
        Matrix< DDRMat > mVal;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aPropertyPointera GEN property pointer
         */
        GEN_Pdv( std::shared_ptr< GEN_Property > aPropertyPointer )
        {
            // assign pdv value from the property pointer
            mVal = aPropertyPointer->val();
        };

//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aPdvVal a value for the pdv
         */
        GEN_Pdv( moris::real aPdvVal )
        {
            // assign pdv value directly
            mVal.resize( 1, 1 );
            mVal( 0, 0 ) = aPdvVal;
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~GEN_Pdv(){};

//------------------------------------------------------------------------------
        /**
         * set index
         * @param[ in ] aPdvIndex an index for the pdv
         */
        void set_index( moris_index aPdvIndex )
        {
            mIndex = aPdvIndex;
        }

//------------------------------------------------------------------------------
        /**
         * get index
         * @param[ out ] aPdvIndex an index for the pdv
         */
        moris_index get_index()
        {
            return mIndex;
        }

//------------------------------------------------------------------------------
        /**
         * get value
         * @param[ out ] mVal a value for the pdv
         */
        Matrix< DDRMat > & get_val()
        {
            return mVal;
        }

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------

}   // end ge namespace
}   // end moris namespace

#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HPP_ */
