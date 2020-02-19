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
    private:
        moris_index      mIndex;

        Matrix< DDRMat > mVal;
        //------------------------------------------------------------------------------
    public :
        GEN_Pdv( std::shared_ptr< GEN_Property > aPropertyPointer )
        {
            // assign pdv value from the property pointer
            mVal = aPropertyPointer->val();
        };
        //------------------------------------------------------------------------------
        GEN_Pdv( moris::real aPdvVal )
        {
            // assign pdv value directly
            mVal.resize( 1, 1 );
            mVal(0,0) = aPdvVal;
        };
        //------------------------------------------------------------------------------
        ~GEN_Pdv(){};
        //------------------------------------------------------------------------------
        void set_index( moris_index aPdvIndex )
        {
            mIndex = aPdvIndex;
        }
        //------------------------------------------------------------------------------
        moris_index get_index()
        {
            return mIndex;
        }
        //------------------------------------------------------------------------------
        Matrix< DDRMat > get_val()
        {
            return mVal;
        }
        //------------------------------------------------------------------------------
    };
}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HPP_ */
