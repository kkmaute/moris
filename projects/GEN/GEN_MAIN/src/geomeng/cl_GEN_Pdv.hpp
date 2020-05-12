/*
 * cl_GEN_Pdv.hpp
 *
 *  Created on: Jan 14, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_Type_HPP_
#define PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_Type_HPP_

// GEN_MAIN
#include "cl_GEN_Field.hpp"
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

        // flag to tell if this DV is changing (i.e. was it generated from a field)
        bool mIsChanging = true;

//------------------------------------------------------------------------------
    public :
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aFieldPointer a GEN Field pointer
         * @param[ in ] aEntityIndex  an index to the associated entity (so the Field returns the correct value)
         */
        GEN_Pdv( std::shared_ptr< GEN_Field > aFieldPointer,
                 moris_index                  aEntityIndex )
        {
            mVal.resize( 1, 1 );
            // assign pdv value from the property pointer
            mVal(0,0) = aFieldPointer->get_field_val_at_vertex( aEntityIndex );

            // flag this PDV_Type so the interface knows it is from a Field and therefore not changing
            mIsChanging = false;
        };
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aPropertyPointer a GEN property pointer
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
         * get value
         * @param[ out ] mVal a value for the pdv
         */
        Matrix< DDRMat > & get_val()
        {
            return mVal;
        }
//------------------------------------------------------------------------------
        /*
         * return flag for if this PDV_Type is changing or not
         */
        bool is_pdv_changing()
        {
            return mIsChanging;
        }
//------------------------------------------------------------------------------
        /*
         * set flag such that the PDV_Type is not changing
         */
        void flag_as_unchanging()
        {
            mIsChanging = false;
        }

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------

}   // end ge namespace
}   // end moris namespace

#endif /* PROJECTS_GEN_GEN_MAIN_SRC_GEOMENG_CL_GEN_PDV_HPP_ */
