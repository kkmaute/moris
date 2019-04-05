/*
 * cl_GE_Interface.hpp
 *
 *  Created on: Apr 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_

// GE includes


namespace moris
{
namespace ge
{

class Geometry;
//------------------------------------------------------------
class Geometry_Engine_Interface
{
public:
    // constuctor
    Geometry_Engine_Interface(){};

    // destructor
    ~Geometry_Engine_Interface(){};
    //------------------------------------------------------------
    virtual void set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
                                                      real   aThreshold   = 0.0)
    {
        { MORIS_ERROR( false, "Geometry_Engine_Interface::set_geometry: not set."); };
    }

    //------------------------------------------------------------

private:


    //------------------------------------------------------------

protected:

    //------------------------------------------------------------
};
//------------------------------------------------------------
} // end ge namespace
} // end moris namespace

#endif /* PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_ */
