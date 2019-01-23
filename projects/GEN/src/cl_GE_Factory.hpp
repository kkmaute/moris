/*
 * cl_GE_Factory.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */
#ifndef SRC_DISTLINALG_CL_GE_FACTORY_HPP_
#define SRC_DISTLINALG_CL_GE_FACTORY_HPP_

#include <memory>
#include <iostream>
#include "cl_GE_Enums.hpp"
#include "cl_GE_Geometry.hpp"
#include "cl_GE_SDF.hpp"
#include "cl_GE_LS.hpp"
#include "assert.hpp"
#include "cl_GE_Analytical.hpp"

namespace moris
{
    namespace ge
    {

    class Ge_Factory
    {
    	private:
    		int flagType = 0;

    	protected:

    	public:
    		Ge_Factory(){
    			std::cout<<"factory constructor"<<std::endl;
    		};
    		~Ge_Factory(){
    		};

    	    /**
    	     * @brief factory member function building GE types
    	     *
    	     * @param[in] aFlagType            Dermines the method to be used for flagging.
    	     * @param[out] tFlagPointer		   GE pointer to base class.
    	     *
    	     */
    		Geometry* pick_flag(enum flagType aFlagType = flagType::LS)
    		{
    			Geometry* tFlagPointer;
    			switch(aFlagType)
    			{
    			case(flagType::SDF):
    					tFlagPointer = new SDF();
    					break;
    			case(flagType::LS):
    					tFlagPointer = new LS();
    					break;
    			case(flagType::Analytical):
    					tFlagPointer = new Analytical();
    					break;
    			default:
    					MORIS_ERROR(false, "Ge_Factory::pickFlag() please input a valid flag method");
    					break;
    			}
    		return tFlagPointer;
    		}

    };
    } /* namespace gen */
} /* namespace moris */

#endif /* SRC_DISTLINALG_CL_GE_FACTORY_HPP_ */
