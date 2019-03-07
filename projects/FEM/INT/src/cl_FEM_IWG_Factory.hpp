/*
 * cl_FEM_IWG_Factory.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_FACTORY_HPP_
#define SRC_FEM_CL_FEM_IWG_FACTORY_HPP_

#include "assert.hpp"

#include "cl_Matrix.hpp"
#include "cl_FEM_IWG.hpp" //FEM/INT/src

#include "cl_FEM_IWG_L2.hpp"                     //FEM/INT/src
#include "cl_FEM_IWG_Spatial_Diffusion_Bulk.hpp" //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"         //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Bulk2.hpp"        //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Interface.hpp"    //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"   //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp"  //FEM/INT/src
#include "cl_FEM_IWG_LSNormal_Bulk.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"        //FEM/INT/src
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"   //FEM/INT/src


namespace moris
{
//------------------------------------------------------------------------------
    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief IWG factory
     */
    class IWG_Factory
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * constructor
         */
        IWG_Factory();

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~IWG_Factory();

//------------------------------------------------------------------------------
        /**
         * create IWGs
         */
        IWG * create_IWGs( IWG_Type aIWGType );

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_FACTORY_HPP_ */
