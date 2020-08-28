/*
 * cl_FEM_IWG_Compressible_NS_Mass_Flux_Neumann.hpp
 *
 *  Created on: Aug 26, 2020
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_MASS_FLUX_NEUMANN_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_MASS_FLUX_NEUMANN_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Compressible_NS_Mass_Flux_Neumann : public IWG
        {
                //------------------------------------------------------------------------------
            public:
                enum class IWG_Property_Type
                {
                    MASS_FLUX,
                    MAX_ENUM
                };

                // Local string to property enum map
                std::map< std::string, IWG_Property_Type > mPropertyMap;

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IWG_Compressible_NS_Mass_Flux_Neumann();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Mass_Flux_Neumann(){};

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_jacobian( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_jacobian_and_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the residual wrt design variables
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dRdp( real aWStar );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_MASS_FLUX_NEUMANN_HPP_ */
