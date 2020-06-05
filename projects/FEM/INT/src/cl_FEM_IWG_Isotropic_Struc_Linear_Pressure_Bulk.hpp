#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_PRESSURE_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_PRESSURE_BULK_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include <map>

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Pressure_Bulk : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                enum class IWG_Constitutive_Type
                {
                    ELAST_LIN_ISO,
                    MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IWG_Isotropic_Struc_Linear_Pressure_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Isotropic_Struc_Linear_Pressure_Bulk(){};

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model
                 * @param[ in ] aConstitutiveModel  a constitutive model pointer
                 * @param[ in ] aConstitutiveString a string defining the constitutive model
                 * @param[ in ] aIsMaster           an enum for master or slave
                 */
                void set_constitutive_model(
                        std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                        std::string                           aConstitutiveString,
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_jacobian( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
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

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_PRESSURE_BULK_HPP_ */
