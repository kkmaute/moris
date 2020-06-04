/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.hpp
 *
 *  Created on: Jan 21, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_

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

        class IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost : public IWG
        {

            protected:
                // order of Shape functions
                uint mOrder;

                //------------------------------------------------------------------------------
            public:

                enum class IWG_Constitutive_Type
                {
                    ELAST_LIN_ISO,
                    MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Constitutive_Type > mConstitutiveMap;

                enum class IWG_Stabilization_Type
                {
                        GHOST_VW_1,
                        MAX_ENUM
                };

                // Local string to constitutive enum map
                std::map< std::string, IWG_Stabilization_Type > mStabilizationMap;

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost(){};

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
                 * set stabilization parameter
                 * @param[ in ] aStabilizationParameter a stabilization parameter pointer
                 * @param[ in ] aStabilizationString    a string defining the stabilization parameter
                 */
                void set_stabilization_parameter(
                        std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                        std::string                                aStabilizationString );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_residual(  real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_jacobian( real tWStar );

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
            private:

                //------------------------------------------------------------------------------
                /**
                 * method to assemble "normal matrix" from normal vector needed for
                 * 2nd and 3rd order Ghost formulations
                 * @param[ in ] aOrderGhost Order of derivatives and ghost formulation
                 */
                void get_normal_matrix (
                        uint               aOrderGhost,
                        Matrix< DDRMat > & aNormalMatrix );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_VIRTUAL_WORK_GHOST_HPP_ */
