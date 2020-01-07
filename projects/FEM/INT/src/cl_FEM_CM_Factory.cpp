#include "assert.hpp"
#include "cl_FEM_CM_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp" //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp" //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear_Isotropic_Pressure.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        std::shared_ptr< Constitutive_Model > CM_Factory::create_CM( fem::Constitutive_Type aConstitutiveType )
        {

            switch( aConstitutiveType )
            {
                case ( Constitutive_Type::DIFF_LIN_ISO ):
                    return std::make_shared< CM_Diffusion_Linear_Isotropic >();

                case ( Constitutive_Type::STRUC_LIN_ISO ):
                    return std::make_shared< CM_Struc_Linear_Isotropic >();

                case (Constitutive_Type::STRUC_LIN_ISO_PRESSURE):
                    return std::make_shared< CM_Struc_Linear_Isotropic_Pressure >();

                default:
                    MORIS_ERROR( false, " CM_Factory::create_CM - No constitutive type specified. " );
                    return nullptr;
            }
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
