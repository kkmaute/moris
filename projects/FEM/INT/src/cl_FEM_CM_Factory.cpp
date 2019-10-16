#include "assert.hpp"
#include "cl_FEM_CM_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp" //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Constitutive_Model* CM_Factory::create_CM( Constitutive_Type aConstitutiveType )
        {

            switch( aConstitutiveType )
            {
                case ( Constitutive_Type::DIFF_LIN_ISO ):
                    return new CM_Diffusion_Linear_Isotropic();
                case ( Constitutive_Type::STRUC_LIN_ISO ):
                    return new CM_Struc_Linear_Isotropic();

                default:
                    MORIS_ERROR( false, " CM_Factory::create_CM - No constitutive type specified. " );
                    return nullptr;
            }
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
