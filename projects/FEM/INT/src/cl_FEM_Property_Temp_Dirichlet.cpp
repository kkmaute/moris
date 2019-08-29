#include "cl_FEM_Property_Temp_Dirichlet.hpp"

#include "op_times.hpp" //LINALG/src
#include "fn_norm.hpp"  //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_dot.hpp"   //LINALG/src
#include "fn_print.hpp" //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    Property_Temp_Dirichlet::Property_Temp_Dirichlet()
    {
        // set the property type
        mPropertyType = fem::Property_Type::TEMP_DIRICHLET;

        // set the useCoeff flag
        mUseCoeff = true;

    }

//------------------------------------------------------------------------------

    void Property_Temp_Dirichlet::val_coeff( Matrix< DDRMat > & aCoeff )
    {
        // FIXME set coeff to 5
        aCoeff.set_size( 1, 1, 5 );
    }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
