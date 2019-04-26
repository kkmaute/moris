
#include "catch.hpp"
#include "fn_equal_to.hpp"

#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Element_Block.hpp"              //FEM/INT/src
#undef protected
#undef private

namespace moris
{
    namespace fem
    {
        TEST_CASE( "Element_Block", "[moris],[fem],[ElementBlock]" )
        {
//           MSI::Equation_Block tElementBlock = new Element_Block();
//
//           tElementBlock.mIWGs = ;
//
//           tElementBlock.create_unique_dof_type_lists();

//            CHECK( equal_to( (tElementBlock.mEqnObjDofTypeList( 0 ), 0 ) );



        }/* TEST_CASE */
    }/* namespace fem */
}/* namespace moris */
