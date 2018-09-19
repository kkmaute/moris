
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Block.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Field::Field(
                const std::string & aLabel,
                      hmr::Block  * aBlock ) :
                        mtk::Field(
                                aLabel,
                                aBlock,
                                nullptr )
        {
            // fixme: Mesh is not supposed to store data
            //        this is preliminary until new MTK can write Exodus properly
            mNodeValues = aBlock->get_lagrange_mesh()->create_field_data( aLabel );
        }

//------------------------------------------------------------------------------

        Field::~Field()
        {
        }

//------------------------------------------------------------------------------


    } /* namespace hmr */
} /* namespace moris */
