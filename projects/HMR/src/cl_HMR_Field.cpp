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
                const moris_id      aID,
                      hmr::Block  * aBlock ) :
                        mtk::Field(
                                aLabel,
                                aID,
                                aBlock,
                                nullptr ),
                        mMesh( aBlock->get_lagrange_mesh() )
        {
            mNodeValues = mMesh->create_field_data( aLabel );
        }

//------------------------------------------------------------------------------

        Field::~Field()
        {
        }

//------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
