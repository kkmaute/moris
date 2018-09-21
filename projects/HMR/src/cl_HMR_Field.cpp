
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
                                aBlock ),
                        mMesh( aBlock->get_lagrange_mesh() ),
                        mFieldIndex( aBlock->get_lagrange_mesh()->create_field_data( aLabel ) )
        {

        }

//------------------------------------------------------------------------------

        Field::~Field()
        {
        }

//------------------------------------------------------------------------------


        Matrix< DDRMat > &
        Field::get_node_values()
        {
            return mMesh->get_field_data( mFieldIndex );
        }


        const Matrix< DDRMat > &
        Field::get_node_values() const
        {
            return mMesh->get_field_data( mFieldIndex );
        }

//------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
