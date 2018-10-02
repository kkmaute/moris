
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Block.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Field::Field(
                const std::string             & aLabel,
                std::shared_ptr< mtk::Mesh >    aMesh,
                std::shared_ptr< Database >     aDatabase,
                Lagrange_Mesh_Base *            aLagrangeMesh ) :
                        mtk::Field( aLabel, aMesh ),
                        mDatabase( aDatabase ),
                        mLagrangeMesh( aLagrangeMesh ),
                        mFieldIndex( aLagrangeMesh->create_field_data( aLabel ) )

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
            return mLagrangeMesh->get_field_data( mFieldIndex );
        }


        const Matrix< DDRMat > &
        Field::get_node_values() const
        {
            return mLagrangeMesh->get_field_data( mFieldIndex );
        }

//------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */
