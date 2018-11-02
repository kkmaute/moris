#include "cl_MTK_Abstract_Field.hpp"
#include "cl_MTK_Mesh.hpp"
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Mesh;
//------------------------------------------------------------------------------

        Abstract_Field::Abstract_Field(
                            const std::string & aLabel,
                            const EntityRank  & aRank,
                            const Mesh        * aMesh,
                            const uint          aDimension ) :
                                    mLabel( aLabel ),
                                    mRank( aRank ),
                                    mMesh( aMesh ),
                                    mDimension( aDimension )
        {

        }

//------------------------------------------------------------------------------

        Abstract_Field::Abstract_Field(
                            const std::string &           aLabel,
                            const EntityRank  &           aRank,
                            const std::shared_ptr< Mesh > aMesh,
                            const uint                    aDimension) :
                Abstract_Field( aLabel, aRank, aMesh.get(), aDimension )
        {

        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
