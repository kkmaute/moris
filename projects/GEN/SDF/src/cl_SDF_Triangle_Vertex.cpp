#include "cl_SDF_Triangle_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

            Triangle_Vertex::Triangle_Vertex(
                const moris_index        aIndex,
                const Matrix< DDRMat > & aNodeCoords ) :
                mIndex( aIndex ),
                mNodeCoords( aNodeCoords )
        {

        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

