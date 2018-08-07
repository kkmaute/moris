#include "cl_HMR_Element.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
//------------------------------------------------------------------------------

        Element::Element(
                Background_Element_Base*  aElement ) :
                    mElement( aElement )
        {

        }
//------------------------------------------------------------------------------

        Element * Element::get_neighbor(
                moris::Cell< Element * > & aAllElementsOnProc,
                const luint              & aNeighborNumber )
        {
            // get neighbor of background element
            Background_Element_Base* tElement
                = mElement->get_neighbor( aNeighborNumber );

            // test if neighbor exists
            if ( tElement != NULL )
            {
                // test if neighbor is on the same level
                if ( tElement->get_level() == mElement->get_level() )
                {
                    return aAllElementsOnProc( tElement->get_memory_index() );
                }
                else
                {
                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------

        /**
         * set the T-Matrix flag
         */
        void
        Element::set_t_matrix_flag()
        {
            mElement->set_t_matrix_flag();
        }

//-------------------------------------------------------------------------------

        /**
         * unset the T-Matrix flag
         */
        void
        Element::unset_t_matrix_flag()
        {
            mElement->unset_t_matrix_flag();
        }

//-------------------------------------------------------------------------------

        /**
         * query the T-Matrix flag
         */
        bool
        Element::get_t_matrix_flag() const
        {
            return mElement->get_t_matrix_flag();
        }

//-------------------------------------------------------------------------------

     } /* namespace hmr */
} /* namespace moris */

