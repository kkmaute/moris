#ifndef MORIS_LINALG_OP_OSTREAM_HPP_
#define MORIS_LINALG_OP_OSTREAM_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

namespace moris
{
template< typename T >
std::ostream &
operator<<(
std::ostream &          aOut,
moris::Mat< T > const & aA )
{
return aOut << aA.data();
}
}

#endif  /* MORIS_LINALG_OP_OSTREAM_HPP_ */

//#ifndef MORIS_LINALG_OP_OSTREAM_HPP_
//#define MORIS_LINALG_OP_OSTREAM_HPP_
//
//// MORIS library header files.
//#include "cl_Mat.hpp"
//#include "cl_Sp_Mat.hpp"
//
//// ----------------------------------------------------------------------------
//
//namespace moris
//{
//	template< typename T >
//	std::ostream &
//	operator<<(
//			std::ostream               & aOut,
//			moris::Base_Mat< T > const & aA )
//	{
//		return aOut << aA.data();
//	}
//}
//
//#endif  /* MORIS_LINALG_OP_OSTREAM_HPP_ */
