/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_to_integer.hpp
 *
 */

#ifndef MORIS_STRING_FN_TO_INTEGER_HPP_
#define MORIS_STRING_FN_TO_INTEGER_HPP_

// C++ header files.
#include <string>

// MORIS header files.
#include "core.hpp"

namespace moris
{
namespace string
{
	/**
	 * Returns an integer from a string.
	 *
	 * @note To be used in switch statements that require a constexpr std::string.
	 * @see [Evaluate a string with a switch in C++]
	 * (http://stackoverflow.com/questions/16388510/evaluate-a-string-with-a-switch-in-c)
	 *
	 * @since 1.0
	 *
	 * @param[in] string String.
	 * @param[in] h
	 * @return Integer equivalent of a string.
	 */
	constexpr
	moris::size_t
	to_integer(
			char const *        string,
			moris::sint const & h = 0)
	{
		return ! string[h] ? 5381 : (moris::string::to_integer(string, h + 1) * 33) ^ string[h];
	}

	/**
	 * Overloaded moris::string::to_integer.
	 *
	 * @since 1.0
	 */
//	constexpr
//	moris::size_t
//	to_integer(
//			std::string const & string,
//			moris::sint const & h = 0)
//	{
//		return moris::string::to_integer(string.c_str(), h);
//	}

}	// namespace string
}	// namespace moris

#endif	/* MORIS_STRING_FN_TO_INTEGER_HPP_ */

