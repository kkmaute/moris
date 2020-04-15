#include "cl_Cell.hpp"

// C++ header files.
#include <vector>
#include <algorithm> // for unique
#include <iostream>

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "assert.hpp"

namespace moris
{

// ------------------------------------------------------------------------------------ //
// --- FREE FUNCTIONS ----------------------------------------------------------------- //
// ------------------------------------------------------------------------------------ //

// ------------------------------------------------------------------------------------ //
//template< typename T >
//void
//unique( Cell< T > & aCell )
//{
//    // get ref to data
//    std::vector< T > & tVec = aCell.data();
//
//    // sort data
//    std::sort( tVec.begin(), tVec.end() );
//
//    // trim vector
//    tVec.erase( std::unique( tVec.begin(), tVec.end() ), tVec.end() );
//}


// ------------------------------------------------------------------------------------ //
//https://stackoverflow.com/questions/25921706/creating-a-vector-of-indices-of-a-sorted-vector
// extended to create a unique with the first index of unique values
//template< typename T >
//Cell<moris::moris_index>
//unique_index( Cell< T > & aCell )
//{
//    std::vector<T> x = aCell.data();
//
//    std::vector<int> y(x.size());
//    std::size_t n(0);
//    std::generate(std::begin(y), std::end(y), [&]{ return n++; });
//
//    std::sort(  std::begin(y),
//                std::end(y),
//                [&](int i1, int i2) { return x[i1] < x[i2]; } );
//
//    Cell<moris::moris_index> tUniqueInds;
//    for(moris::uint  i = 0; i < aCell.size(); i++)
//    {
//        if(i == 0)
//        {
//            tUniqueInds.push_back( y[i] );
//        }
//
//        else if(aCell(y[i-1]) != aCell(y[i]) )
//        {
//            tUniqueInds.push_back( y[i] );
//        }
//    }
//
//    return tUniqueInds;
//
//}


// ------------------------------------------------------------------------------------ //
/*!
 * Iterates through cell and prints each cell.
 * Will only work on data types that allow std::cout calls
 */
//template< typename T >
//void
//print(Cell< T > const & aCell,
//      std::string aStr)
//{
//    std::cout<<"Cell Name: "<<aStr<<"\n";
//    std::cout<<"Number of entries = "<<aCell.size()<<"\n";
//    for(moris::uint  i = 0; i <aCell.size(); i++)
//    {
//        std::cout<<aCell(i)<<"\n";
//    }
//
//    std::cout<<std::endl;
//}


// ------------------------------------------------------------------------------------ //
moris::Cell<char>
string_to_char(moris::Cell<std::string>& strings)
{
    moris::Cell<char> cstrings;
    cstrings.reserve(strings.size());
    for(std::string s: strings)
    {
        for(size_t i = 0; i < strlen(s.c_str()); ++i)
        {
            cstrings.push_back(s.c_str()[i]);
        }
    }

    return cstrings;
}

} // namespace moris


