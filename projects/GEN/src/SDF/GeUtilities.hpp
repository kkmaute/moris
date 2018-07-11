//
// Created by messe on 1/8/18.
//

#ifndef MORIS_GEUTILITIES_HPP
#define MORIS_GEUTILITIES_HPP

#include <fstream>
#include <iostream>

#include <cstdlib>
#include <cmath>
//#include <ctime>
#include <limits>

#include "typedefs.hpp" // COR/src
#include "assert.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "fn_sort.hpp" // LNA/src
#include "fn_dot.hpp" // LNA/src
#include "cl_Hierarchical_Mesh_Main.hpp" // STK/src/Heirarchical

#ifdef  MORIS_USE_32BIT
#define MORIS_GE_HUGE std::numeric_limits<double>::max()
#elif   MORIS_USE_64BIT
#define MORIS_GE_HUGE std::numeric_limits<long double>::max()
#endif
#define MORIS_GE_EPSILON 1e-7
namespace ge {

// =============================================================================
// LEXICAL FUNCTIONS
// =============================================================================

    template <typename T> T min(const T& aA, const T& aB)
    {
        return (aA < aB) ? (aA) : (aB);
    }

// -----------------------------------------------------------------------------

    template <typename T> T max(const T& aA, const T& aB)
    {
        return (aA > aB) ? (aA) : (aB);
    }

// -----------------------------------------------------------------------------

    template <typename T> T max(const T& aA, const T& aB, const T& aC)
    {
        return ge::max(ge::max(aA, aB), aC);
    }

// -----------------------------------------------------------------------------

    template <typename T> T min(const T& aA, const T& aB, const T& aC)
    {
        return ge::min(ge::min(aA, aB), aC);
    }

// -----------------------------------------------------------------------------

    template <typename T> T max(const T& aA, const T& aB, const T& aC, const T& aD)
    {
        return ge::max(ge::max(aA, aB), ge::max(aC, aD));
    }

// =============================================================================
// Basic Math
// =============================================================================

    template <typename T> T abs(const T& aValue)
    {
        return std::abs(aValue);
    }

// -----------------------------------------------------------------------------

    template <typename T> T round(const T aValue)
    {
        return std::round(aValue);
    }

// -----------------------------------------------------------------------------

    /**
     * @brief returns Pi
     *
     */
    constexpr moris::real pi(){
        return std::atan2(0,-1);
    }

// =============================================================================
// Vector Functions
// =============================================================================

    template <typename T> T dot(const moris::Mat<T>& aA, const moris::Mat<T>& aB)
    {
        //return moris::dot(aA, aB);
#if !defined(NDEBUG) || defined(DEBUG)
        moris::uint tLengthA = ge::max(aA.n_rows(), aA.n_cols());
        moris::uint tLengthB = ge::max(aB.n_rows(), aB.n_cols());
        moris::uint tSmallA  = ge::min(aA.n_rows(), aA.n_cols());
        moris::uint tSmallB  = ge::min(aB.n_rows(), aB.n_cols());

        MORIS_ASSERT ((tLengthA == tLengthB ) && (tSmallA == 1) &&  (tSmallB == 1),
                      "ge::dot: input vectors must be of same length and either row or col vector");

#else
        moris::uint tLengthA =  ge::max(aA.n_rows(),aA.n_cols());
#endif

        T aValue = 0;
        for(moris::uint i=0; i< tLengthA; ++i)
        {
            aValue += aA(i)*aB(i);
        }
        return aValue;
    }

// -----------------------------------------------------------------------------

    template <typename T> T norm(const moris::Mat<T>& aMatrix) {
        auto tMagnitude = ge::dot(aMatrix, aMatrix);
        return std::sqrt(tMagnitude);
        //return aMatrix.norm();
    }

// -----------------------------------------------------------------------------

    /**
    * @brief Returns the cross product of two vectors. If both vectors are of
    *        different datatype, return value corresponds to first datatype.
    *
    */
    template <typename T, typename U> moris::Mat<T> cross(const moris::Mat<T>& aA, const moris::Mat<U>& aB)
    {

#if !defined(NDEBUG) || defined(DEBUG)
        moris::uint tLengthA = ge::max(aA.n_rows(), aA.n_cols());
        moris::uint tLengthB = ge::max(aB.n_rows(), aB.n_cols());
        moris::uint tSmallA  = ge::min(aA.n_rows(), aA.n_cols());
        moris::uint tSmallB  = ge::min(aB.n_rows(), aB.n_cols());

        MORIS_ASSERT ((tLengthA == 3 ) && (tLengthB == 3 ) && (tSmallA == 1) &&  (tSmallB == 1),
                      "ge::cross: both vectors must and either a row or col vector of length 3");

#endif
        moris::Mat<T> aOut(3,1);
        aOut(0) = aA(1)*aB(2) - aA(2)*aB(1);
        aOut(1) = aA(2)*aB(0) - aA(0)*aB(2);
        aOut(2) = aA(0)*aB(1) - aA(1)*aB(0);
        return aOut;
    }

// -----------------------------------------------------------------------------

    /**
     * @brief returns the rotation matrix around the rotation vector and an angle
     *
     */
    template <typename T> moris::Mat<T> rotationMatrix(const moris::Mat<T>& aRotationVector, const T aAngle){
        T tCos = std::cos(aAngle);
        T tSin = std::sin(aAngle);
        T tCos_minusOne = tCos  -1;

        moris::Mat<T> aT(3,3);

        aT(0,0) = tCos-aRotationVector(0)*aRotationVector(0)*tCos_minusOne;
        aT(1,0) = aRotationVector(2)*tSin-aRotationVector(0)*aRotationVector(1)*tCos_minusOne;
        aT(2,0) = -aRotationVector(1)*tSin-aRotationVector(0)*aRotationVector(2)*tCos_minusOne;
        aT(0,1) = -aRotationVector(2)*tSin-aRotationVector(0)*aRotationVector(1)*tCos_minusOne;
        aT(1,1) = tCos-aRotationVector(1)*aRotationVector(1)*tCos_minusOne;
        aT(2,1) = aRotationVector(0)*tSin-aRotationVector(1)*aRotationVector(2)*tCos_minusOne;
        aT(0,2) = aRotationVector(1)*tSin-aRotationVector(0)*aRotationVector(2)*tCos_minusOne;
        aT(1,2) = -aRotationVector(0)*tSin-aRotationVector(1)*aRotationVector(2)*tCos_minusOne;
        aT(2,2) = tCos-aRotationVector(2)*aRotationVector(2)*tCos_minusOne;
        return aT;
    }

// -----------------------------------------------------------------------------

    /**
     * @brief returns the Distance between two points
     *
     * @param[in] aA    First Point
     * @param[in] aB    Second Point
     */

    template <typename T> T point_to_point_distance(const moris::Mat<T>& aA, const moris::Mat<T>& aB)
    {
#if !defined(NDEBUG) || defined(DEBUG)
        moris::uint tLengthA = ge::max(aA.n_rows(), aA.n_cols());
        moris::uint tLengthB = ge::max(aB.n_rows(), aB.n_cols());
        moris::uint tSmallA  = ge::min(aA.n_rows(), aA.n_cols());
        moris::uint tSmallB  = ge::min(aB.n_rows(), aB.n_cols());

        MORIS_ASSERT ((tLengthA == tLengthB ) && (tSmallA == 1) &&  (tSmallB == 1),
                      "ge::dot: input vectors must be of same length and either row or col vector");

#else
        moris::uint tLengthA =  ge::max(aA.n_rows(),aA.n_cols());
#endif

        moris::Mat<T> tDeltaX(tLengthA,1);
        for(moris::uint i=0; i< tLengthA; ++i){
            tDeltaX(i) = aA(i)-aB(i);
        }
        return ge::norm(tDeltaX);
    }

// =============================================================================
// Random Stuff
// =============================================================================
    moris::uint randomSeed() {
        std::ifstream file ("/dev/urandom", std::ios::binary);
        moris::uint tSeed;
        if (file.is_open())
        {
            char * memblock;
            int size = sizeof(moris::uint);
            memblock = new char [size];
            file.read (memblock, size);
            file.close();
            tSeed = *reinterpret_cast<int*>(memblock);
            delete[] memblock;
        }
        else
        {
            tSeed = time(NULL);
        }
        return tSeed;

    }

// -----------------------------------------------------------------------------

    /**
     * @brief returns a normalized pseudorandom vector
     *
     * @param[out] aVector     The random vector. Must be initialized already.
     */

    template <typename T> void randomVector(moris::Mat<T>& aVector){
#if !defined(NDEBUG) || defined(DEBUG)
        moris::uint tLength = ge::max(aVector.n_rows(), aVector.n_cols());
        moris::uint tSmall  = ge::min(aVector.n_rows(), aVector.n_cols());

        MORIS_ASSERT(((tLength == 3) && (tSmall == 1)),"randomVector: argument must be a 3x1 or 1x3 vector");
#endif
        std::srand(ge::randomSeed());

        for(moris::uint i=0; i<3; ++i){
            aVector(i) = std::rand();
        }
        T tNorm = ge::norm(aVector);
        for(moris::uint i=0; i<3; ++i){
            aVector(i) /= tNorm;
        }
    }

// -----------------------------------------------------------------------------

    /**
     * @brief returns a pseudorandom angle between -Pi and Pi
     *
     * @param[out] aAngle: The returned angle in rad.
     *
     */
    template <typename T> void randomAngle(T& aAngle){
        std::srand(ge::randomSeed());
        aAngle = ( ((T) std::rand())/RAND_MAX - 0.5)*2*ge::pi();
    }

// =============================================================================
// Special Functions
// =============================================================================

    // functionality similar to moris::unique, but regards epsilon environment
   template <typename T> moris::Mat<T> unique(moris::Mat<T>& aMatrix){
        // step 1: sort

        moris::uint tNRows = aMatrix.n_rows();
        moris::uint tNCols = aMatrix.n_cols();

        MORIS_ASSERT(((tNRows==1) || (tNCols==1)),
                     "ge::unique: Matrix must be a single row or single column");

        moris::uint tSize = tNRows*tNCols;

        moris::Mat<T> tSorted = moris::sort(aMatrix);

        // step 2: count unique entries
        moris::uint tUniqueCount = 1;
        for (moris::uint k=1; k<tSize; ++k){
            if((ge::abs(tSorted(k)-tSorted(k-1)) > MORIS_GE_EPSILON) &&
               (tSorted(k) != MORIS_GE_HUGE))
                ++tUniqueCount;
        }

        // step 3: write values in returned matrix

        if (tNRows==1)
            tNCols = tUniqueCount;
        else
            tNRows = tUniqueCount;
        moris::Mat<T> aUnique(tNRows, tNCols);

        tUniqueCount = 1;
        aUnique(0) = tSorted(0);

        for (moris::uint k=1; k<tSize; ++k){
             if((ge::abs(tSorted(k)-tSorted(k-1)) > MORIS_GE_EPSILON) &&
                (tSorted(k) != MORIS_GE_HUGE)){
                    aUnique(tUniqueCount) = tSorted(k);
                    ++tUniqueCount;
            }
        }

        return aUnique;

    }

// -----------------------------------------------------------------------------

    void TrianglePermutation(const moris::uint aZ, moris::uint& aX, moris::uint&aY){
        if(aZ==0){  // y-z
            aX = 1;
            aY = 2;
        } else if(aZ==1) { // x-z
            aX = 2;
            aY = 0;
        } else { // x-y
            aX = 0;
            aY = 1;
        }

    }

// -----------------------------------------------------------------------------

     template <typename T> T swap_byte_endian(T aValue)
     {
         T aOutValue;
         auto *tPointer = (char*) &aValue;
         auto *tOutPointer = (char*)&aOutValue;
         int size = sizeof(T);
         for(int i=0; i<size; i++)
         {
             tOutPointer[size - 1 - i] = tPointer[i];
         }
         return aOutValue;
     }

// -----------------------------------------------------------------------------

   /**
    * @brief convert a bitset to a vector which contains flagged indices
    *
    * @param[in] aBitset: BoostBitset to be converted
    *
    */
    moris::Mat< moris:: uint > BitsetToMat(const moris::BoostBitset aBitset)
    {
        moris::Mat< moris::uint > aVector(aBitset.count(), 1);

        moris::uint tCounter = 0;
        for(moris::uint k=0; k<aBitset.size(); ++k)
        {
            if (aBitset.test( k ))
            {
                aVector( tCounter ) = k;
                ++tCounter;
            }
        }
        return aVector;
    }

// -----------------------------------------------------------------------------
}
#endif //MORIS_GEUTILITIES_HPP
