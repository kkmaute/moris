/*
 * fn_cylinder_with_end_caps.hpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_FN_CYLINDER_WITH_END_CAPS_HPP_
#define PROJECTS_GEN_SRC_FN_CYLINDER_WITH_END_CAPS_HPP_

#include "HDF5_Tools.hpp"

#include "fn_dot.hpp"
#include "fn_norm.hpp"
//#include "fn_bounding_box.hpp"

namespace moris{
namespace ge{

/*
 * @brief creates a cylinder with endcaps
 * @param[in] aPoint - point at which to compute LS value at
 * @param[in] aConst - list of parameters
 *                     *aConst(0) to aConst(2) = coordinates of point A
 *                     *aConst(3) to aConst(5) = coordinates of point B
 *                     *aConst(6) = number of 'bars' i.e. pairs of points
 *                     *aConst(7) = radius of point A
 *                     *aConst(8) = radius 0f point B
 *
 * @param[out] tLSVal - LS val
 */
inline
//real cylinderWithEndCaps( const moris::Matrix< moris::DDRMat > & aPoint, const moris::Cell< moris::real > aConst )
real cylinderWithEndCaps( const moris::Matrix< moris::DDRMat > & aPoint )
{
uint tNumFiber = 100;
Matrix<DDRMat> lsValue(1,tNumFiber);   // n_cols = number of fibers
for(uint k=0; k<tNumFiber; ++k)
{
    //------------------------------------------------------------------------------
    std::string tMorisRoot = std::getenv("MORISROOT");     // get root from environment
    std::string tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/allFibers.hdf5" ;
    hid_t tFileID = open_hdf5_file( tHdf5FilePath );
    herr_t tStatus = 0; // error handler
    std::string tIter = std::to_string( k );
    std::string tSetLabel = "fiber"+tIter;     // create label
    Matrix<DDRMat> tFiber;
    load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber, tStatus );  // read solution from file
    close_hdf5_file( tFileID ); // close file
    //------------------------------------------------------------------------------
    uint tNumBars = tFiber.n_cols()-1;

    Matrix<DDRMat> xa(1,3);
    Matrix<DDRMat> xb(1,3);
    Matrix<DDRMat> xn(1,3);
    Matrix<DDRMat> xx(1,3);

    const real KSbeta=-10.0;
    real KSvalue = 0.0;
    uint tPointCount = 0;

    for (uint ib=0; ib<tNumBars; ++ib)
    {
        xa.set_row(0,tFiber.get_column(tPointCount));

        real rada = 0.18;

        xb.set_row(0,tFiber.get_column(tPointCount+1));
    //====================================================================
    // bounding box check
//        if( check_if_in_bounding_box( aPoint,xa,xb,10 ) )
//        {
//            continue;
//        }
    //====================================================================
        real radb = 0.18;

        xn.set_row(0,xb);

        xn.set_row(0,xn-xa);

        real len = norm(xn);

        real tTemp=1.0/len;
        xn(0) = tTemp*xn(0);
        xn(1) = tTemp*xn(1);
        xn(2) = tTemp*xn(2);

        real pc = dot(xn,aPoint);

        real s = pc - dot(xn,xa);

        if (s >= 0 && s <= len)
        {
            xx.set_row(0,aPoint);

            xx.set_row(0,xn-xa);
            //------------------------------------------------------------------------------
            xx(0) = -s*xn(0)+xx(0);
            xx(1) = -s*xn(1)+xx(1);
            xx(2) = -s*xn(2)+xx(2);

            real rad = rada + (radb-rada)*s/len;

            real lscyl = norm(xx)-rad;

            KSvalue += std::exp(KSbeta*lscyl);
        }
        else if (s < 0.0)
        {
            xx.set_row(0,aPoint);

            xx.set_row(0,xx-xa);

            real lscap1=norm(xx)-rada;

            KSvalue += std::exp(KSbeta*lscap1);
        }
        else
        {
            xx.set_row(0,aPoint);

            xx(0) = (xx(0)-xb(0));
            xx(1) = (xx(1)-xb(1));
            xx(2) = (xx(2)-xb(2));

            real lscap2=norm(xx)-radb;

            KSvalue += std::exp(KSbeta*lscap2);
       }
    tPointCount++;
    }
    lsValue(k) = -1.0/KSbeta*std::log(KSvalue);
}
return lsValue.max();

}

}   // ge
}   // moris


#endif /* PROJECTS_GEN_SRC_FN_CYLINDER_WITH_END_CAPS_HPP_ */
