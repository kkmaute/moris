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
#include "fn_bounding_box.hpp"

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
uint tNumFiber = 448;
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






/*
//------------------------------------------------------------------------------
// fiber 0
std::string tMorisRoot = std::getenv("MORISROOT");     // get root from environment
//std::string tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/fiber0.hdf5" ;
std::string tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/fibers.hdf5" ;
hid_t tFileID = open_hdf5_file( tHdf5FilePath );
herr_t tStatus = 0; // error handler
std::string tSetLabel = "fiber0";     // create label
Matrix<DDRMat> tFiber_0;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_0, tStatus );  // read solution from file
close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Matrix<DDRMat> lsValue(1,11);   // n_cols = number of fibers

//uint tNumBars = aConst(6);
    uint tNumBars = tFiber_0.n_cols()-1;

Matrix<DDRMat> xa(1,3);

Matrix<DDRMat> xb(1,3);

Matrix<DDRMat> xn(1,3);

Matrix<DDRMat> xx(1,3);

const real KSbeta=-10.0;

real KSvalue = 0.0;
uint tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{

    // coordinates of point A

//    xa(0,0) = aConst(0);
//    xa(1,0) = aConst(1);
//    xa(2,0) = aConst(2);
    xa.set_row(0,tFiber_0.get_column(tPointCount));

    // radius at point A

//    real rada = aConst(7);
    real rada = 0.18;

    // coordinates of point B

//    xb(0,0) = aConst(3);
//    xb(1,0) = aConst(4);
//    xb(2,0) = aConst(5);
    xb.set_row(0,tFiber_0.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if( check_if_in_bounding_box( aPoint,xa,xb,10 ) )
    {
        continue;
    }
//====================================================================
    // radius at point B

//    real radb = aConst(8);
    real radb = 0.18;

    // xn = xb

    xn.set_row(0,xb);

    // xn = xn - xa

    xn.set_row(0,xn-xa);

    // length of xn

    real len = norm(xn);

    // normalize xn
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

//real lsValue = 1.0/KSbeta*std::log(KSvalue);
//lsValue(0) = -1.0/KSbeta*std::log(KSvalue);
real lsValue = 1.0/KSbeta*std::log(KSvalue);
*/
/*
//------------------------------------------------------------------------------
// fiber 1
//tMorisRoot = std::getenv("MORISROOT");     // get root from environment
//tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/fiber1.hdf5" ;
//tFileID = open_hdf5_file( tHdf5FilePath );
tSetLabel = "fiber1";     // create label
Matrix<DDRMat> tFiber_1;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_1, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_1.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_1.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_1.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(1) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 2
//tMorisRoot = std::getenv("MORISROOT");     // get root from environment
//tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/fiber2.hdf5" ;
//tFileID = open_hdf5_file( tHdf5FilePath );
tSetLabel = "fiber2";     // create label
Matrix<DDRMat> tFiber_2;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_2, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_2.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_2.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_2.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(2) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 3
tSetLabel = "fiber3";     // create label
Matrix<DDRMat> tFiber_3;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_3, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_3.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_3.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_3.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(3) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 4
tSetLabel = "fiber4";     // create label
Matrix<DDRMat> tFiber_4;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_4, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_4.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_4.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_4.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(4) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 5
tSetLabel = "fiber5";     // create label
Matrix<DDRMat> tFiber_5;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_5, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_5.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_5.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_5.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(5) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 6
tSetLabel = "fiber6";     // create label
Matrix<DDRMat> tFiber_6;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_6, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_6.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_6.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_6.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(6) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 7
tSetLabel = "fiber7";     // create label
Matrix<DDRMat> tFiber_7;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_7, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_7.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_7.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_7.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(7) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 8
tSetLabel = "fiber8";     // create label
Matrix<DDRMat> tFiber_8;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_8, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_8.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_8.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_8.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(8) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 9
tSetLabel = "fiber9";     // create label
Matrix<DDRMat> tFiber_9;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_9, tStatus );  // read solution from file
//close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_9.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_9.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_9.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(9) = -1.0/KSbeta*std::log(KSvalue);
//------------------------------------------------------------------------------
// fiber 10
tSetLabel = "fiber10";     // create label
Matrix<DDRMat> tFiber_10;
load_matrix_from_hdf5_file( tFileID, tSetLabel, tFiber_10, tStatus );  // read solution from file
close_hdf5_file( tFileID ); // close file
//------------------------------------------------------------------------------
tNumBars = tFiber_10.n_cols()-1;
KSvalue = 0.0;
tPointCount = 0;
for (uint ib=0; ib<tNumBars; ++ib)
{
xa.set_row(0,tFiber_10.get_column(tPointCount));
real rada = 0.18;
xb.set_row(0,tFiber_10.get_column(tPointCount+1));
//====================================================================
// bounding box check
    if(!check_if_in_bounding_box(aPoint,xa,xb,75))
    {
        continue;
    }
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
lsValue(10) = -1.0/KSbeta*std::log(KSvalue);

//------------------------------------------------------------------------------
//return lsValue;

return lsValue.max();
*/
}

}   // ge
}   // moris


#endif /* PROJECTS_GEN_SRC_FN_CYLINDER_WITH_END_CAPS_HPP_ */
