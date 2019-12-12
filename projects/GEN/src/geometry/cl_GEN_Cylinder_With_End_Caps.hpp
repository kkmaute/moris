/*
 * cl_GEN_Cylinder_With_End_Caps.hpp
 *
 *  Created on: Nov 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CYLINDER_WITH_END_CAPS_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CYLINDER_WITH_END_CAPS_HPP_


#include "../geometry/cl_GEN_Geometry.hpp"
#include "../geometry/fn_bounding_box.hpp"
#include "HDF5_Tools.hpp"

#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace ge
{

class GEN_CylinderWithEndCaps : public GEN_Geometry
{
private:    // member data
    uint mNumberOfFibers;
    moris::Cell< Matrix< DDRMat > > mAllFibers;
    bool mFileHasBeenRead;

    //------------------------------------------------------------------------------
public:
    GEN_CylinderWithEndCaps(  ) :
        mNumberOfFibers( 448 ), // total number of fibers in the data file = 448
        mFileHasBeenRead( false )
    {
            mAllFibers.resize( mNumberOfFibers );
    }

    GEN_CylinderWithEndCaps( uint aNumberOfFibers) :
        mNumberOfFibers( aNumberOfFibers ),
        mFileHasBeenRead( false )
    {
            mAllFibers.resize( aNumberOfFibers );
    }

    ~GEN_CylinderWithEndCaps()
    {

    }
    //------------------------------------------------------------------------------
    bool is_analytic() const
    {
        return true;
    }

    void get_dphi_dp_size(moris::size_t & aNumRows, moris::size_t & aNumCols) const
    {
        aNumRows = 1;
        aNumCols = 1;
    }
    //------------------------------------------------------------------------------
    real cylinderWithEndCaps( const moris::Matrix< moris::DDRMat > & aPoint )   // full function with smoothing
    {
        if ( mFileHasBeenRead == false )
        {
            readFiberData();
            mFileHasBeenRead=true;
        }

        //------------------------------------------------------------------------------
        Matrix<DDRMat> xa(1,3);
        Matrix<DDRMat> xb(1,3);
        Matrix<DDRMat> xn(1,3);
        Matrix<DDRMat> xx(1,3);

        const real KSbeta=-10.0;
        real KSvalue = 0.0;

        for( uint k=0; k<mNumberOfFibers; ++k)
        {
            uint tNumBars = mAllFibers(k).n_cols()-1;

            uint tPointCount = 0;

            for (uint ib=0; ib<tNumBars; ++ib)
            {
                xa.set_row(0,mAllFibers(k).get_column(tPointCount));
                xb.set_row(0,mAllFibers(k).get_column(tPointCount+1));

                if( outside_bounding_box_check( aPoint,xa,xb,0.5) == false )    // bounding box check
                {
                    real rada = 0.18;
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

                        xx.set_row(0,xx-xa);
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
                }
                tPointCount++;
            }
        }

        return  std::abs(KSvalue) > 1e-9 ? 1.0/KSbeta*std::log(KSvalue) : 1.0;
    }
    //------------------------------------------------------------------------------

    real create_cylinder( const moris::Matrix< moris::DDRMat > & aPoint,
                          const uint aFiberIndex,
                          const uint aCylinderIndex )
    {
        if ( mFileHasBeenRead == false )
        {
            readFiberData();
            mFileHasBeenRead = true;
        }
        //------------------------------------------------------------------------------
        Matrix<DDRMat> xa(1,3);
        Matrix<DDRMat> xb(1,3);
        Matrix<DDRMat> xn(1,3);
        Matrix<DDRMat> xx(1,3);

//        const real KSbeta=-10.0;
        real KSvalue = 0.0;

            xa.set_row(0,mAllFibers(aFiberIndex).get_column(aCylinderIndex));
            xb.set_row(0,mAllFibers(aFiberIndex).get_column(aCylinderIndex+1));

            real rada = 0.18;
            real radb = 0.18;

            xn.set_row(0,xb);

            xn.set_row(0,xn-xa);

            real len = norm(xn);

            real tTemp = 1.0/len;
            xn(0) = tTemp*xn(0);
            xn(1) = tTemp*xn(1);
            xn(2) = tTemp*xn(2);

            real pc = dot(xn,aPoint);

            real s = pc - dot(xn,xa);

            if (s >= 0 && s <= len)
            {
                xx.set_row(0,aPoint);

                xx.set_row(0,xx-xa);
                //------------------------------------------------------------------------------
                xx(0) = -s*xn(0)+xx(0);
                xx(1) = -s*xn(1)+xx(1);
                xx(2) = -s*xn(2)+xx(2);

                real rad = rada + (radb-rada)*s/len;

                KSvalue = norm(xx)-rad;

            }
            else if (s < 0.0)
            {
                xx.set_row(0,aPoint);

                xx.set_row(0,xx-xa);

                KSvalue = norm(xx)-rada;

            }
            else
            {
                xx.set_row(0,aPoint);

                xx(0) = (xx(0)-xb(0));
                xx(1) = (xx(1)-xb(1));
                xx(2) = (xx(2)-xb(2));

                KSvalue = norm(xx)-radb;

            }

        return KSvalue;
    }

    void midPoint_and_BB_dims( const uint aFiberNumber,
                               const uint aCylinderIndex,
                               Matrix<DDRMat> & aMidPoint,
                               Matrix<DDRMat> & aLength,
                               const real aRadius = 0.18,
                               const real aBuffer = 0.5 )
    {
        if ( mFileHasBeenRead == false )
        {
            readFiberData();
            mFileHasBeenRead = true;
        }
        Matrix< DDRMat > tFiberEnd00 = mAllFibers(aFiberNumber).get_column( aCylinderIndex );
        Matrix< DDRMat > tFiberEnd01 = mAllFibers(aFiberNumber).get_column( aCylinderIndex+1 );

        aMidPoint.set_size( 1, 3, 0.0 );
        aLength.set_size( 1, 3, 0.0 );
        aMidPoint(0) = (tFiberEnd00(0)+tFiberEnd01(0))/2;
        aMidPoint(1) = (tFiberEnd00(1)+tFiberEnd01(1))/2;
        aMidPoint(2) = (tFiberEnd00(2)+tFiberEnd01(2))/2;

    //    aLength(0) = aMidPoint(0) - tFiberEnd00(0) + aRadius + aBuffer;
    //    aLength(1) = aMidPoint(1) - tFiberEnd00(1) + aRadius + aBuffer;
    //    aLength(2) = aMidPoint(2) - tFiberEnd00(2) + aRadius + aBuffer;
        aLength(0) = -aMidPoint(0) + tFiberEnd01(0) + aRadius + aBuffer;
        aLength(1) = -aMidPoint(1) + tFiberEnd01(1) + aRadius + aBuffer;
        aLength(2) = -aMidPoint(2) + tFiberEnd01(2) + aRadius + aBuffer;
    }

    uint get_number_of_cylinders( const uint aFiberNumber )
    {
        if ( mFileHasBeenRead == false )
        {
            readFiberData();
            mFileHasBeenRead = true;
        }
        return mAllFibers(aFiberNumber).n_cols()-1;
    }

    //------------------------------------------------------------------------------

private:    // private functions
    //------------------------------------------------------------------------------

    void readFiberData()
    {
        std::string tMorisRoot = std::getenv("MORISROOT");     // get root from environment
        std::string tHdf5FilePath = tMorisRoot + "/projects/GEN/test/hdf5_files/allFibers.hdf5" ;
        hid_t tFileID = open_hdf5_file( tHdf5FilePath );
        herr_t tStatus = 0; // error handler

        for( uint k=0; k<mNumberOfFibers; ++k)
        {
            std::string tK = std::to_string(k);
            std::string tSetLabel = "fiber" + tK;     // create label
            load_matrix_from_hdf5_file( tFileID, tSetLabel, mAllFibers(k), tStatus );  // read solution from file
        }

        close_hdf5_file( tFileID ); // close file
    }

};

}   // end ge namespace
}   // end moris namespace


#endif /* PROJECTS_GEN_SRC_NEW_GEOMETRY_CL_GEN_CYLINDER_WITH_END_CAPS_HPP_ */
