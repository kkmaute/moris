/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Debug.cpp
 *
 */

#include "cl_Debug.hpp"
namespace moris
{

    Matrix< DDUMat >
    Debug::duplicate_row_check(Matrix< DDRMat >  & aCoord)
    {
        //explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf

        // moris::uint tz;
        moris::uint tNumInt = aCoord.n_rows();
        moris::uint tdim = aCoord.n_cols();
        moris::uint tCombination =  tNumInt*tNumInt; //boost::math::binomial_coefficient<double>(tNumInt, tNumInt-2); // Maximum combinations of duplicates
        Matrix< DDUMat >  duplicate_list(tCombination,2,UINT_MAX);
        moris::real tx;
        moris::uint tz = 0;
        moris::real ttol = 1.E-16;   // tolerance for checking nodes

        if (tCombination>UINT_MAX)
        {
            std::cout << "Size of the vector for duplicate check is to big" << '\n';
        }

        for (moris::uint  i = 0; i<tNumInt;i++){
            tx = aCoord(i,0);

            for (moris::uint  j = i+1; j<tNumInt;j++){
                if(std::abs(tx-aCoord(j,0))<ttol){

                    if(tdim==0){
                        duplicate_list(tz,0)=i; duplicate_list(tz,1)=j;
                        tz++;
                    }

                    if(std::abs(aCoord(i,1)-aCoord(j,1))<ttol && tdim>0){

                        if(tdim==1){
                            duplicate_list(tz,0)=i; duplicate_list(tz,1)=j;
                            tz++;
                        }

                        if(std::abs(aCoord(i,2)-aCoord(j,2))<ttol && tdim>1){
                            duplicate_list(tz,0)=i; duplicate_list(tz,1)=j;
                            tz++;
                        }
                    }
                }
            }
        }
        duplicate_list.resize(tz,2);
        return duplicate_list;
    }

    Matrix< DDUMat >
    Debug::duplicate_row_check(Matrix< DDUMat >  & aId)
    {
        //explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf

        // moris::uint tz;
        moris::uint tNumInt = aId.n_rows();
        moris::uint tCombination =  tNumInt*tNumInt;//boost::math::binomial_coefficient<double>(tNumInt, tNumInt-2); // Maximum combinations of duplicates
        Matrix< DDUMat >  duplicate_list(tCombination,2,UINT_MAX);
        moris::real tx;
        moris::uint tz = 0;
        moris::real ttol = 1.E-16;   // tolerance for checking nodes

        if (tCombination>UINT_MAX)
        {
            std::cout << "Size of the vector for duplicate check is to big" << '\n';
        }

        for (moris::uint  i = 0; i<tNumInt;i++){
            tx = aId(i);

            for (moris::uint  j = i+1; j<tNumInt;j++){
                if(std::abs(tx-aId(j))<ttol){

                    duplicate_list(tz,0)=i; duplicate_list(tz,1)=j;
                    tz++;
                }
            }
        }
        duplicate_list.resize(tz,2);
        return duplicate_list;
    }

    Matrix< DDUMat >
    Debug::duplicate_col_check(Matrix< DDUMat >  & aId)
    {
        //explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf

        // moris::uint tz;
        moris::uint tNumInt = aId.n_cols();
        moris::uint tCombination =  tNumInt*tNumInt;//boost::math::binomial_coefficient<double>(tNumInt, tNumInt-2); // Maximum combinations of duplicates
        Matrix< DDUMat >  duplicate_list(tCombination,2,UINT_MAX);
        moris::real tx;
        moris::uint tz = 0;
        moris::real ttol = 1.E-16;   // tolerance for checking nodes

        if (tCombination>UINT_MAX)
        {
            std::cout << "Size of the vector for duplicate check is to big" << '\n';
        }

        for (moris::uint  i = 0; i<tNumInt;i++){
            tx = aId(i);

            for (moris::uint  j = i+1; j<tNumInt;j++){
                if(std::abs(tx-aId(j))<ttol){

                    duplicate_list(tz,0)=i; duplicate_list(tz,1)=j;
                    tz++;
                }
            }
        }
        duplicate_list.resize(tz,2);
        return duplicate_list;
    }

    Matrix< DDUMat >
    Debug::duplicate_row_check(Matrix< DDUMat >  & aId1,
            Matrix< DDUMat >  & aId2)
    {
        //explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf

        // moris::uint tz;
        moris::uint tId1 = aId1.n_rows();
        moris::uint tId2 = aId2.n_rows();
        moris::uint tCombination =  (tId1+tId2)*(tId1+tId2);//boost::math::binomial_coefficient<double>(tId1+tId2, tId1+tId2-2); // Maximum combinations of duplicates
        Matrix< DDUMat >  duplicate_list(tCombination,2,UINT_MAX);
        Matrix< DDUMat >  complete_list(tId1+tId2,2,UINT_MAX);
        moris::real tx;
        moris::uint tz = 0;
        moris::real ttol = 1.E-16;   // tolerance for checking nodes

        if (tCombination>UINT_MAX)
        {
            std::cout << "Size of the vector for duplicate check is to big" << '\n';
        }

        for (moris::uint  i = 0; i<tId1;i++){
            complete_list(i,0)=aId1(i,0); complete_list(i,1)=aId1(i,1);
        }
        for (moris::uint  i = 0; i<tId2;i++){
            complete_list(i+tId1,0)=aId2(i,0); complete_list(i+tId1,1)=aId2(i,1);
        }

        for (moris::uint i = 0; i<tId1+tId2;i++){
            tx = complete_list(i,0);

            for (moris::uint  j = i+1; j<tId1+tId2;j++){
                if(std::abs(tx-complete_list(j,0))<ttol){
                    if(complete_list(i,1)-complete_list(j,1)==0){
                        duplicate_list(tz,0)=complete_list(i,0); duplicate_list(tz,1)=complete_list(i,1);
                        tz++;
                    }
                }
            }

        }
        duplicate_list.resize(tz,2);
        return duplicate_list;
    }

    Matrix< DDUMat >
    Debug::duplicate_row_check_problems(Matrix< DDUMat >  & aId1,
            Matrix< DDUMat >  & aId2)
    {
        moris::uint tId1 = aId1.n_rows();
        moris::uint tId2 = aId2.n_rows();
        moris::uint tCombination =  (tId1+tId2)*(tId1+tId2);//boost::math::binomial_coefficient<double>(tId1, tId1-2); // Maximum combinations of duplicates
        Matrix< DDUMat >  duplicate_list(tCombination,2,UINT_MAX);
        Matrix< DDUMat >  problem_list(tId1,2,UINT_MAX);
        Matrix< DDUMat >  position(tCombination,1,UINT_MAX);
        moris::uint tz = 0;
        moris::uint ty = 0;

        if (tCombination>UINT_MAX)
        {
            std::cout << "Size of the vector for duplicate check is to big" << '\n';
        }

        for (moris::uint i = 0; i<tId1;i++)
        {
            ty = 1;
            for (moris::uint  j = 0; j<tId2;j++)
            {
                if( aId1(i,0) == aId2(j,0) && aId1(i,1) == aId2(j,1) )
                {
                    position(tz,0) = i;
                    tz++;
                }
            }
        }
        position.resize(tz,1);
        tz = 0;
        for (moris::uint i = 0; i<tId1;i++)
        {
            ty = 1;
            for (moris::uint  j = 0; j<position.n_rows();j++)
            {
                if(position(j) == i)
                {
                    ty = 0;
                }
            }
            if( ty == 1)
            {
                problem_list(tz,0) = aId1(i,0);   problem_list(tz,1) = aId1(i,1);
                tz++;
            }
        }
        problem_list.resize(tz,2);
        return problem_list;
    }
}

