/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_properties_and_constitutive_model.hpp
 *
 */

#ifndef PROJECTS_FEM_MDL_TUTORIALS_FN_PROPERTIES_AND_CONSTITUTIVE_MODEL_HPP_
#define PROJECTS_FEM_MDL_TUTORIALS_FN_PROPERTIES_AND_CONSTITUTIVE_MODEL_HPP_

namespace moris
{

moris::real
Plane4MatMDL1(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.1;
    moris::real mYC = 0.1;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return ( mNx*( aPoint(0)-mXC ) + mNy*( aPoint(1)-mYC ) );
}

moris::real
Circle4MatMDL(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXCenter = 0.01;
    moris::real mYCenter = 0.01;
    moris::real mRadius = 0.47334;

    return  (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);
}

Matrix< DDRMat > tConstValFunction2MatMDL( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                           moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                           moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                           fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_bulk_properties(fem::Property_User_Defined_Info & aConductivityProp,
                       fem::Property_User_Defined_Info & aHeatLoadProp)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 2 );
    tProps( 0 )( 0 ) = aConductivityProp;
    tProps( 0 )( 1 ) = aHeatLoadProp;

    return tProps;
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_dirichlet_properties(fem::Property_User_Defined_Info & aConductivity,
                            fem::Property_User_Defined_Info & aDirchletTemp)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 2 );
    tProps( 0 )( 0 ) = aConductivity;
    tProps( 0 )( 1 ) = aDirchletTemp;

    return tProps;

}

Cell< Cell< fem::Property_User_Defined_Info > >
create_neumann_properties(fem::Property_User_Defined_Info & aNeumannFlux)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(1);
    tProps.resize( 1 );
    tProps( 0 ).resize( 1 );
    tProps( 0 )( 0 ) = aNeumannFlux;

    return tProps;
}

Cell< Cell< fem::Property_User_Defined_Info > >
create_interface_properties(fem::Property_User_Defined_Info & aMasterCond,
                            fem::Property_User_Defined_Info & aSlaveCond)
{
    Cell< Cell< fem::Property_User_Defined_Info > > tProps(2);
    tProps( 0 ).resize( 1 );
    tProps( 0 )( 0 ) = aMasterCond;
    tProps( 1 ).resize( 1 );
    tProps( 1 )( 0 ) = aSlaveCond;

    return tProps;
}

fem::Constitutive_User_Defined_Info
create_diff_lin_constitutive_info()
{
    return fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO, {{ MSI::Dof_Type::TEMP }}, { fem::Property_Type::CONDUCTIVITY } );
}

Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_bulk_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aDiffLinConst )
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(1);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aDiffLinConst;

    return tConstitutiveUserDefInfo;
}

Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_dbc_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aDiffLinConst )
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(1);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aDiffLinConst;

    return tConstitutiveUserDefInfo;
}

Cell< Cell< fem::Constitutive_User_Defined_Info > >
create_interface_diff_lin_constitutive( fem::Constitutive_User_Defined_Info & aMasterDiffLinConst,
                                        fem::Constitutive_User_Defined_Info & aSlaveDiffLinConst)
{
    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefInfo(2);
    tConstitutiveUserDefInfo(0).resize(1);
    tConstitutiveUserDefInfo(0)(0) = aMasterDiffLinConst;

    tConstitutiveUserDefInfo(1).resize(1);
    tConstitutiveUserDefInfo(1)(0) = aSlaveDiffLinConst;

    return tConstitutiveUserDefInfo;
}

}

#endif /* PROJECTS_FEM_MDL_TUTORIALS_FN_PROPERTIES_AND_CONSTITUTIVE_MODEL_HPP_ */

