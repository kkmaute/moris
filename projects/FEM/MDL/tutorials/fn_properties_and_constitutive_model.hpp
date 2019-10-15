/*
 * fn_properties_and_constitutive_model.hpp
 *
 *  Created on: Oct 14, 2019
 *      Author: noel
 */

#ifndef PROJECTS_FEM_MDL_TUTORIALS_FN_PROPERTIES_AND_CONSTITUTIVE_MODEL_HPP_
#define PROJECTS_FEM_MDL_TUTORIALS_FN_PROPERTIES_AND_CONSTITUTIVE_MODEL_HPP_

namespace moris
{


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
