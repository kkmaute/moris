/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tet_56.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_56_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_56_HPP_

// MRS/COR/src
#include "typedefs.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Coeffs.hpp"

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_56 >::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_56 >::get_number_of_points()
        {
            return 56;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_56 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { //
                        0.9551438045408220,
                        0.0149520651530592,
                        0.0149520651530592,
                        0.0149520651530592,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.5965649956210170,
                        0.1344783347929940,
                        0.1344783347929940,
                        0.1344783347929940 },
                { //
                        0.0149520651530592,
                        0.9551438045408220,
                        0.0149520651530592,
                        0.0149520651530592,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1344783347929940,
                        0.5965649956210170,
                        0.1344783347929940,
                        0.1344783347929940 },
                { //
                        0.0149520651530592,
                        0.0149520651530592,
                        0.0149520651530592,
                        0.9551438045408220,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.0340960211962615,
                        0.0340960211962615,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.1518319491659370,
                        0.7799760084415400,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.0462051504150017,
                        0.0462051504150017,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.5526556431060170,
                        0.3549340560639790,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.0055147549744775,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.2281904610687610,
                        0.2281904610687610,
                        0.5381043228880020,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.0992057202494530,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.3523052600879940,
                        0.3523052600879940,
                        0.1961837595745600,
                        0.1344783347929940,
                        0.1344783347929940,
                        0.1344783347929940,
                        0.5965649956210170 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_56 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { //
                        0.0010373112336140,
                        0.0010373112336140,
                        0.0010373112336140,
                        0.0010373112336140,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0096016645399480,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0164493976798232,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0153747766513310,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0293520118375230,
                        0.0366291366405108,
                        0.0366291366405108,
                        0.0366291366405108,
                        0.0366291366405108 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_56_HPP_ */