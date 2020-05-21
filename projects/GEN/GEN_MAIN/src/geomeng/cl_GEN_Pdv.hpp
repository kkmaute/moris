#ifndef MORIS_CL_GEN_PDV_Type_HPP_
#define MORIS_CL_GEN_PDV_Type_HPP_

// GEN_MAIN
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv
        {
        public:
            bool mIsActive = true;

        private:
            real mValue; // PDV value

        public :
            /**
             * constructor
             * @param[ in ] aFieldPointer a GEN Field pointer
             * @param[ in ] aEntityIndex  an index to the associated entity (so the Field returns the correct value)
             */
            Pdv(std::shared_ptr< GEN_Field > aFieldPointer,
                moris_index                  aEntityIndex );

            /**
             * constructor
             * @param[ in ] aPdvVal a value for the pdv
             */
            Pdv(moris::real aPdvVal );

            /**
             * trivial destructor
             */
            ~Pdv();

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            virtual real get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * Get the PDV sensitivity with respect to ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @param aSensitivities Matrix of sensitivities to be returned
             */
            virtual void get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }   // end ge namespace
}   // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HPP_ */
