#ifndef MORIS_CL_GEN_PDV_HPP_
#define MORIS_CL_GEN_PDV_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv
        {
        public:
            bool mIsActive = true;

            moris_id mId = gNoID;

        protected:

            /**
             * constructor
             */
            Pdv();

        public:

            /**
             * trivial destructor
             */
            ~Pdv();

            void set_id( const moris_id & tId)
            {
                MORIS_ASSERT( mId == gNoID, "Pdv::set_id(), id of this pdv was set before" );
                mId = tId;
            };

            moris_id get_id()
            {
                return mId;
            };

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            virtual real get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates) = 0;

            virtual void set_value(const moris::real & aValue )
            {
                MORIS_ERROR( false, "Pdv::set_value(), not implemented for this pdv type" );
            }

            /**
             * Get the PDV sensitivity with respect to ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @param aSensitivities Matrix of sensitivities to be returned
             */
            virtual void get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities) = 0;

        };
    }   // end ge namespace
}   // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HPP_ */
