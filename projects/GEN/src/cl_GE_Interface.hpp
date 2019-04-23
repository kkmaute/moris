/*
 * cl_GE_Interface.hpp
 *
 *  Created on: Apr 4, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_

// GE includes
//------------------------------------------------------------
// MTK includes
#include "cl_MTK_Mesh.hpp"
//------------------------------------------------------------


namespace moris
{
namespace ge
{

class Geometry;

//------------------------------------------------------------
class Geometry_Engine_Interface
{
public:
    // constuctor
    Geometry_Engine_Interface(){};

    // destructor
    ~Geometry_Engine_Interface(){};
    //------------------------------------------------------------
    virtual void
    set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
                                                      real   aThreshold   = 0.0)
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::set_geometry: not set. ");
    }

    //------------------------------------------------------------
    virtual moris::Cell< uint >
    flag_element_list_for_refinement( moris::Cell< mtk::Cell* > & aElementList ,
                                      moris::Cell< real >       & aConstant,
                                      uint                      aWhichGeometry )
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::flag_element_list_for_refinement: not implemented. ");
        moris::Cell< uint > tReturn(1);
        tReturn(0) = 1.0;

        return tReturn;
    }
    //------------------------------------------------------------
    virtual moris::Cell< uint >
    check_for_intersection( moris::Cell< real > &   aConstant,
                            mtk::Mesh*          &   aMeshPointer,
                            uint                    aWhichGeometry)
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::check_for_intersection: not implemented. ");
        moris::Cell< uint > tReturn(1);
        tReturn(0) = 1.0;

        return tReturn;
    }
    //------------------------------------------------------------
    virtual Matrix< DDRMat >
    get_edge_normal_for_straight_edge_quad4(    uint const & aElemGlobInd,
                                                uint const & aEdgeSideOrd,
                                                mtk::Mesh* & aMeshPointer   )
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::get_edge_normal_for_straight_edge_quad4: not implemented. ");
        Matrix< DDRMat > tNormal;
        tNormal(0,0) = 1.0;
        return tNormal;
    }
    //------------------------------------------------------------
    virtual void
    find_cells_within_levelset(
                  Cell< mtk::Cell * >      & aCells,
                  Cell< mtk::Cell * >      & aCandidates,
                  const  Matrix< DDRMat >  & aVertexValues,
                  const  uint                aUpperBound = 0.0 )
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::find_cells_within_levelset: not implemented. ");
    }
    //------------------------------------------------------------
    virtual void
    find_cells_intersected_by_levelset(
            Cell< mtk::Cell * >         & aCells,
            Cell< mtk::Cell * >         & aCandidates,
            const  Matrix< DDRMat >     & aVertexValues,
            const  real                   aLowerBound = -0.0001,
            const  real                   aUpperBound =  0.0001)
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::find_cells_within_levelset: not implemented. ");
    }
    //------------------------------------------------------------
    virtual Matrix< DDRMat >
    determine_phi_hat( Matrix< DDRMat >      & aTMatrix,
                       mtk::Mesh*            & aMeshPointer,
                       moris::Cell< real >   & aGeomParams,
                       uint                  aWhichGeom )
    {
        MORIS_ERROR( false, "Geometry_Engine_Interface::get_edge_normal_for_straight_edge_quad4: not implemented. ");
        Matrix< DDRMat > tPhiHat;
        tPhiHat(0,0) = 1.0;
        return tPhiHat;
    }
    //------------------------------------------------------------



    //------------------------------------------------------------
private:


    //------------------------------------------------------------

protected:

    //------------------------------------------------------------
};
//------------------------------------------------------------
} // end ge namespace
} // end moris namespace

#endif /* PROJECTS_GEN_SRC_CL_GE_INTERFACE_HPP_ */
