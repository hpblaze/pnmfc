// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief Spatial parameters for a pore-network model (properties of the porous matrix)
 * 
 */
#ifndef DUMUX_PNM_CREEPING_SPATIAL_PARAMS_1P_HH
#define DUMUX_PNM_CREEPING_SPATIAL_PARAMS_1P_HH

#include <dumux/porenetwork/common/spatialparams.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar>
class CreepingSpatialParams : public SpatialParams<GridGeometry, Scalar, CreepingSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = SpatialParams<GridGeometry, Scalar, CreepingSpatialParams<GridGeometry, Scalar>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    template<class ElementSolutionVector>
    auto temperature(const Element& element,
                     const SubControlVolume& scv,
                     const ElementSolutionVector& elemSol) const
    { return 273.15 + 10.0; }

    template<class ElementSolutionVector>
    auto poreShapeFactor(const Element& element,
                         const SubControlVolume& scv,
                         const ElementSolutionVector& elemSol) const
    {
        return Throat::shapeFactorCircle<Scalar>();
    }

    template<class ElementSolutionVector>
    Scalar poreCrossSectionalArea(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        const auto poreRadius = this->gridGeometry().poreInscribedRadius(scv.dofIndex());
        return poreRadius*poreRadius*M_PI;
    }

    template<class ElementSolutionVector>
    Scalar poreLength(const Element& element,
                      const SubControlVolume& scv,
                      const ElementSolutionVector& elemSol) const
    {
        return this->gridGeometry().poreInscribedRadius(scv.dofIndex());
    }

};
} // end namespace Dumux::PoreNetwork

#endif
