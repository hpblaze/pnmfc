// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A problem file for the one-phase pore network model.
 */
#ifndef DUMUX_PNM1P2C_PROBLEM_HH
#define DUMUX_PNM1P2C_PROBLEM_HH
// standard C++ library
#include <ctime> //can include time related stuff
#include <iostream> //can print on the terminal and to interact
// grid
#include <dune/foamgrid/foamgrid.hh>  //creates 1d grid for pore network
// spatial params
#include "files/spatialparams.hh" //We can change the properties of the pores and throats like contact angle and advection flux
// problem
#include <dumux/porousmediumflow/problem.hh> // base problem 
#include <dumux/porenetwork/1pnc/model.hh> // Pore network model
#include <dumux/common/boundarytypes.hh> //useful to specify boundary type
//properties
#include <dumux/common/properties.hh>
#include <dumux/material/fluidsystems/1padapter.hh> //provides all fluid properties details for human readable
#include <dumux/material/fluidsystems/h2oair.hh> //provides required fluid properties data for jacobian

#ifndef ISOTHERMAL //helps to choose between isothermal and non-isothermal case
#define ISOTHERMAL 1
#endif

 // namespace is used to organize code into logical groups within a library makes it easy to manage.
namespace Dumux {      
// templates to write generic code that works with different data types without having to rewrite the same code for each type
template <class TypeTag>
class PNMOnePTwoCProblem;

namespace Properties{
//Create new type tags
namespace TTag{
#if ISOTHERMAL
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNC>; }; //tuples provide a flexible and convenient way to work with collections of heterogeneous elements and to pass around groups of values in a concise manner.
#else
struct PNMOnePTwoCProblem { using InheritsFrom = std::tuple<PNMOnePNCNI>; };
#endif
} //end namespace type tag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePTwoCProblem> { using type = PNMOnePTwoCProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePTwoCProblem>
{
    //using Policy = FluidSystems::H2OAirDefaultPolicy</*simple*/true>;
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    using type = FluidSystems::OnePAdapter<H2OAir>;
};

// Set the Spatial Params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOnePTwoCProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::CompositionalSpatialParams<GridGeometry, Scalar>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePTwoCProblem> { using type = Dune::FoamGrid<1, 3>; };

// Set the advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePTwoCProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type =  PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

} //end namespace Properties

template <class TypeTag>
class PNMOnePTwoCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    //using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    //using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    //using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    //using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
        transportEqIdx = 1,
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
    
    };

    
// defining get params from input file
public:
    template<class SpatialParams>
    PNMOnePTwoCProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        //get some parameters from the input file
        verbose_ = getParam<bool>("Problem.Verbose", false);
        VtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        inletPressure_ = getParam<Scalar>("Problem.InletPressure");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure");
        inletMoleFraction_ = getParam<Scalar>("Problem.InletMoleFraction");
        outletMoleFraction_ = getParam<Scalar>("Problem.OutletMoleFraction");
        sourceFluxH2O_ = getParam<Scalar>("Problem.sourceFluxH2O");

        FluidSystem::init();
    }

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const //define output
    {
        if(VtpOutputFrequency_ < 0)
            return true;
        if(VtpOutputFrequency_ == 0)
            return ( timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return ( timeStepIndex % VtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
    }

#if ISOTHERMAL
    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 80; } // 10°C

#endif

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv) || isOutletPore_(scv))
            bcTypes.setAllDirichlet();
        else // neuman for the remaining boundaries
            bcTypes.setAllNeumann();
#if !ISOTHERMAL
        bcTypes.setDirichlet(Indices::temperatureIdx);
#endif
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        control volume.
     *
     * \param values The Dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        
        if(isInletPore_(scv)){
            values[Indices::pressureIdx] = inletPressure_;
            values[transportEqIdx] = inletMoleFraction_;
            std::cout << "Actual.AirFlux value: " << values[transportEqIdx] << "moles"<< std::endl;

        }
        /*
        if(isSourcePore_(scv)){
            values[Indices::pressureIdx] = inletPressure_;
        }
        */
        if(isOutletPore_(scv)){
            values[Indices::pressureIdx] = outletPressure_;
            values[transportEqIdx] = outletMoleFraction_;
            
        }       
#if !ISOTHERMAL
         if(isInletPore_(scv))
            values[Indices::temperatureIdx] = 273.15 +25;
         else
            values[Indices::temperatureIdx] = 273.15 +20;
#endif
         return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        //return PrimaryVariables(0);
        PrimaryVariables values(0.0);
        if(isSourcePore_(scv)){
            std::cout << "Problem.SourceFluxH2O value: " << sourceFluxH2O_/(this->gridGeometry().poreVolume(scv.dofIndex())) << "[kg/(m^3.s)]"<< std::endl;
            values[Indices::conti0EqIdx] += sourceFluxH2O_/(this->gridGeometry().poreVolume(scv.dofIndex()));
            std::cout << "Actual.sourceFluxH2O value: " << values[Indices::conti0EqIdx] << "[kg/(m^3.s)]"<< std::endl;
            std::cout << "Pore volume: " << this->gridGeometry().poreVolume(scv.dofIndex()) << std::endl;
        }
        return values;
    }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
            values[Indices::pressureIdx] = 1e5;
#if !ISOTHERMAL
            values[Indices::temperatureIdx] = 273.15 +20;
#endif
            return values;
    }

    const bool verbose() const
    { return verbose_; }

    bool simulationFinished() const
    { return true ; }


private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::inlet;
    }

    bool isSourcePore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::source;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }


    bool verbose_;
    int VtpOutputFrequency_;
    Scalar inletPressure_;
    Scalar outletPressure_;
    Scalar inletMoleFraction_;
    Scalar outletMoleFraction_;
    Scalar sourceFluxH2O_;
};
} //end namespace Dumux

#endif
