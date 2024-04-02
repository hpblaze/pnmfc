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
#ifndef DUMUX_PNM1P_PROBLEM_HH
#define DUMUX_PNM1P_PROBLEM_HH
// standard C++ library
#include <ctime> //can include time related stuff
#include <iostream> //can print on the terminal and to interact
// grid
#include <dune/foamgrid/foamgrid.hh>  //creates 1d grid for pore network
// spatial params
#include "files/spatialparams.hh" //We can change the properties of the pores and throats like contact angle and advection flux
// problem
#include <dumux/porousmediumflow/problem.hh> // base problem 
#include <dumux/porenetwork/1p/model.hh> // Pore network model
#include <dumux/common/boundarytypes.hh> //useful to specify boundary type
//properties
#include <dumux/common/properties.hh>
#include <dumux/material/components/simpleh2o.hh> // all required details of fluid component properties
#include <dumux/porenetwork/common/utilities.hh> // useful functions for pnm (e.g. for the calculation of fluxes at the boundary)
#include <dumux/material/fluidsystems/1pliquid.hh> //provides required fluid properties data for jacobian
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh> //Element-wise calculation of the residual and its derivatives for a single-phase, incompressible, test problem.

#ifndef ISOTHERMAL //helps to choose between isothermal and non-isothermal case
#define ISOTHERMAL 1
#endif

 // namespace is used to organize code into logical groups within a library makes it easy to manage.
namespace Dumux {      
// templates to write generic code that works with different data types without having to rewrite the same code for each type
template <class TypeTag>
class PNMOnePProblem;

namespace Properties{
//Create new type tags
namespace TTag{
#if ISOTHERMAL
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOneP>; }; //tuples provide a flexible and convenient way to work with collections of heterogeneous elements and to pass around groups of values in a concise manner.
#else
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOnePNI>; };
#endif
} 

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePProblem> { using type = PNMOnePProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePProblem>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePProblem> { using type = Dune::FoamGrid<1, 3>; };

// Set the advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type =  PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMOnePProblem> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

} //end namespace Properties

template <class TypeTag>
class PNMOnePProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
    

    };

    
// defining get params from input file
public:
    template<class SpatialParams>
    PNMOnePProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        verbose_ = getParam<bool>("Problem.Verbose", false);
        VtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        testHydrostaticPressure_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableGravity");
        inletPressure_ = getParam<Scalar>("Problem.InletPressure");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure");
        sourceFluxH2O_ = getParam<Scalar>("Problem.sourceFluxH2O");
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
        if (isOutletPore_(scv))
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
        /*
        if(isInletPore_(scv)){
            values[Indices::pressureIdx] = inletPressure_;
        }
        
        if(isSourcePore_(scv)){
            values[Indices::pressureIdx] = inletPressure_;
        }
        */
        if(isOutletPore_(scv)){
            values[Indices::pressureIdx] = outletPressure_;
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
            values[Indices::pressureIdx] = 1e6;
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
        if (testHydrostaticPressure_)
            return false;

        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }


    bool verbose_;
    int VtpOutputFrequency_;
    bool testHydrostaticPressure_;
    Scalar inletPressure_;
    Scalar outletPressure_;
    Scalar sourceFluxH2O_;
};
} //end namespace Dumux

#endif
