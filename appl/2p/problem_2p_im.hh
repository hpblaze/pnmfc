// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the two-phase pore network model.
 */
#ifndef DUMUX_IMBIBITION_PROBLEM_HH
#define DUMUX_IMBIBITION_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/2p/model.hh>
#include <dumux/porenetwork/2p/spatialparams.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include <dumux/common/properties.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/porenetwork/common/utilities.hh>

#include "files/h2oair.hh"
#include "files/spatialparams_im.hh" // spatial params

#ifndef ISOTHERMAL
#define ISOTHERMAL 1
#endif


namespace Dumux 
{
template <class TypeTag>
class ImbibitionProblem;

namespace Properties {

// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct ImbibitionProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
#else
struct ImbibitionProblem { using InheritsFrom = std::tuple<PNMTwoPNI>; };
#endif
}  // end namespace TTag

// Set formulation (pw and sn or pn and sw)
template<class TypeTag>
struct Formulation<TypeTag, TTag::ImbibitionProblem>
{ static constexpr auto value = TwoPFormulation::p1s0; }; //p1s0 is pnsw

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ImbibitionProblem> { using type = ImbibitionProblem<TypeTag>; };

// Set the fluid system: H2O and Air as phases
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ImbibitionProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ImbibitionProblem> { using type = Dune::FoamGrid<1, 3>; };

// Set the spatial params  
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ImbibitionProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPImbibitionSpatialParams<GridGeometry, Scalar, LocalRules>;
};

//decision to plot
#if PLOTCONSTITUTIVERELATIONS
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ImbibitionProblem> { static constexpr bool value = !PLOTCONSTITUTIVERELATIONS; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ImbibitionProblem> { static constexpr bool value = !PLOTCONSTITUTIVERELATIONS; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ImbibitionProblem> { static constexpr bool value = !PLOTCONSTITUTIVERELATIONS; };
#endif

} //end namespace Dumux::Properties

template <class TypeTag>
class ImbibitionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // phase indices
        wPhaseIdx = FluidSystem::phase0Idx, //H2O
        nPhaseIdx = FluidSystem::phase1Idx, //Air

        // primary variable indices
        pnIdx = Indices::pressureIdx,
        swIdx = Indices::saturationIdx,


#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

// defining get params from input file
public:
    template<class SpatialParams>
    ImbibitionProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        //get some parameters from the input file
        verbose_ = getParam<bool>("Problem.Verbose", true);
        VtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        distributeByVolume_ = getParam<bool>("Problem.DistributeByVolume", true);
        pnOutlet_ = getParam<Scalar>("Problem.pnOutlet");
        swOutlet_ = getParam<Scalar>("Problem.swOutlet");
        sigma_ = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725);
        sourceWetFluxH2O_ = getParam<Scalar>("Problem.sourceWetFluxH2O");
        sourceNonWetFluxAir_ = getParam<Scalar>("Problem.sourceNonWetFluxAir", 0);
        logfile_.open("time_steps_" + this->name() + ".txt");


        //get Indices from components and eq's --> output on console
        std::cout << "---------------- START INDICES OUTPUT ----------------" << std::endl;
        std::cout << "wetting phase index:" << wPhaseIdx << std::endl;
        std::cout << "non-wetting phase index:" << nPhaseIdx << std::endl;
        std::cout << "pressure index:" << pnIdx << std::endl;
        std::cout << "saturation index:" << swIdx << std::endl;
        //std::cout << "Temperature index:" << temperatureIdx << std::endl;
        //std::cout << "energy eq index:" << energyEqIdx << std::endl;
        std::cout << "---------------- END INDICES OUTPUT ----------------"<< std::endl;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        if (VtpOutputFrequency_ < 0)
            return true;

        if (VtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return (timeStepIndex % VtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
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
        /*
        if (onLeftBoundary_(scv) || onRightBoundary_(scv) || onBottomBoundary_(scv))
        {
            bcTypes.setAllNeumann();
            //bcTypes.setNeumann(pnIdx); //set pressure constant
            //bcTypes.setNeumann(swIdx); //set saturation constant
        }
        else 
        */
        if (isInletPore_(scv) || isOutletPore_(scv))
            bcTypes.setAllDirichlet();
        else // neuman for the remaining boundaries
            bcTypes.setAllNeumann();

#if !ISOTHERMAL
        bcTypes.setDirichlet(temperatureIdx);
#endif
        return bcTypes;
    }
   
    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        
        if (isInletPore_(scv))
        {
            values[pnIdx] = 1.0e5;
            values[swIdx] = 0.5;
        }
        else if (isOutletPore_(scv))
        {
            values[pnIdx] = pnOutlet_;
            values[swIdx] = swOutlet_;
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

    //! Evaluate the source term for all phases within a given sub-control-volume.
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);

        // The total inlet mass flux is distributed according to the ratio of pore volume
        if(isSourcePore_(scv)){
            std::cout << "Problem.SourceWetFluxH2O value: " << sourceWetFluxH2O_/(this->gridGeometry().poreVolume(scv.dofIndex())) << "[kg/(m^3.s)]"<< std::endl;
            values[nPhaseIdx] -= sourceWetFluxH2O_/(this->gridGeometry().poreVolume(scv.dofIndex()));
            std::cout << "Actual.sourceWetFluxH2O value: " << values[nPhaseIdx] << "[kg/(m^3.s)]"<< std::endl;
            //std::cout << "Pore volume: " << this->gridGeometry().poreVolume(scv.dofIndex()) << std::endl;
        }
        return values;
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pnIdx] = 1e5;
        values[swIdx] = 0.0;

        // get global index of pore
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if (isInletPore_(dofIdxGlobal))
        {
            values[pnIdx] = 1e5;
            values[swIdx] = 0.0;
        }

#if !ISOTHERMAL
        values[temperatureIdx] = 273.15 + 10;
#endif
        return values;
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    { return true; }

    // \}

    //! Loop over the scv in the domain to calculate the sum volume of inner inlet pores
    void calculateSumInletVolume()
    {
        sumInletPoresVolume_ = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (isInletPore_(scv))
                {
                    if (distributeByVolume_)
                        sumInletPoresVolume_ += this->gridGeometry().poreVolume(scv.dofIndex())/this->gridGeometry().coordinationNumber(scv.dofIndex());
                    else
                        sumInletPoresVolume_ += std::pow(this->gridGeometry().poreInscribedRadius(scv.dofIndex()), 4.0);
                }
            }
        }
    }

    /*!
     * \brief Called at the end of each time step
     */
    void postTimeStep(const Scalar time)
    {
        logfile_ << std::fixed << std::left << std::setw(20)
                 << std::setfill(' ') << time << std::endl;
    }
    
private:

    bool onLeftBoundary_(const SubControlVolume& scv) const
    {
        if(scv.dofPosition()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            return true;
        else
            return false;

    }

    bool onRightBoundary_(const SubControlVolume& scv) const
    {
        if(scv.dofPosition()[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            return true;
        else
            return false;

    }

    bool onBottomBoundary_(const SubControlVolume& scv) const
    {
        if(scv.dofPosition()[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            return true;
        else
            return false;
    }

    bool onTopBoundary_(const SubControlVolume& scv) const
    {
        if(scv.dofPosition()[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            return true;
        else
            return false;
    }

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    bool isSourcePore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) ==Labels::source;
    }
    bool isSourcePore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) ==Labels::source;
    }

    static constexpr Scalar eps_ = 1e-6; 
    int VtpOutputFrequency_;
    bool verbose_;
    bool useFixedPressureAndSaturationBoundary_;
    bool distributeByVolume_;
    Scalar pnOutlet_;
    Scalar swOutlet_;
    Scalar sigma_;
    Scalar sourceWetFluxH2O_;
    Scalar sourceNonWetFluxAir_;
    Scalar sumInletPoresVolume_;
    std::ofstream logfile_;
};
} //end namespace Dumux

#endif
