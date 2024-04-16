// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the two-phase pore network model (Drainage).
 */
#ifndef DUMUX_Drainage_PROBLEM_HH
#define DUMUX_Drainage_PROBLEM_HH

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
#include "files/spatialparams_dr.hh" // spatial params


namespace Dumux 
{
template <class TypeTag>
class DrainageProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
}  // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = DrainageProblem<TypeTag>; };

// Set the fluid system: H2O and Air as phases
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };

// Set formulation (pw and sn or pn and sw)
template<class TypeTag>
struct Formulation<TypeTag, TTag::DrainageProblem>
{ static constexpr auto value = TwoPFormulation::p0s1; }; //p1s0 is pnsw

// Set the spatial params  
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DrainageProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using LocalRules = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<Scalar>;
public:
    using type = PoreNetwork::TwoPDrainageSpatialParams<GridGeometry, Scalar, LocalRules>;
};

} //end namespace Dumux::Properties

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
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
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::saturationIdx,

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

// defining get params from input file
public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        //get some parameters from the input file
        verbose_ = getParam<bool>("Problem.Verbose", true);
        VtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        sink_k_ = getParam<Scalar>("Problem.sink_k");
        pnInlet_ = getParam<Scalar>("Problem.pnInlet");
        swInlet_ = getParam<Scalar>("Problem.swInlet");
        sigma_ = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725);
        logfile_.open("time_steps_" + this->name() + ".txt");


        //get Indices from components and eq's --> output on console
        std::cout << "---------------- START INDICES OUTPUT ----------------" << std::endl;
        std::cout << "wetting phase index:" << wPhaseIdx << std::endl;
        std::cout << "non-wetting phase index:" << nPhaseIdx << std::endl;
        std::cout << "pressure index:" << pnIdx << std::endl;
        std::cout << "saturation index:" << swIdx << std::endl;
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

        return bcTypes;
    }
   
    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);

        if (isInletPore_(scv))
        {
            values[pwIdx] = pwInlet_;
            values[snIdx] = snInlet_;
        }

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
        // We fix the mass flux of non-wetting injection at inlet of pore-network
        // The total inlet mass flux is distributed according to the ratio
        // of pore volume at inlet to the total volumess
        if (isSourcePore_(scv))
        {   

            //std::cout << "water saturation in sink pore: " << elemVolVars[scv].saturation(wPhaseIdx) << std::endl;
            //std::cout << "Air saturation in source pore: " << elemVolVars[scv].saturation(nPhaseIdx) << std::endl;
            using MaterialLaw = typename GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw;
            values[wPhaseIdx] -= 0.001*(/this->gridGeometry().poreVolume(scv.dofIndex()));; //sink_k_ * (1-(numThroatsInvaded/this->gridGeometry().gridView().size(0))); //static_cast<Scalar>(this->gridGeometry().size(0))
            std::cout << "Actual.sinkNonWetFluxAir value: " << values[wPhaseIdx] << std::endl;
            //std::cout << "invaded pores: " << numThroatsInvaded << std::endl;
            //std::cout << "total throats: " << this->gridGeometry().gridView().size(0) << std::endl;
            //std::cout << "scv.volume(): " << scv.volume() << std::endl;
        }
        
        return values;
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        /*
        values[pnIdx] = pnOutlet_;

        //get global index of pore
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if(isInletPore_(dofIdxGlobal))
            values[swIdx] = swInlet_;
	    else if (isSourcePore_(dofIdxGlobal))
	        values[swIdx] = 1.0;
        else
            values[swIdx] = swInlet_;
        return values;
        */
        values[pnIdx] = 1e4;
        values[swIdx] = 1.0;

        return values;

    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    { return false; }

    // \}

    //! Loop over the scv in the domain to calculate the sum volume of inner inlet pores
    void calculateSumSourceVolume()
    {
        sumSourcePoresVolume_ = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (isSourcePore_(scv))
                {
                    sumSourcePoresVolume_ += this->gridGeometry().poreVolume(scv.dofIndex())/this->gridGeometry().coordinationNumber(scv.dofIndex());
                }
            }
        }
    }

    /*!
     * \brief Called at the end of each time step
     */
    template<class AveragedValues>
    void postTimeStep(const Scalar time, const AveragedValues& avgValues, std::size_t numThroatsInvaded, const Scalar dt)
    {
        //const Scalar avgSw = avgValues["avgSat"];

        std::cout << "invaded throats: " << numThroatsInvaded << std::endl;
        logfile_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                 << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                 << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                 << std::endl;
    }
    
private:

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
    Scalar sink_k_;
    Scalar pnInlet_;
    Scalar swInlet_;
    Scalar sigma_;
    Scalar sumSourcePoresVolume_;
    Scalar sumInletPoresVolume_;
    std::ofstream logfile_;
    std::size_t numThroatsInvaded;
};
} //end namespace Dumux

#endif
