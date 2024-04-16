// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the two-phase pore network model (imbibition).
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


namespace Dumux 
{
template <class TypeTag>
class ImbibitionProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct ImbibitionProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
}  // end namespace TTag

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

// Set formulation (pw and sn or pn and sw)
template<class TypeTag>
struct Formulation<TypeTag, TTag::ImbibitionProblem>
{ static constexpr auto value = TwoPFormulation::p1s0; }; //p1s0 is pnsw

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

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

// defining get params from input file
public:
    template<class SpatialParams>
    ImbibitionProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        //get some parameters from the input file
        verbose_ = getParam<bool>("Problem.Verbose", true);
        VtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        sourceWetFluxH2O_ = getParam<Scalar>("Problem.sourceWetFluxH2O");
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

        if (isOutletPore_(scv)) //|| isSourcePore_(scv)
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
        
        if (isOutletPore_(scv))
        {
            values[pnIdx] = 1e5;
            if (swOutlet_ < 1 && count_ == 34)
            {
                values[swIdx] = swOutlet_/9;
                std::cout << "Updated swIdx: " << swOutlet_ << std::endl;
                std::cout << "Updated sumSourcePoresVolume: " << sumSourcePoresVolume_ << std::endl;
                std::cout << "count: " << count_ << std::endl;
            }
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
            values[wPhaseIdx] = sourceWetFluxH2O_ * ( scv.volume() / sumSourcePoresVolume_ );
            //std::cout << "Actual.sourceWetFluxH2O value: " << values[swIdx] << std::endl;
            //std::cout << "scv.volume(): " << scv.volume() << std::endl;
        }

        return values / scv.volume();
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pnIdx] = 1e5;
        values[swIdx] = 0.0;
        /*
        //get global index of pore
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);

        if (isSourcePore_(dofIdxGlobal))
        {
            values[snIdx] = 0;
        }
        else
        */
            

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

    //! Loop over the scv in the outlet to update the saturation with magic pore data
    template<class Solvector>
    void ObtainMagicPoreData(const Solvector& x)
    {
        //Obtaining Wet-Sat from magic pores from previous timestep
        //Scalar swOutlet_(0.0);
        swOutlet_ = 0.0;
        count_ = 0;
        auto fvGeometry = localView(this->gridGeometry());
        for(const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {   
                const auto dofIdx = scv.dofIndex();
                if (isInletPore_(scv))
                    {
                        swOutlet_ += x[dofIdx][swIdx];
                        count_++;
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
    Scalar pnOutlet_;
    Scalar swOutlet_;
    Scalar sigma_;
    Scalar sourceWetFluxH2O_;
    Scalar sumSourcePoresVolume_;
    Scalar AvgNonWetSat_;
    Scalar count_;
    std::ofstream logfile_;
};
} //end namespace Dumux

#endif
