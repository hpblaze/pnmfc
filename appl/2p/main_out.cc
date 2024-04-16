// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the two-phase pore-network model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/porenetwork/2p/newtonsolver.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/porenetwork/dgfwriter.hh>
#include "problem_2p_im_magic_out.hh"

#include <dumux/common/initialize.hh>


int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::ImbibitionProblem;

    // maybe initialize MPI and/or multithreading backend, finalize is done automatically on exit
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = Dumux::PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();
    std::cout << "Grid initialized" << std::endl;

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    PoreNetwork::writeDgf("2p_im-grid.dgf", leafGridView, *gridData);

    // initialize marker vector for hydrophilic BC
    //std::vector<int> marker(leafGridView.size(1),0);

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);
    std::cout << "Finite Volume Grid created" << std::endl;

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);
    std::cout << "Spatial parameters set" << std::endl;

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);
    std::cout << "Problem set" << std::endl;

    // the solution vector
    using GridView = typename GridGeometry::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));
    problem->applyInitialSolution(x);
    problem->calculateSumSourceVolume();
    auto xOld = x;
    std::cout << "Solution vector set" << std::endl;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    std::cout << "Grid variables initialized" << std::endl;
    
    /*
    const auto poresToPlot = getParam<std::vector<std::size_t>>("PlotMaterialLaw.PoresToPlot", std::vector<std::size_t>{});
    const auto sRange = getParam<std::array<double, 2>>("PlotMaterialLaw.SaturationRange", std::array<double, 2>{0, 1});
    spatialParams->plotPcSw<typename GridVariables::VolumeVariables>(poresToPlot, *problem, sRange[0], sRange[1]);
    */
    
    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    /*
    // specify whether to drainage shall occur with step-wise applied pressures or using constant BCs
    const bool stepWiseDrainage = getParam<int>("Problem.NumSteps") > 1;
    */

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");
    std::cout << "set restart time" << std::endl;

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); //! Add model specific output fields
    std::cout << "Vtk output files initialized" << std::endl;

    vtkWriter.addField(gridGeometry->poreVolume(), "poreVolume", Vtk::FieldType::vertex);
    vtkWriter.addField(gridGeometry->throatShapeFactor(), "throatShapeFactor", Vtk::FieldType::element);
    vtkWriter.addField(gridGeometry->throatCrossSectionalArea(), "throatCrossSectionalArea", Vtk::FieldType::element);
    vtkWriter.write(0.0);

    Dumux::PoreNetwork::AveragedValues<GridVariables, SolutionVector> avgValues(*gridVariables, x);
    using FS = typename GridVariables::VolumeVariables::FluidSystem;
    avgValues.addAveragedQuantity([](const auto& v){ return v.saturation(FS::phase0Idx); }, [](const auto& v){ return v.poreVolume(); }, "avgSat"); // Function to extract saturation values
    avgValues.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase0Idx); }, [](const auto& v){ return v.saturation(FS::phase0Idx)*v.poreVolume(); }, "avgPw"); 
    avgValues.addAveragedQuantity([](const auto& v){ return v.pressure(FS::phase1Idx); }, [](const auto& v){ return v.saturation(FS::phase1Idx)*v.poreVolume(); }, "avgPn");
    std::vector<std::size_t> dofsToNeglect;

    for (const auto& vertex : vertices(leafGridView))
    {
        using Labels = GetPropType<TypeTag, Properties::Labels>;
        const auto vIdx = gridGeometry->vertexMapper().index(vertex);
        if (gridGeometry->poreLabel(vIdx) == Labels::inlet || gridGeometry->poreLabel(vIdx) == Labels::outlet)
            dofsToNeglect.push_back(vIdx);
    }

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    std::cout << "time loop init check" << std::endl;

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
    std::cout << "assembler check" << std::endl;

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    std::cout << "linear solver check" << std::endl;

    // the non-linear solver
    using NewtonSolver = PoreNetwork::TwoPNewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    std::cout << "non-linear solver check" << std::endl;

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // try solving the non-linear system
        nonLinearSolver.solve(x, *timeLoop);
        std::cout << "non-linear solver check time loop" << std::endl;

        // make the new solution the old solution
        xOld = x;

        // calculate the averaged values
        avgValues.eval();
        problem->postTimeStep(timeLoop->time(), avgValues, gridVariables->gridFluxVarsCache().invasionState().numThroatsInvaded(), timeLoop->timeStepSize());

        gridVariables->advanceTimeStep();

        //Update the outlet pores data with magic pores
        problem->ObtainMagicPoreData(x);

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if(problem->shouldWriteOutput(timeLoop->timeStepIndex(), *gridVariables))
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    problem->postTimeStep(timeLoop->time(), avgValues, gridVariables->gridFluxVarsCache().invasionState().numThroatsInvaded(), timeLoop->timeStepSize());

    nonLinearSolver.report();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
}
