# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()
dune_symlink_to_source_files(FILES grids py_plot)

dumux_add_test(NAME 2p_magic_het_reg
               SOURCES main_het.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=0
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME 2p_magic_het_noreg
               SOURCES main_het.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS REGULARIZATIONWITHPRESSURE=0
               COMPILE_DEFINITIONS ALLOWBREAKTHROUGH=1
               COMPILE_DEFINITIONS NOREGULARIZTAION=1
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")




