# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(static)

add_input_file_links()
dune_symlink_to_source_files(FILES grids)
dune_symlink_to_source_files(FILES extractresults.py plot_constitutive_laws.py)

dumux_add_test(NAME 2p_im
               SOURCES main.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME 2p_im_bc
               SOURCES main_bc.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME 2p_im_bc_sink
               SOURCES main_sink.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME 2p_im_magic_out
               SOURCES main_out.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME 2p_dr_bc_sink
               SOURCES main_dr_sink.cc
               LABELS porenetwork
               COMPILE_DEFINITIONS ISOTHERMAL=1
               CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
               TIMEOUT "1000"
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p  params.input")

dumux_add_test(NAME plot_constitutive_relations
               SOURCES main_plot_constitutive_laws.cc
               COMPILE_DEFINITIONS ISOTHERMAL=1
               COMPILE_DEFINITIONS ENABLECACHE=0
               COMPILE_DEFINITIONS PLOTCONSTITUTIVERELATIONS=1)

