# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_subdirectory(nonisothermal)
add_subdirectory(noncreepingflow)

add_input_file_links()
dune_symlink_to_source_files(FILES grids)

add_executable(pnm_1p EXCLUDE_FROM_ALL main.cc)

#set(CMAKE_BUILD_TYPE Debug)