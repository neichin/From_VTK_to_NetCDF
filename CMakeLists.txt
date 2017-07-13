
cmake_minimum_required(VERSION 3.5.0)

PROJECT(fromVTKtoElmer)

set(VTK_DIR "/home/imerino/CODES/VTK/VTK-7.1.1/BUILD")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(fromVTKtoElmer MACOSX_BUNDLE fromVTKtoElmer)

include_directories("/home/imerino/CODES/netcdf-cxx-4.2/BUILD/include/")
#link_directories("/home/imerino/CODES/netcdf-cxx-4.2/BUILD/lib/")

#link_directories("/../../../CODES/netcdf-cxx-4.2/BUILD/lib/")

find_library(NETCDF_LIB NAMES netcdf_c++ PATHS /home/imerino/CODES/netcdf-cxx-4.2/BUILD/lib/)

target_link_libraries(fromVTKtoElmer ${NETCDF_LIB})

if(VTK_LIBRARIES)
  target_link_libraries(fromVTKtoElmer ${VTK_LIBRARIES})
else()
  target_link_libraries(fromVTKtoElmer vtkHybrid vtkWidgets)
endif()
