cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(fibergen)

find_package(OpenMP)
find_package(VTK 9.3.0 COMPONENTS
  CommonColor
  CommonCore
  CommonDataModel
  FiltersCore
  FiltersSources
  IOImage
  IOGeometry
  ImagingHybrid
  InteractionStyle
  InteractionWidgets
  RenderingAnnotation
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "MarchingCubes: Unable to find the VTK build folder.")
endif()

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")
add_compile_options(-O3 -ffast-math -g)
add_executable(fibergen MACOSX_BUNDLE src/fibergen.cpp)
  target_link_libraries(fibergen PRIVATE ${VTK_LIBRARIES} OpenMP::OpenMP_CXX
)


include_directories(include)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS fibergen
  MODULES ${VTK_LIBRARIES}
)
