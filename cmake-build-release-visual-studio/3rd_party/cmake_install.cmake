# Install script for directory: D:/Q4/photogrammetry/assignment2/A2_Triangulation/A2_Triangulation_Code/3rd_party

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/A2_Triangulation_Code")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/Q4/photogrammetry/assignment2/A2_Triangulation/A2_Triangulation_Code/cmake-build-release-visual-studio/3rd_party/glew/cmake_install.cmake")
  include("D:/Q4/photogrammetry/assignment2/A2_Triangulation/A2_Triangulation_Code/cmake-build-release-visual-studio/3rd_party/glfw/cmake_install.cmake")
  include("D:/Q4/photogrammetry/assignment2/A2_Triangulation/A2_Triangulation_Code/cmake-build-release-visual-studio/3rd_party/glog/cmake_install.cmake")
  include("D:/Q4/photogrammetry/assignment2/A2_Triangulation/A2_Triangulation_Code/cmake-build-release-visual-studio/3rd_party/cminpack/cmake_install.cmake")

endif()

