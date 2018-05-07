CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Name of the project
PROJECT(FEM)

# FOR LAB MACHINES ONLY!
# DO NOT EDIT
 SET(DEF_DIR_GLM "D:\\c++\\glm-0.9.8.4")
 SET(DEF_DIR_GLFW "D:\\c++\\glfw-3.2.1")
 SET(DEF_DIR_GLEW "D:\\c++\\glew-2.1.0")
 SET(DEF_DIR_EIGEN "D:\\c++\\eigen-3.3.4")
 SET(DEF_DIR_MOSEK "C:\\Program Files\\Mosek\\8\\tools\\platform\\win64x86")
 SET(DEF_DIR_IGL "D:\\c++\\libigl-1.2.1\\include")
#SET(DEF_DIR_GLM "C:\\c++\\glm-0.9.8.5")
#SET(DEF_DIR_GLFW "C:\\c++\\glfw-3.2.1")
#SET(DEF_DIR_GLEW "C:\\c++\\glew-2.1.0")
#SET(DEF_DIR_EIGEN "C:\\c++\\eigen-3.3.4")
#SET(DEF_DIR_MOSEK "C:\\Program Files\\Mosek\\8\\tools\\platform\\win64x86")
#SET(DEF_DIR_IGL "C:\\c++\\libigl-1.2.1\\include")

# Is this the solution?
# Override with `cmake -DSOL=ON ..`
OPTION(SOL "Solution" OFF)

# Use glob to get the list of all source files.
# We don't really need to include header and resource files to build, but it's
# nice to have them also show up in IDEs.
IF(${SOL})
  FILE(GLOB_RECURSE SOURCES "src0/*.cpp" "src0/*.c")
  FILE(GLOB_RECURSE HEADERS "src0/*.h" "src0/*.hpp")
ELSE()
  FILE(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.c")
  FILE(GLOB_RECURSE HEADERS "src/*.h" "src/*.hpp")
ENDIF()
FILE(GLOB_RECURSE GLSL "resources/*.glsl")

# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

# Get the GLM environment variable. Since GLM is a header-only library, we
# just need to add it to the include directory.
SET(GLM_INCLUDE_DIR "$ENV{GLM_INCLUDE_DIR}")
IF(NOT GLM_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLM_INCLUDE_DIR to the root directory of your GLM installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLM in ${DEF_DIR_GLM}")
    IF(IS_DIRECTORY ${DEF_DIR_GLM})
      MESSAGE(STATUS "Found!")
      SET(GLM_INCLUDE_DIR ${DEF_DIR_GLM})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLM_INCLUDE_DIR})

# Get the GLFW environment variable. There should be a CMakeLists.txt in the 
# specified directory.
SET(GLFW_DIR "$ENV{GLFW_DIR}")
IF(NOT GLFW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLFW_DIR to the root directory of your GLFW installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLFW in ${DEF_DIR_GLFW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLFW})
      MESSAGE(STATUS "Found!")
      SET(GLFW_DIR ${DEF_DIR_GLFW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
OPTION(GLFW_BUILD_EXAMPLES "GLFW_BUILD_EXAMPLES" OFF)
OPTION(GLFW_BUILD_TESTS "GLFW_BUILD_TESTS" OFF)
OPTION(GLFW_BUILD_DOCS "GLFW_BUILD_DOCS" OFF)
IF(CMAKE_BUILD_TYPE MATCHES Release)
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/release)
ELSE()
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/debug)
ENDIF()
INCLUDE_DIRECTORIES(${GLFW_DIR}/include)
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} glfw ${GLFW_LIBRARIES})

# Get the GLEW environment variable.
SET(GLEW_DIR "$ENV{GLEW_DIR}")
IF(NOT GLEW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLEW_DIR to the root directory of your GLEW installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLEW in ${DEF_DIR_GLEW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLEW})
      MESSAGE(STATUS "Found!")
      SET(GLEW_DIR ${DEF_DIR_GLEW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLEW_DIR}/include)
IF(WIN32)
  # With prebuilt binaries
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_DIR}/lib/Release/x64/glew32s.lib)
ELSE()
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_DIR}/lib/libGLEW.a)
ENDIF()

# Get the EIGEN environment variable. Since EIGEN is a header-only library, we
# just need to add it to the include directory.
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}")
IF(NOT EIGEN3_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable EIGEN3_INCLUDE_DIR to the root directory of your EIGEN installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for EIGEN in ${DEF_DIR_EIGEN}")
    IF(IS_DIRECTORY ${DEF_DIR_EIGEN})
      MESSAGE(STATUS "Found!")
      SET(EIGEN3_INCLUDE_DIR ${DEF_DIR_EIGEN})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# Get the MOSEK environment variable.
SET(MOSEK_DIR "$ENV{MOSEK_DIR}")
IF(NOT MOSEK_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLEW_DIR to the root directory of your GLEW installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLEW in ${DEF_DIR_MOSEK}")
    IF(IS_DIRECTORY ${DEF_DIR_MOSEK})
      MESSAGE(STATUS "Found!")
      SET(MOSEK_DIR ${DEF_DIR_MOSEK})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${MOSEK_DIR}/h)
IF(WIN32)
  # With prebuilt binaries
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MOSEK_DIR}/bin/mosek64_8_1.lib)
ENDIF()

# Get the IGL environment variable.
SET(IGL_INCLUDE_DIR "$ENV{IGL_INCLUDE_DIR}")
IF(NOT IGL_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable IGL_INCLUDE_DIR to the root directory of your IGL installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for IGL in ${DEF_DIR_IGL}")
    IF(IS_DIRECTORY ${DEF_DIR_IGL})
      MESSAGE(STATUS "Found!")
      SET(IGL_INCLUDE_DIR ${DEF_DIR_IGL})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${IGL_INCLUDE_DIR})

# OS specific options and libraries
IF(WIN32)
  # c++11 is enabled by default.
  # -Wall produces way too many warnings.
  # -pedantic is not supported.
  # Disable warning 4996.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} opengl32.lib)
ELSE()
  # Enable all pedantic warnings.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
  IF(APPLE)
    # Add required frameworks for GLFW.
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  ELSE()
    #Link the Linux OpenGL library
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "GL")
  ENDIF()
ENDIF()