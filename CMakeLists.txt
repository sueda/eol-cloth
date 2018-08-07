# Author: Shinjiro Sueda
# Edited: Nick Weidner
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Name of the project
PROJECT(EOL-Cloth)

#-----------------------------------------------
# USER SETTINGS
#-----------------------------------------------

# COMPILER OPTIONS

OPTION(ONLINE "Build with online rendering supported" ON)

# Include support for one or multiple solvers
OPTION(MOSEK "Build with Mosek support" OFF)
OPTION(GUROBI "Build with Gurobi support" OFF)
OPTION(IGL "Build with libIGL support" OFF)


# LIBRARY PATHS
# These paths are machine specific
# These paths are examples

 SET(DEF_DIR_EIGEN "D:\\c++\\eigen-3.3.4")
 
 IF(${ONLINE})
	 SET(DEF_DIR_GLM "D:\\c++\\glm-0.9.8.4")
	 SET(DEF_DIR_GLFW "D:\\c++\\glfw-3.2.1")
	 SET(DEF_DIR_GLEW "D:\\c++\\glew-2.1.0")
	 ADD_DEFINITIONS(-DEOLC_ONLINE)
 ENDIF()
 
 IF(${MOSEK}) 
	 SET(DEF_DIR_MOSEK "C:\\Program Files\\Mosek\\8\\tools\\platform\\win64x86") # Location of /h and /bin
	 SET(DEF_BINARY_MOSEK "mosek64_8_1.lib") # Platform specific
	 ADD_DEFINITIONS(-DkMosekLicensePath="C:/Users/weid691/mosek/mosek.lic") # Path to license file
	 ADD_DEFINITIONS(-DEOLC_MOSEK)
 ENDIF()
 
 IF(${GUROBI})
	 SET(DEF_DIR_GUROBI "C:\\gurobi752\\win64")
	 SET(DEF_BINARY_GUROBI "gurobi75.lib") # Platform specific
	 SET(DEF_BINARY_GUROBI_RELEASE "gurobi_c++md2017.lib") # Platform specific
	 SET(DEF_BINARY_GUROBI_DEBUG "gurobi_c++mdd2017.lib") # Platform specific
	 ADD_DEFINITIONS(-DEOLC_GUROBI)
 ENDIF()
 
 IF(${IGL})
	 SET(DEF_DIR_IGL "C:\\c++\\libigl-1.2.1\\include") # Does not build static lib
	 ADD_DEFINITIONS(-DEOLC_IGL)
 ENDIF()
 
# Additional defines
 ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
 ADD_DEFINITIONS(-D_CRT_NO_VA_START_VALIDATION)

#-----------------------------------------------
# END USER SETTINGS, DON'T TOUCH BELOW
#-----------------------------------------------

# Use glob to get the list of all source files.
# We don't really need to include header and resource files to build, but it's
# nice to have them also show up in IDEs.
 FILE(GLOB SOURCES "src/*.cpp" "src/*.c")
 FILE(GLOB HEADERS "src/*.hpp" "src/*.h")
 FILE(GLOB_RECURSE GLSL "resources/*.glsl")
 FILE(GLOB ARCSIMC "src/external/ArcSim/*.cpp" "src/external/ArcSim/*.c")
 FILE(GLOB ARCSIMH "src/external/ArcSim/*.hpp" "src/external/ArcSim/*.h")
 FILE(GLOB JSONC "src/external/Json/*.cpp" "src/external/Json/*.c")
 FILE(GLOB JSONH "src/external/Json/*.hpp" "src/external/Json/*.h")
 IF(${ONLINE})
	FILE(GLOB ONLINEC "src/online/*.cpp" "src/online/*.c")
	FILE(GLOB ONLINEH "src/online/*.hpp" "src/online/*.h")
 ENDIF()
 IF(${MOSEK})
	FILE(GLOB MOSEKC "src/external/SolverWrappers/Mosek/*.cpp" "src/external/SolverWrappers/Mosek/*.c")
	FILE(GLOB MOSEKH "src/external/SolverWrappers/Mosek/*.hpp" "src/external/SolverWrappers/Mosek/*.h")
 ENDIF()
 IF(${GUROBI})
	FILE(GLOB GUROBIC "src/external/SolverWrappers/Gurobi/*.cpp" "src/external/SolverWrappers/Gurobi/*.c")
	FILE(GLOB GUROBIH "src/external/SolverWrappers/Gurobi/*.hpp" "src/external/SolverWrappers/Gurobi/*.h")
 ENDIF()
 # FILE(GLOB_RECURSE SOURCES 
	# RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
	# ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)
 # FILE(GLOB_RECURSE EXTERNAL 
	# RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
	# ${CMAKE_CURRENT_SOURCE_DIR}/external/ArcSim/*.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/external/ArcSim/*.hpp)
	
# SOURCE_GROUP(Source FILES ${SOURCES})
# SOURCE_GROUP(External FILES ${EXTERNAL})

# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL} ${ARCSIMC} ${ARCSIMH} ${JSONC} ${JSONH} ${ONLINEC} ${ONLINEH} ${MOSEKC} ${MOSEKH} ${GUROBIC} ${GUROBIH})

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

IF(${ONLINE})
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
ENDIF()

# Get the MOSEK environment variable.
IF(${MOSEK})
	SET(MOSEK_DIR "$ENV{MOSEK_DIR}")
	IF(NOT MOSEK_DIR)
		# The environment variable was not set
		SET(ERR_MSG "Please point the environment variable MOSEK_DIR to the root directory of your MOSEK installation.")
		MESSAGE(STATUS "Looking for MOSEK in ${DEF_DIR_MOSEK}")
		IF(IS_DIRECTORY ${DEF_DIR_MOSEK})
			MESSAGE(STATUS "Found!")
			SET(MOSEK_DIR ${DEF_DIR_MOSEK})
		ELSE()
			MESSAGE(FATAL_ERROR ${ERR_MSG})
		ENDIF()
	ENDIF()
	INCLUDE_DIRECTORIES(${MOSEK_DIR}/h)
	# With prebuilt binaries
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MOSEK_DIR}/bin/${DEF_BINARY_MOSEK})
ENDIF()

# Get the MOSEK environment variable.
IF(${GUROBI})
	SET(GUROBI_DIR "$ENV{GUROBI_DIR}")
	IF(NOT GUROBI_DIR)
		# The environment variable was not set
		SET(ERR_MSG "Please point the environment variable GUROBI_DIR to the root directory of your GUROBI installation.")
		MESSAGE(STATUS "Looking for GUROBI in ${DEF_DIR_GUROBI}")
		IF(IS_DIRECTORY ${DEF_DIR_GUROBI})
			MESSAGE(STATUS "Found!")
			SET(GUROBI_DIR ${DEF_DIR_GUROBI})
		ELSE()
			MESSAGE(FATAL_ERROR ${ERR_MSG})
		ENDIF()
	ENDIF()
	INCLUDE_DIRECTORIES(${GUROBI_DIR}/include)
	# With prebuilt binaries
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GUROBI_DIR}/lib/${DEF_BINARY_GUROBI})
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} debug ${GUROBI_DIR}/lib/${DEF_BINARY_GUROBI_DEBUG})
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} optimized ${GUROBI_DIR}/lib/${DEF_BINARY_GUROBI_RELEASE})
ENDIF()

# Get the IGL environment variable.
IF(${IGL})
	SET(IGL_INCLUDE_DIR "$ENV{IGL_INCLUDE_DIR}")
	IF(NOT IGL_INCLUDE_DIR)
		# The environment variable was not set
		SET(ERR_MSG "Please point the environment variable IGL_INCLUDE_DIR to the root directory of your IGL installation.")
	  		MESSAGE(STATUS "Looking for IGL in ${DEF_DIR_IGL}")
		IF(IS_DIRECTORY ${DEF_DIR_IGL})
			MESSAGE(STATUS "Found!")
			SET(GUROBI_DIR ${DEF_DIR_IGL})
		ELSE()
			MESSAGE(FATAL_ERROR ${ERR_MSG})
		ENDIF()
	ENDIF()
	INCLUDE_DIRECTORIES(${IGL_INCLUDE_DIR})
ENDIF()

# OS specific options and libraries
IF(WIN32)
  # c++11 is enabled by default.
  # -Wall produces way too many warnings.
  # -pedantic is not supported.
  # Disable warning 4996.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
  IF(${ONLINE})
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} opengl32.lib)
  ELSE()
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME})
  ENDIF()
ELSE()
  # Enable all pedantic warnings.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
  IF(${ONLINE})
	  IF(APPLE)
		# Add required frameworks for GLFW.
		TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
	  ELSE()
		#Link the Linux OpenGL library
		TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "GL")
	  ENDIF()
	  ENDIF()
ENDIF()