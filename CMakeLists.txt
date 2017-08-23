#----------------------------------------------------------------------------
# Setup the project
#
#FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
#
project(SDHCAL_G4Simu)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
option(USE_CXX11 "Use cxx11" True)

SET( Geant4_CONFIG_DEBUG ON)


option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
find_package(Geant4 REQUIRED qt ui_all vis_all)
else()
find_package(Geant4 REQUIRED)
endif()
#find_package( Marlin 1.0 REQUIRED ) # minimum required Marlin version
#find_package( ROOT REQUIRED ) # minimum required Marlin version
#find_package( CLHEP REQUIRED )
#find_package( LCIO REQUIRED )
#find_package( MarlinUtil REQUIRED )

#IF(NOT CRY_DIR)
#	MESSAGE( SEND_ERROR "Please define CRY_DIR !" )
#ENDIF()

SET( CRY_DIR "/home/hu/Stage/cry_v1.7")

SET(CRY_INCLUDE_DIRS ${CRY_DIR}/src)
SET(CRY_LIBRARIES ${CRY_DIR}/lib/libCRY.a)

INCLUDE_DIRECTORIES( ${CRY_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${CRY_LIBRARIES} )

# LCIO   MarlinUtil ROOT
#FOREACH( pkg  Geant4)
#    IF( ${pkg}_FOUND )
#        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
#        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
#        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )

#    ENDIF()
#ENDFOREACH()

#
# Add specific modules
#

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

MESSAGE ("G4 def : ${Geant4_DEFINITIONS}")
ADD_DEFINITIONS("-DG4VIS_USE_QT")
include(${Geant4_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

set(PROJECT_SRC
 
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonPrimaryGeneratorAction.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonPrimaryGeneratorMessenger.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/RNGWrapper.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonRunAction.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonAnalysisToolkit.cc        
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonSensitiveDetector.cc  
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonDetectorConstruction.cc   
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonTrackHit.cc 
	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonEventAction.cc
        ${CMAKE_CURRENT_SOURCE_DIR}/src/LyonEventActionMessenger.cc
#-----------------------------Obsolete------------------------------
#	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonGlobalMagneticField.cc
#	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonCalorimeter.cc
#       ${CMAKE_CURRENT_SOURCE_DIR}/src/CalParameterisation.cc  
#       ${CMAKE_CURRENT_SOURCE_DIR}/src/LyonDetectorFieldMessenger.cc
#	${CMAKE_CURRENT_SOURCE_DIR}/src/LyonDetectorFieldSetup.cc	
    )

set(PROJECT_HEADER


        ${CMAKE_CURRENT_SOURCE_DIR}/include/LyonConstants.hh
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonPrimaryGeneratorAction.hh  
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonPrimaryGeneratorMessenger.hh
	${CMAKE_CURRENT_SOURCE_DIR}/include/RNGWrapper.hh
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonRunAction.hh
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonAnalysisToolkit.hh       
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonSensitiveDetector.hh   
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonDetectorConstruction.hh 
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonTrackHit.hh  
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonEventAction.hh
        ${CMAKE_CURRENT_SOURCE_DIR}/include/LyonEventActionMessenger.hh
	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonRunManager.hh
#-----------------------------Obsolete------------------------------
#       ${CMAKE_CURRENT_SOURCE_DIR}/include/CalParameterisation.hh
#	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonDetectorFieldMessenger.hh
#	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonDetectorFieldSetup.hh
#	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonCalorimeter.hh
#	${CMAKE_CURRENT_SOURCE_DIR}/include/LyonGlobalMagneticField.hh 
    )

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries

add_executable(SDHCAL_Simu SimuLyonTest.cc ${PROJECT_HEADER} ${PROJECT_SRC})

#

target_link_libraries(SDHCAL_Simu ${Geant4_LIBRARIES} ${CRY_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

#install(TARGETS SDHCAL_Simu DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)
