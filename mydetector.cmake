SET( CMAKE_PREFIX_PATH "/home/hu/Stage/geant4" CACHE PATH "CMAKE_PREFIX_PATH" FORCE )

option(USE_CXX11 "Use cxx11" True)

SET( CRY_DIR "/home/hu/Stage/cry_v1.7" CACHE PATH "Path to CRY package" FORCE)
MARK_AS_ADVANCED( CRY_DIR )
