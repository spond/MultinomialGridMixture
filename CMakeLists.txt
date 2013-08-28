cmake_minimum_required(VERSION 2.8)
project(MultiMix)

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)

#-------------------------------------------------------------------------------
# make it fast, scotty
#-------------------------------------------------------------------------------
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
    )

    set(GCC45 false)
    set(GCC46 false)

    if(${GCC_VERSION} VERSION_GREATER 4.5 OR ${GCC_VERSION} VERSION_EQUAL 4.5)
        set(GCC45 true)
    endif(${GCC_VERSION} VERSION_GREATER 4.5 OR ${GCC_VERSION} VERSION_EQUAL 4.5)

    if(${GCC_VERSION} VERSION_GREATER 4.6 OR ${GCC_VERSION} VERSION_EQUAL 4.6)
        set(GCC46 true)
    endif(${GCC_VERSION} VERSION_GREATER 4.6 OR ${GCC_VERSION} VERSION_EQUAL 4.6)

    if(${MACOSX_LION})
        set(DEFAULT_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-conversion-null -Wno-dangling-else -Wno-logical-op-parentheses")
    endif(${MACOSX_LION})

    if(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-g -fsigned-char -O3 ")
    else(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-g -fsigned-char -O3  ")
    endif(${GCC45})

    if(${GCC46})
        set(DEFAULT_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-conversion-null")
    endif(${GCC46})
    
    set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")
endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)


set(CMAKE_CONFIGURATION_TYPES Release)

include_directories(
    src/
)

add_executable(
    ngs_cmgs
    src/ComputeMultinomialGridScores.cpp
)

set_target_properties(
    ngs_cmgs
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}"
    LINK_FLAGS "${DEFAULT_LINK_FLAGS}"
)

install(
    TARGETS ngs_cmgs
	RUNTIME DESTINATION bin
	OPTIONAL
)

