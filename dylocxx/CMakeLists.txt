project(project_libdylocxx CXX)

set(DART_IMPLEMENTATIONS_LIST ${DART_IMPLEMENTATIONS_LIST}
    PARENT_SCOPE)

set(ENABLE_HWLOC      ${ENABLE_HWLOC}
    PARENT_SCOPE)
set(ENABLE_LIBNUMA    ${ENABLE_LIBNUMA}
    PARENT_SCOPE)
set(ENABLE_PAPI       ${ENABLE_PAPI}
    PARENT_SCOPE)
set(ENABLE_LOGGING    ${ENABLE_LOGGING}
    PARENT_SCOPE)
set(ENABLE_ASSERTIONS ${ENABLE_ASSERTIONS}
    PARENT_SCOPE)
set(BUILD_TESTS       ${BUILD_TESTS}
    PARENT_SCOPE)

# ---------------------------------------------------------------------------
# Dependencies

include(${CMAKE_SOURCE_DIR}/CMakeExt/InstallHelpers.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeExt/DART.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeExt/hwloc.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeExt/PAPI.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeExt/NUMA.cmake)

set(Boost_USE_STATIC_LIBS     ON)
set(Boost_USE_MULTITHREADED  OFF)
set(Boost_USE_STATIC_RUNTIME  ON)

find_package(Boost 1.53.0 REQUIRED)

message(STATUS "Boost include directories: "
        ${Boost_INCLUDE_DIRS})
message(STATUS "Boost libraries: "
        ${Boost_LIBRARIES})


# ---------------------------------------------------------------------------
# Source Files

file(GLOB_RECURSE DYLOCXX_SOURCES
     "src/*.c" "src/*.h" "src/*.cc" "src/*.cpp" )
file(GLOB_RECURSE DYLOCXX_HEADERS
     "include/*.h")


# ---------------------------------------------------------------------------
# Linker- and Compiler Flags

set (ADDITIONAL_COMPILE_FLAGS
     "${ADDITIONAL_COMPILE_FLAGS} -DDYLOCXX")

if (ENABLE_LOGGING)
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_LOGGING")
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_LOGGING")
endif()

if (ENABLE_ASSERTIONS)
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_ASSERTIONS")
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_ASSERTIONS")
endif()


set (ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} rt)

if (PAPI_FOUND AND ENABLE_PAPI)
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_PAPI")
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_PAPI")
  set (ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES}
       ${PAPI_INCLUDE_DIRS})
  set (ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES}
       ${PAPI_LIBRARIES})
endif()
if (HWLOC_FOUND AND ENABLE_HWLOC)
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_HWLOC")
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_HWLOC")
  set (ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES}
       ${HWLOC_INCLUDE_DIRS})
  set (ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES}
       ${HWLOC_LIBRARIES})
  if (ENABLE_HWLOC_PCI)
    set (ADDITIONAL_COMPILE_FLAGS
         "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_HWLOC_PCI")
    set (ADDITIONAL_COMPILE_FLAGS
         "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_HWLOC_PCI")
  endif()
endif()
if (NUMA_FOUND AND ENABLE_LIBNUMA)
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDYLOC_ENABLE_NUMA")
  set (ADDITIONAL_COMPILE_FLAGS
       "${ADDITIONAL_COMPILE_FLAGS} -DDART_ENABLE_NUMA")
  set (ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES}
       ${NUMA_INCLUDE_DIRS})
  set (ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES}
       ${NUMA_LIBRARIES})
endif()


message (STATUS "libdylocxx additional compile flags:")
set(ADDITIONAL_COMPILE_FLAGS_STR "")
foreach (ADDITIONAL_FLAG ${ADDITIONAL_COMPILE_FLAGS})
  message (STATUS "    " ${ADDITIONAL_FLAG})
  set(ADDITIONAL_COMPILE_FLAGS_STR
      "${ADDITIONAL_COMPILE_FLAGS_STR} ${ADDITIONAL_FLAG}")
endforeach()
message (STATUS "libdylocxx additional libraries:")
foreach (ADDITIONAL_LIB ${ADDITIONAL_LIBRARIES})
  message (STATUS "    " ${ADDITIONAL_LIB})
endforeach()

# ---------------------------------------------------------------------------

set(DYLOCXX_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(DYLOCXX_SOURCE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(DYLOC_COMMON_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/common/include)

include_directories(
  ${DYLOCXX_INCLUDE_DIRS}
  ${DYLOC_COMMON_INCLUDE_DIRS}
  ${DART_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${ADDITIONAL_INCLUDES}
)

foreach (dart_variant ${DART_IMPLEMENTATIONS_LIST})

  set(DYLOCXX_LIBRARY      "dylocxx-${dart_variant}")
  set(DYLOC_COMMON_LIBRARY "dyloc-common-${dart_variant}")

  add_library(
    ${DYLOCXX_LIBRARY}
    ${DYLOCXX_HEADERS}
    ${DYLOCXX_SOURCES}
  )

  target_link_libraries(
    ${DYLOCXX_LIBRARY}
    ${DYLOC_COMMON_LIBRARY}
    ${DART_LIBRARIES}
    ${Boost_LIBRARIES}
    ${ADDITIONAL_LIBRARIES}
  )

  set_target_properties(
    ${DYLOCXX_LIBRARY} PROPERTIES
    COMPILE_FLAGS ${ADDITIONAL_COMPILE_FLAGS}
  )

  set_target_properties(
    ${DYLOCXX_LIBRARY} PROPERTIES
    CXX_STANDARD 14
  )
  set_target_properties(
    ${DYLOCXX_LIBRARY} PROPERTIES
    CXX_STANDARD_REQUIRED 14
  )

  # Library Installation
  install(TARGETS ${DYLOCXX_LIBRARY}
          DESTINATION lib
          EXPORT "${DYLOCXX_LIBRARY}-targets")

  if(${CMAKE_VERSION} VERSION_GREATER 3.0.0 )
    include(CMakePackageConfigHelpers)
    target_include_directories("${DYLOCXX_LIBRARY}" INTERFACE
                               $<INSTALL_INTERFACE:include>)

    configure_package_config_file(
      "dylocxx-config.cmake.in"
      "${DYLOCXX_LIBRARY}-config.cmake"
      INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake"
      PATH_VARS CMAKE_INSTALL_PREFIX)

    # Exports
    install(EXPORT "${DYLOCXX_LIBRARY}-targets"
            DESTINATION share/cmake)

  	# install custom config
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${DYLOCXX_LIBRARY}-config.cmake"
            DESTINATION share/cmake)
  endif()

endforeach() # dart_variant

# ---------------------------------------------------------------------------
# Test Build Target

if (BUILD_TESTS)
  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
## Installation

# Headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/dylocxx
        DESTINATION include FILES_MATCHING PATTERN "*.h")

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/dylocxx.h
        DESTINATION include)

