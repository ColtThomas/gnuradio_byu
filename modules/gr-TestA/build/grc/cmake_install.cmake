# Install script for directory: /usr/src/gnuradio/gnuradio/gr-TestA/grc

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gnuradio/grc/blocks" TYPE FILE FILES
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_my_qpsk_demod_cb.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_cleanslate.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_soqpsk_demod_cc.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_derp.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_soqpsk_df.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_matfile.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_qpsk_modulate_custom.xml"
    "/usr/src/gnuradio/gnuradio/gr-TestA/grc/TestA_lookup_table.xml"
    )
endif()

