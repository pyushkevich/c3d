PROJECT(CONVERT3D)

# CMake compatibility code
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
if(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0077 NEW)
  CMAKE_POLICY(SET CMP0042 NEW)
endif(COMMAND cmake_policy)

#--------------------------------------------------------------------------------
# VERSION
#--------------------------------------------------------------------------------
INCLUDE(CMake/PYVersion.cmake)

# Set the semantic version and version date
VERSION_VARS(1 4 0 "" "20220519" "May 19, 2022")

# Enable languages
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

# Specify the C++ standard
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED True)

# Option allowing to build c3d as a subproject of another project
OPTION(CONVERT3D_BUILD_AS_SUBPROJECT "Build Convert3D as a sub-project of a larger project" OFF)
MARK_AS_ADVANCED(CONVERT3D_BUILD_AS_SUBPROJECT)

# Options corresponding to remote ITK modules
OPTION(CONVERT3D_USE_ITK_REMOTE_MODULES "Use remote modules from ITK" ON)

# If we are building as a sub-project we skip all of this extra stuff
IF(NOT CONVERT3D_BUILD_AS_SUBPROJECT)

  # Do we want to build the UI
  OPTION(BUILD_GUI "Do you want to build the Qt-based C3D user interface?" OFF)

  # List of ITK modules
  SET(ITK_MODULE_LIST
    ITKAnisotropicSmoothing
    ITKAntiAlias
    ITKBiasCorrection
    ITKBinaryMathematicalMorphology
    ITKColormap
    ITKCommon
    ITKConnectedComponents
    ITKConvolution
    ITKDistanceMap
    ITKFFT
    ITKGDCM
    ITKIOGDCM
    ITKIOImageBase
    ITKIOTransformBase
    ITKImageAdaptors
    ITKImageCompare
    ITKImageCompose
    ITKImageFeature
    ITKImageFilterBase
    ITKImageFunction
    ITKImageGradient
    ITKImageGrid
    ITKImageIntensity
    ITKImageLabel
    ITKImageNoise
    ITKImageStatistics
    ITKLabelMap
    ITKLabelVoting
    ITKLevelSets
    ITKMathematicalMorphology
    ITKRegistrationCommon
    ITKSmoothing
    ITKStatistics
    ITKThresholding
    ITKTransform
    )

  # Add ITK version specific modules
  FIND_PACKAGE(ITK "4.7" CONFIG REQUIRED)
  MESSAGE(STATUS "ITK")
  MESSAGE(STATUS "  Version: ${ITK_VERSION}")
  IF(ITK_VERSION VERSION_LESS "4.13")
    FOREACH(mod IN LISTS ITK_MODULES_ENABLED)
      IF(${mod} MATCHES "IO")
        LIST(APPEND ITK_MODULE_LIST ${mod})
      ENDIF()
    ENDFOREACH()
  ELSE(ITK_VERSION VERSION_LESS "4.13")
  SET(ITK_MODULE_LIST ${ITK_MODULE_LIST}
    ITKTransformFactory
    ITKImageIO
    ITKTransformIO)
  ENDIF(ITK_VERSION VERSION_LESS "4.13")

  IF(ITK_VERSION VERSION_LESS "4.11")
    IF(CONVERT3D_USE_ITK_REMOTE_MODULES)
      SET(CONVERT3D_USE_ITK_REMOTE_MODULES OFF)
      MESSAGE(WARNING "CONVERT3D_USE_ITK_REMOTE_MODULES has been automatically disabled because ITK version is less than 4.11 which introduced the MorphologicalContourInterpolation that it requires, see:
      https://github.com/KitwareMedical/ITKMorphologicalContourInterpolation")
    ENDIF(CONVERT3D_USE_ITK_REMOTE_MODULES)
  ENDIF(ITK_VERSION VERSION_LESS "4.11")

  # Add the remote modules
  IF(CONVERT3D_USE_ITK_REMOTE_MODULES)
    ADD_DEFINITIONS(-DCONVERT3D_USE_ITK_REMOTE_MODULES)
    SET(ITK_MODULE_LIST ${ITK_MODULE_LIST} 
      MorphologicalContourInterpolation)
  ENDIF(CONVERT3D_USE_ITK_REMOTE_MODULES)

  # Get ITK
  FIND_PACKAGE(ITK REQUIRED COMPONENTS ${ITK_MODULE_LIST})
  INCLUDE(${ITK_USE_FILE})

  # Command-line tools are being built
  SET(BUILD_CLI ON)
  SET(INSTALL_CLI ON)

ELSE()

  # Should the CLI tools be compiled or only the libraries?
  OPTION(CONVERT3D_SUBPROJECT_BUILD_CLI_TOOLS "Should the Convert3D command-line tools be built?" OFF)

  # Should the CLI tools be installed? Then also provide CONVERT3D_SUBPROJECT_CLI_INSTALL_PATH
  OPTION(CONVERT3D_SUBPROJECT_INSTALL_CLI_TOOLS "Should the Convert3D command-line tools be installed?" OFF)
  MARK_AS_ADVANCED(CONVERT3D_SUBPROJECT_BUILD_CLI_TOOLS CONVERT3D_SUBPROJECT_INSTALL_CLI_TOOLS)

  SET(BUILD_CLI ${CONVERT3D_SUBPROJECT_BUILD_CLI_TOOLS})
  SET(INSTALL_CLI ${CONVERT3D_SUBPROJECT_INSTALL_CLI_TOOLS})

ENDIF()

# Include the library file
INCLUDE(${CONVERT3D_SOURCE_DIR}/ConvertNDLibrary.cmake)

# Build the command-line executables
IF(BUILD_CLI)

  # Command-line executables
  ADD_EXECUTABLE(c3d Convert3DMain.cxx)
  ADD_EXECUTABLE(c2d Convert2DMain.cxx)
  ADD_EXECUTABLE(c4d Convert4DMain.cxx)
  ADD_EXECUTABLE(c3d_affine_tool utilities/AffineTransformTool.cxx)

  FOREACH(target c3d c2d c4d c3d_affine_tool)
    TARGET_INCLUDE_DIRECTORIES(${target} PRIVATE ${CONVERT3D_INCLUDE_DIRS})
    TARGET_INCLUDE_DIRECTORIES(${target} 
      INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>)

    TARGET_LINK_LIBRARIES(${target} PUBLIC ${C3D_LINK_LIBRARIES})
  ENDFOREACH()

ENDIF(BUILD_CLI)

# Define the main executable
SET(C3DGUI_BUNDLE_NAME "Convert3DGUI")

# Install command-line executables
IF(INSTALL_CLI)
  IF(CONVERT3D_SUBPROJECT_CLI_INSTALL_PATH)
    SET(CLI_INSTALL_PATH ${CONVERT3D_SUBPROJECT_CLI_INSTALL_PATH})
    SET(LIB_INSTALL_PATH lib)
    SET(INC_INSTALL_PATH include)
  ELSEIF(APPLE AND BUILD_GUI)
    SET(CLI_INSTALL_PATH "${C3DGUI_BUNDLE_NAME}.app/Contents/bin")
    SET(LIB_INSTALL_PATH "${C3DGUI_BUNDLE_NAME}.app/Contents/lib")
    SET(INC_INSTALL_PATH "${C3DGUI_BUNDLE_NAME}.app/Contents/include")
  ELSE()
    SET(CLI_INSTALL_PATH bin)
    SET(LIB_INSTALL_PATH lib)
    SET(INC_INSTALL_PATH include)
  ENDIF()

  # Install the executables
  INSTALL(TARGETS c2d c3d c4d c3d_affine_tool DESTINATION ${CLI_INSTALL_PATH} COMPONENT Runtime)

  # Install the libraries
  INSTALL(TARGETS cnd_driver cnd_api cnd_maxflow
          DESTINATION ${LIB_INSTALL_PATH}
          EXPORT Convert3DTargets)

  # Install the API header
  INSTALL(FILES api/ConvertAPI.h DESTINATION ${INC_INSTALL_PATH})

ENDIF(INSTALL_CLI)

IF(NOT CONVERT3D_BUILD_AS_SUBPROJECT)

  IF(WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
    SOURCE_GROUP("Adapter Sources" REGULAR_EXPRESSION "adapters/*cxx")
    SOURCE_GROUP("Adapter Headers" REGULAR_EXPRESSION "adapters/*h")
    
    # On windows we want to expand wildcards
    SET_TARGET_PROPERTIES(c3d PROPERTIES LINK_FLAGS "/link setargv.obj")
  ENDIF(WIN32)

  # =====================================================
  # OPTIONAL GUI BUILD
  # =====================================================
  IF(BUILD_GUI)

    # Turn off warning
    CMAKE_POLICY(SET CMP0080 OLD)

    # CMake modules taken from the GUI dir
    SET(CMAKE_MODULE_PATH ${CONVERT3D_SOURCE_DIR}/gui/CMake)

    #--------------------------------------------------------------------------------
    # Find Qt6
    #--------------------------------------------------------------------------------
    FIND_PACKAGE(Qt6Widgets)
    SET(CONVERT3D_QT_INCLUDE_DIRS ${Qt6Widgets_INCLUDE_DIRS})
    SET(CONVERT3D_QT_LIBRARIES Qt6::Widgets)
    GET_FILENAME_COMPONENT(QT_BINARY_DIR "${Qt6Core_DIR}/../../../bin" ABSOLUTE)
    GET_FILENAME_COMPONENT(QT_LIBRARY_DIR "${Qt6Core_DIR}/../../" ABSOLUTE)

    #--------------------------------------------------------------------------------
    # Specify source files and headers
    #--------------------------------------------------------------------------------
    SET(UI_QT_CXX 
      gui/CommandEditor.cxx 
      gui/ConvertSyntaxHighlighter.cxx 
      gui/HistoryDialog.cxx
      gui/MainWindow.cxx
      gui/SettingsDialog.cxx
      gui/main.cxx)

    SET(UI_MOC_HEADERS
      gui/CommandEditor.h
      gui/ConvertSyntaxHighlighter.h
      gui/HistoryDialog.h
      gui/MainWindow.h
      gui/SettingsDialog.h)

    SET(UI_FORMS 
      gui/HistoryDialog.ui
      gui/MainWindow.ui
      gui/SettingsDialog.ui)

    SET(UI_NONMOC_HEADERS)

    # Wrap the QT input files
    QT_WRAP_UI(UI_FORM_HEADERS ${UI_FORMS})
    QT_WRAP_CPP(UI_WRAPPED_MOC_HEADERS ${UI_MOC_HEADERS})

    # Configure the OS-specific parts of the GUI exe
    IF(APPLE)
      SET(C3DGUI_OSX_ICON ${CONVERT3D_SOURCE_DIR}/gui/resources/macos/c3dgui.icns)
      SET(UI_OS_EXTRAS ${C3DGUI_OSX_ICON})
    ENDIF(APPLE)

    #--------------------------------------------------------------------------------
    # Define main GUI executable
    #--------------------------------------------------------------------------------
    SET(C3DGUI_MAIN_SRC ${UI_QT_CXX} ${UI_WRAPPED_MOC_HEADERS} ${UI_MOC_HEADERS} 
      ${UI_NONMOC_HEADERS} ${UI_FORM_HEADERS} ${UI_OS_EXTRAS})

    # Configure the executable's sources and libraries
    ADD_EXECUTABLE(${C3DGUI_BUNDLE_NAME} WIN32 MACOSX_BUNDLE ${C3DGUI_MAIN_SRC})
    TARGET_LINK_LIBRARIES(${C3DGUI_BUNDLE_NAME} ${CONVERT3D_QT_LIBRARIES} ${C3D_LINK_LIBRARIES})

    # Configure the include path for the GUI
    TARGET_INCLUDE_DIRECTORIES(
      ${C3DGUI_BUNDLE_NAME} PUBLIC
      ${CONVERT3D_INCLUDE_DIRS}
      ${CONVERT3D_SOURCE_DIR}/gui
      ${CONVERT3D_BINARY_DIR}
      ${CONVERT3D_QT_INCLUDE_DIRS}
    )

    #--------------------------------------------------------------------------------
    # Install the application
    #--------------------------------------------------------------------------------

    # On Apple, configure the application icon
    IF(APPLE)

      # set how it shows up in the Info.plist file
      SET(CPACK_BUNDLE_ICON ${C3DGUI_OSX_ICON})
      SET(CPACK_PACKAGE_ICON ${C3DGUI_OSX_ICON})

      # Code signing certificate - optional
      SET(CONVERT3D_MACOSX_CODESIGN_CERT "" CACHE STRING "Name of the Apple Developer Certificate to sign application")
      MARK_AS_ADVANCED(CONVERT3D_MACOSX_CODESIGN_CERT)
      IF(CONVERT3D_MACOSX_CODESIGN_CERT)
        SET(CPACK_BUNDLE_APPLE_CERT_APP ${CONVERT3D_MACOSX_CODESIGN_CERT})
      ENDIF()

      # set where in the bundle to put the icns file
      SET_SOURCE_FILES_PROPERTIES(${C3DGUI_OSX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

      # Set up the Info.plist file
      SET_TARGET_PROPERTIES(${C3DGUI_BUNDLE_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST
        ${CONVERT3D_SOURCE_DIR}/gui/resources/macos/Info.plist)

      # on Apple, the bundle is at the root of the
      # install tree, and on other platforms it'll go into the bin directory.
      INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} BUNDLE DESTINATION . COMPONENT Runtime)

      # Include the qt4 dependent libraries
      include(DeployQt5)
      get_filename_component(QT_COCOA_PLUGIN "${QT_BINARY_DIR}/../plugins/platforms/libqcocoa.dylib" ABSOLUTE)
      install_qt5_executable(${C3DGUI_BUNDLE_NAME}.app "${QT_COCOA_PLUGIN}")

    ELSEIF(WIN32)

      # Install to the bin directory
      INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION bin)

      # Include the qt4 dependent libraries
      include(DeployQt5)

      # Make sure the GIF plugin is included
      get_filename_component(QT_WIN_PLUGIN "${QT_BINARY_DIR}/../plugins/platforms/qwindows.dll" ABSOLUTE)

      # Install with the plugin
      install_qt5_executable(bin/${C3DGUI_BUNDLE_NAME}.exe "${QT_WIN_PLUGIN}")

      # On windows, we have to configure NSIS
      SET(CPACK_NSIS_INSTALLED_ICON_NAME "${C3DGUI_BUNDLE_NAME}.exe")
      SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Convert3D")
      SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.itksnap.org/c3d")
      SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.itksnap.org/c3d")

      # CMake does not yet know to install into (x64) program files or not
      IF(CMAKE_CL_64)
        SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
      ENDIF(CMAKE_CL_64)

      # Give it a windowsy directory name
      SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Convert3D")
      
      # On Win32, the executable is the actual exe
      SET(CPACK_PACKAGE_EXECUTABLES ${C3DGUI_BUNDLE_NAME} "Convert3D GUI")

    ELSE()

      SET(C3DGUI_EXE "c3d_gui")
      SET(C3DGUI_MAIN_INSTALL_DIR lib/c3d_gui-${CONVERT3D_VERSION_FULL})
      SET(CPACK_PACKAGE_EXECUTABLES ${C3DGUI_EXE} ${C3DGUI_BUNDLE_NAME})
      
      INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION ${C3DGUI_MAIN_INSTALL_DIR})

      # On Linux, we generate forward shared executable
      SUBDIRS(gui/Utilities/Forwarding)

      include(DeployQt5)
      install_qt5_executable(${C3DGUI_MAIN_INSTALL_DIR}/${C3DGUI_BUNDLE_NAME} "qgif")

    ENDIF()

  ENDIF(BUILD_GUI)

  # Do the packaging
  INCLUDE(CMake/Package.cmake)

  # CDash Support
  ENABLE_TESTING()
  INCLUDE(CTest)

  # This must be at the end
  IF(INSTALL_CLI)
    install(EXPORT Convert3DTargets
            FILE Convert3DTargets.cmake
            DESTINATION ${LIB_INSTALL_PATH}/cmake/Convert3D)

    include(CMakePackageConfigHelpers)

    set ( Convert3D_INCLUDE_DIRS "include" )

    # generate the config file that is includes the exports
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
      "${CMAKE_CURRENT_BINARY_DIR}/Convert3DConfig.cmake"
      INSTALL_DESTINATION "${LIB_INSTALL_PATH}/cmake/example"
      NO_SET_AND_CHECK_MACRO
      NO_CHECK_REQUIRED_COMPONENTS_MACRO
      PATH_VARS Convert3D_INCLUDE_DIRS)

    # generate the version file for the config file
    write_basic_package_version_file(
      "${CMAKE_CURRENT_BINARY_DIR}/Convert3DConfigVersion.cmake"
      VERSION "${CONVERT3D_VERSION_FULL}"
      COMPATIBILITY AnyNewerVersion
    )

    # install the configuration file
    install(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/Convert3DConfig.cmake
      DESTINATION ${LIB_INSTALL_PATH}/cmake/Convert3D
      )

  ENDIF(INSTALL_CLI)



ENDIF(NOT CONVERT3D_BUILD_AS_SUBPROJECT)
