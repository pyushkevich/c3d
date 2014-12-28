PROJECT(CONVERT3D)

# Nonsense for cmake 2.6 compatibility
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Versioning information
SET(C3D_VERSION_MAJOR 1)
SET(C3D_VERSION_MINOR 0)
SET(C3D_VERSION_PATCH 0)
SET(C3D_VERSION_FULL "${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}.${C3D_VERSION_PATCH}")

# Do we want to build the UI
OPTION(BUILD_GUI "Do you want to build the Qt-based C3D user interfaace?" OFF)

# Include directories
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/adapters)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/VoxBoIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/PovRayIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/N3ITK)

# Get ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  SOURCE_GROUP("Adapter Sources" REGULAR_EXPRESSION "adapters/*cxx")
  SOURCE_GROUP("Adapter Headers" REGULAR_EXPRESSION "adapters/*h")
ENDIF(WIN32)

SET(SOURCES
  adapters/AddImages.cxx
  adapters/AlignByLandmarks.cxx
  adapters/AntiAliasImage.cxx
  adapters/ApplyMetric.cxx
  adapters/BinaryHoleFill.cxx
  adapters/BiasFieldCorrection.cxx
  adapters/BiasFieldCorrectionN4.cxx
  adapters/BinaryImageCentroid.cxx
  adapters/BinaryMathOperation.cxx
  adapters/ClipImageIntensity.cxx
  adapters/ComputeFFT.cxx
  adapters/ComputeOverlaps.cxx
  adapters/ConnectedComponents.cxx
  adapters/CoordinateMap.cxx
  adapters/CopyTransform.cxx
  adapters/CreateImage.cxx
  adapters/CreateInterpolator.cxx
  adapters/ExtractRegion.cxx
  adapters/ExtractSlice.cxx
  adapters/FlipImage.cxx
  adapters/HessianObjectness.cxx
  adapters/HistogramMatch.cxx
  adapters/ImageERF.cxx
  adapters/ImageLaplacian.cxx
  adapters/GeneralLinearModel.cxx
  adapters/LabelOverlapMeasures.cxx
  adapters/LabelStatistics.cxx
  adapters/LandmarksToSpheres.cxx
  adapters/LevelSetSegmentation.cxx
  adapters/MathematicalMorphology.cxx
  adapters/MixtureModel.cxx
  adapters/MRFVote.cxx
  adapters/MultiplyImages.cxx
  adapters/NormalizeLocalWindow.cxx
  adapters/NormalizedCrossCorrelation.cxx
  adapters/OverlayLabelImage.cxx
  adapters/PadImage.cxx
  adapters/PeronaMalik.cxx
  adapters/PrintImageInfo.cxx
  adapters/Rank.cxx
  adapters/ReadImage.cxx
  adapters/ReciprocalImage.cxx
  adapters/ReorderStack.cxx
  adapters/ReplaceIntensities.cxx
  adapters/ResampleImage.cxx
  adapters/ResliceImage.cxx
  adapters/SampleImage.cxx
  adapters/ScaleShiftImage.cxx
  adapters/ScalarToRGB.cxx
  adapters/SetSform.cxx
  adapters/SetOrientation.cxx
  adapters/SignedDistanceTransform.cxx
  adapters/SLICSuperVoxel.cxx
  adapters/SmoothImage.cxx
  adapters/SplitMultilabelImage.cxx
  adapters/StapleAlgorithm.cxx
  adapters/TestImage.cxx
  adapters/ThresholdImage.cxx
  adapters/TileImages.cxx
  adapters/TrimImage.cxx
  adapters/UnaryMathOperation.cxx
  adapters/UpdateMetadataKey.cxx
  adapters/Vote.cxx
  adapters/VoxelwiseRegression.cxx
  adapters/WarpImage.cxx
  adapters/WarpLabelImage.cxx
  adapters/WeightedSum.cxx
  adapters/WeightedSumVoxelwise.cxx
  adapters/WrapDimension.cxx
  adapters/WriteImage.cxx
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx)

# Pull in FFTW
FIND_LIBRARY(FFTW_LIB fftw3f)

# Configure the version file
CONFIGURE_FILE(
  ${CONVERT3D_SOURCE_DIR}/ConvertImageVersion.cxx.in
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx @ONLY IMMEDIATE)

# Get the extra stuff compiled
SUBDIRS(itkextras)

ADD_LIBRARY(cnd_adapters ${SOURCES})
ADD_LIBRARY(cnd_driver ConvertImageND.cxx)

LINK_LIBRARIES(
  cnd_adapters cnd_driver ${ITK_LIBRARIES} ITKVoxBoIO ITKPovRayIO)

ADD_EXECUTABLE(c3d Convert3DMain.cxx)

ADD_EXECUTABLE(c2d Convert2DMain.cxx)

ADD_EXECUTABLE(c4d Convert4DMain.cxx)

ADD_EXECUTABLE(c3d_affine_tool utilities/AffineTransformTool.cxx)

SET(CMDLINE_TARGETS c2d c3d c4d c3d_affine_tool)

# =====================================================
# OPTIONAL GUI BUILD
# =====================================================
IF(BUILD_GUI)

  # CMake modules taken from the GUI dir
  SET(CMAKE_MODULE_PATH ${CONVERT3D_SOURCE_DIR}/gui/CMake)

  #--------------------------------------------------------------------------------
  # Find Qt5
  #--------------------------------------------------------------------------------
  FIND_PACKAGE(Qt5Widgets)
  SET(CONVERT3D_QT5_INCLUDE_DIRS ${Qt5Widgets_INCLUDE_DIRS})
  SET(CONVERT3D_QT5_LIBRARIES Qt5::Widgets)

  #--------------------------------------------------------------------------------
  # Specify source files and headers
  #--------------------------------------------------------------------------------
  INCLUDE_DIRECTORIES(
    ${CONVERT3D_SOURCE_DIR}/gui
    ${CONVERT3D_BINARY_DIR}
    ${CONVERT3D_QT5_INCLUDE_DIRS}
  )

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
  QT5_WRAP_UI(UI_FORM_HEADERS ${UI_FORMS})
  QT5_WRAP_CPP(UI_WRAPPED_MOC_HEADERS ${UI_MOC_HEADERS})

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

  # Define the main executable
  SET(C3DGUI_BUNDLE_NAME "Convert3DGUI")

  # Configure the executable's sources and libraries
  ADD_EXECUTABLE(${C3DGUI_BUNDLE_NAME} WIN32 MACOSX_BUNDLE ${C3DGUI_MAIN_SRC})
  TARGET_LINK_LIBRARIES(${C3DGUI_BUNDLE_NAME} ${CONVERT3D_QT5_LIBRARIES} cnd_adapters cnd_driver)

  #--------------------------------------------------------------------------------
  # Install the application
  #--------------------------------------------------------------------------------

  # On Apple, configure the application icon
  IF(APPLE)

    # set how it shows up in the Info.plist file
    SET(MACOSX_BUNDLE_ICON_FILE c3dgui.icns) 
    
    # set where in the bundle to put the icns file
    SET_SOURCE_FILES_PROPERTIES(${C3DGUI_OSX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

    # Set up the Info.plist file
    SET_TARGET_PROPERTIES(${C3DGUI_BUNDLE_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST
      ${CONVERT3D_SOURCE_DIR}/gui/resources/macos/Info.plist)

    # on Apple, the bundle is at the root of the
    # install tree, and on other platforms it'll go into the bin directory.
    INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} BUNDLE DESTINATION . COMPONENT Runtime)

    # Install the command line programs inside of the application
    INSTALL(TARGETS ${CMDLINE_TARGETS} DESTINATION ${C3DGUI_BUNDLE_NAME}.app/Contents/bin)

    # Include the qt4 dependent libraries
    include(DeployQt5)
    install_qt5_executable(${C3DGUI_BUNDLE_NAME}.app)

  ELSEIF(WIN32)

    # Set the PATH to include Qt libraries
	GET_FILENAME_COMPONENT(QT_BINARY_DIR "${Qt5Core_DIR}/../../../bin" ABSOLUTE)
	GET_FILENAME_COMPONENT(QT_LIBRARY_DIR "${Qt5Core_DIR}/../../" ABSOLUTE)

    # Install to the bin directory
    INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION bin)

    # Install the command line programs inside of the application
    INSTALL(TARGETS ${CMDLINE_TARGETS} DESTINATION bin)

    # Include the qt4 dependent libraries
    include(DeployQt5)

    # Make sure the GIF plugin is included
    get_property(QT_WIN_PLUGIN TARGET Qt5::QWindowsIntegrationPlugin PROPERTY LOCATION_RELEASE)

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
    SET(C3DGUI_MAIN_INSTALL_DIR lib/c3d_gui-${C3D_VERSION_FULL})
    SET(CPACK_PACKAGE_EXECUTABLES ${C3DGUI_EXE} ${C3DGUI_BUNDLE_NAME})
    
    INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION ${C3DGUI_MAIN_INSTALL_DIR})

    # Install the command line programs inside of the application
    INSTALL(TARGETS ${CMDLINE_TARGETS} DESTINATION bin)

    # On Linux, we generate forward shared executable
    SUBDIRS(gui/Utilities/Forwarding)

    include(DeployQt5)
    install_qt5_executable(${C3DGUI_MAIN_INSTALL_DIR}/${C3DGUI_BUNDLE_NAME} "qgif")

  ENDIF()

ELSE(BUILD_GUI)

  INSTALL(TARGETS c3d c2d c4d RUNTIME DESTINATION bin)

ENDIF(BUILD_GUI)

# Do the packaging
INCLUDE(cmake/Package.cmake)

# CDash Support
ENABLE_TESTING()
INCLUDE(CTest)
