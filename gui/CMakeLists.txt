PROJECT(C3DGUI)

#--------------------------------------------------------------------------------
# CMAKE PRELIMINARIES                       
#--------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.8)

#############################################
# REQUIRE QT                                #
#############################################
FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

INCLUDE_DIRECTORIES(
    ${C3DGUI_SOURCE_DIR}
    ${C3DGUI_BINARY_DIR}
)

SET(UI_QT_CXX
    CommandEditor.cxx
    ConvertSyntaxHighlighter.cxx
    MainWindow.cxx
    main.cxx
)

SET(UI_MOC_HEADERS
    CommandEditor.h
    ConvertSyntaxHighlighter.h
    MainWindow.h
)

SET(UI_NONMOC_HEADERS
)

SET(UI_FORMS
    MainWindow.ui
)

# Wrap the QT input files
QT4_WRAP_UI(UI_FORM_HEADERS ${UI_FORMS})
QT4_WRAP_CPP(UI_WRAPPED_MOC_HEADERS ${UI_MOC_HEADERS})

#--------------------------------------------------------------------------------
# Define main GUI executable
#--------------------------------------------------------------------------------
SET(MY_MAIN_SRC
    ${UI_QT_CXX}
    ${UI_WRAPPED_MOC_HEADERS}
    ${UI_MOC_HEADERS}
    ${UI_NONMOC_HEADERS}
    ${UI_FORM_HEADERS})

# On Apple, configure the application icon
IF(APPLE)

  # the icon file
  SET(MY_OSX_ICON ${C3DGUI_SOURCE_DIR}/resources/macos/c3dgui.icns)

  # set how it shows up in the Info.plist file
  SET(MACOSX_BUNDLE_ICON_FILE c3dgui.icns) 
  
  # set where in the bundle to put the icns file
  SET_SOURCE_FILES_PROPERTIES(${MY_OSX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

  # include the icns file in the target
  SET(MY_MAIN_SRC ${MY_MAIN_SRC} ${MY_OSX_ICON})

ENDIF(APPLE)

# Define the main executable
SET(MY_BUNDLE_NAME "c3d_gui")

# Configure the executable's sources and libraries
ADD_EXECUTABLE(${MY_BUNDLE_NAME} MACOSX_BUNDLE ${MY_MAIN_SRC})
TARGET_LINK_LIBRARIES(${MY_BUNDLE_NAME} ${QT_LIBRARIES} cnd_adapters cnd_driver)

# On apple, we also need our custom plist.info file to be configured
IF(APPLE)
  SET_TARGET_PROPERTIES(${MY_BUNDLE_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST
    ${C3DGUI_SOURCE_DIR}/resources/macos/Info.plist)
ENDIF(APPLE)

#--------------------------------------------------------------------------------
# Install the application

IF(APPLE)

  # on Apple, the bundle is at the root of the
  # install tree, and on other platforms it'll go into the bin directory.
  INSTALL(TARGETS ${MY_BUNDLE_NAME} 
      BUNDLE DESTINATION . COMPONENT Runtime)

  # Include the qt4 dependent libraries
  include(DeployQt4)
  install_qt4_executable(${MY_BUNDLE_NAME}.app "qsqlite;qgif")

ELSEIF(WIN32)

  # Install to the bin directory
  INSTALL(TARGETS ${MY_BUNDLE_NAME} RUNTIME DESTINATION bin)

  # Include the qt4 dependent libraries
  include(DeployQt4)
  install_qt4_executable(bin/${MY_BUNDLE_NAME}.exe "qsqlite;qgif")

  # CMake does not yet know to install into (x64) program files or not
  IF(CMAKE_CL_64)
    SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
  ENDIF(CMAKE_CL_64)

  # Give it a windowsy directory name
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "C3D GUI ${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}")
  
  # On Win32, the executable is the actual exe
  SET(CPACK_PACKAGE_EXECUTABLES c3d_gui "C3D_GUI")

ELSE()

  # On Linux, we generate forward shared executable
  #SUBDIRS(Utilities/Forwarding)

  #SET(SNAP_EXE "ITK-SNAP")
  #SET(SNAP_MAIN_INSTALL_DIR lib/snap-${SNAP_VERSION_FULL})
  #SET(SNAP_DATA_INSTALL_DIR ${SNAP_MAIN_INSTALL_DIR})
  #SET(CPACK_PACKAGE_EXECUTABLES "itksnap" "ITK-SNAP")
  #INSTALL(TARGETS ${SNAP_BUNDLE_NAME} RUNTIME DESTINATION ${SNAP_MAIN_INSTALL_DIR})

  include(DeployQt4)
  install_qt4_executable(${SNAP_MAIN_INSTALL_DIR}/${SNAP_EXE} "qgif")

ENDIF()

