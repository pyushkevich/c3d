#--------------------------------------------------------------------------------
# Find Qt4
#--------------------------------------------------------------------------------
FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

#--------------------------------------------------------------------------------
# Specify source files and headers
#--------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CONVERT3D_SOURCE_DIR}/gui
  ${CONVERT3D_BINARY_DIR}/gui
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
SET(C3DGUI_MAIN_SRC
    ${UI_QT_CXX}
    ${UI_WRAPPED_MOC_HEADERS}
    ${UI_MOC_HEADERS}
    ${UI_NONMOC_HEADERS}
    ${UI_FORM_HEADERS})

# On Apple, configure the application icon
IF(APPLE)

  # the icon file
  SET(C3DGUI_OSX_ICON ${CONVERT3D_SOURCE_DIR}/gui/resources/macos/c3dgui.icns)

  # set how it shows up in the Info.plist file
  SET(MACOSX_BUNDLE_ICON_FILE c3dgui.icns) 
  
  # set where in the bundle to put the icns file
  SET_SOURCE_FILES_PROPERTIES(${C3DGUI_OSX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

  # include the icns file in the target
  SET(C3DGUI_MAIN_SRC ${C3DGUI_MAIN_SRC} ${C3DGUI_OSX_ICON})

ENDIF(APPLE)

# Define the main executable
SET(C3DGUI_BUNDLE_NAME "Convert3DGUI")

# Configure the executable's sources and libraries
ADD_EXECUTABLE(${C3DGUI_BUNDLE_NAME} MACOSX_BUNDLE ${C3DGUI_MAIN_SRC})
TARGET_LINK_LIBRARIES(${C3DGUI_BUNDLE_NAME} ${QT_LIBRARIES} cnd_adapters cnd_driver)

# On apple, we also need our custom plist.info file to be configured
IF(APPLE)
  SET_TARGET_PROPERTIES(${C3DGUI_BUNDLE_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST
    ${CONVERT3D_SOURCE_DIR}/gui/resources/macos/Info.plist)
ENDIF(APPLE)

#--------------------------------------------------------------------------------
# Install the application
#--------------------------------------------------------------------------------

IF(APPLE)

  # on Apple, the bundle is at the root of the
  # install tree, and on other platforms it'll go into the bin directory.
  INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} 
      BUNDLE DESTINATION . COMPONENT Runtime)

  # Include the qt4 dependent libraries
  include(DeployQt4)
  install_qt4_executable(${C3DGUI_BUNDLE_NAME}.app "qsqlite;qgif")

ELSEIF(WIN32)

  # Install to the bin directory
  INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION bin)

  # Include the qt4 dependent libraries
  include(DeployQt4)
  install_qt4_executable(bin/${C3DGUI_BUNDLE_NAME}.exe "qsqlite;qgif")

  # CMake does not yet know to install into (x64) program files or not
  IF(CMAKE_CL_64)
    SET(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES64")
  ENDIF(CMAKE_CL_64)

  # Give it a windowsy directory name
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "C3D GUI ${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}")
  
  # On Win32, the executable is the actual exe
  SET(CPACK_PACKAGE_EXECUTABLES ${C3DGUI_BUNDLE_NAME})

ELSE()

  SET(C3DGUI_EXE "c3d_gui")
  SET(C3DGUI_MAIN_INSTALL_DIR lib/c3d_gui-${C3D_VERSION_FULL})
  SET(CPACK_PACKAGE_EXECUTABLES ${C3DGUI_EXE} ${C3DGUI_BUNDLE_NAME})
  
  INSTALL(TARGETS ${C3DGUI_BUNDLE_NAME} RUNTIME DESTINATION ${C3DGUI_MAIN_INSTALL_DIR})

  # On Linux, we generate forward shared executable
  SUBDIRS(Utilities/Forwarding)

  include(DeployQt4)
  install_qt4_executable(${C3DGUI_MAIN_INSTALL_DIR}/${C3DGUI_BUNDLE_NAME} "qgif")

ENDIF()

