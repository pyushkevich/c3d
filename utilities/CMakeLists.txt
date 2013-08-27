ADD_EXECUTABLE(c3d_affine_tool AffineTransformTool.cxx)
TARGET_LINK_LIBRARIES(c3d_affine_tool ${ITK_LIBRARIES})
INSTALL_TARGETS("/bin" c3d_affine_tool)
INSTALL(FILES bashcomp.sh DESTINATION share)
