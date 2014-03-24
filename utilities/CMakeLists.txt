ADD_EXECUTABLE(c3d_affine_tool AffineTransformTool.cxx)
TARGET_LINK_LIBRARIES(c3d_affine_tool ${ITK_LIBRARIES})

INSTALL(TARGETS c3d_affine_tool RUNTIME DESTINATION bin)
INSTALL(FILES bashcomp.sh DESTINATION share)
