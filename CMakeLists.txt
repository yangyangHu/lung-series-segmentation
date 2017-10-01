#################################################################
# Edit by HuYangYang and JiGuohua
# 2014-3-18
#################################################################

cmake_minimum_required( VERSION 2.8 )
project( lung_seg )

find_package (ITK)
if (ITK_FOUND)
  include ( ${ITK_USE_FILE} )
else (ITK_FOUND)
  message (FATAL_ERROR "Cannot build without ITK, Please set ITK_DIR.")
endif (ITK_FOUND)

add_executable( lung_seg lung_seg.cxx )

target_link_libraries( lung_seg
  ${ITK_LIBRARIES}
 )
