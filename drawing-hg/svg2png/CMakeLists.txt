cmake_minimum_required(VERSION 2.6)

find_package(Qt4 4.1 COMPONENTS QtCore QtGui QtXml QtSvg REQUIRED )

SET( QT_USE_QTWEBKIT TRUE )

include_directories(
  ${QT_HEADERS_DIR}
  ${QT_QTCORE_INCLUDE_DIR}
  ${QT_QTGUI_INCLUDE_DIR}
  ${QT_QTSVG_INCLUDE_DIR}
)

set(svg2png_SRCS
  svg2png.cpp
)



add_executable(svg2png
  ${svg2png_SRCS}
)

target_link_libraries(svg2png
  ${QT_LIBRARIES}
  ${QT_QTCORE_LIBRARY}
  ${QT_QTGUI_LIBRARY}
  ${QT_QTSVG_LIBRARY}
)
