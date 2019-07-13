TEMPLATE = app                                                                                                                                     
TARGET = paris_heatmap
INCLUDEPATH += src
INCLUDEPATH += RNA_class

QT += core
QT += gui widgets

# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG -= app_bundle

# Input
HEADERS += pan_type.h paris.h paris_plot.h sam.h string_split.h param.h fasta.h shape.h sstructure.h align.h fold.h
SOURCES += paris_heatmap.cpp paris.cpp paris_plot.cpp sam.cpp string_split.cpp param.cpp pan_type.cpp fasta.cpp shape.cpp sstructure.cpp align.cpp fold.cpp

LIBS += -lhybrid -L"lib"

 
 