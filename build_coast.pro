INCLUDEPATH += /c/PFM/compile/include
LIBS += -L /c/PFM/compile/lib -lnvutility -lgdal -lxml2 -lpoppler -lz -lm -liconv
DEFINES += NVWIN3X
CONFIG += console
CONFIG -= qt
QMAKE_LFLAGS += 
######################################################################
# Automatically generated by qmake (2.01a) Fri Dec 11 16:04:08 2020
######################################################################

TEMPLATE = app
TARGET = build_coast
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += version.h
SOURCES += main.c