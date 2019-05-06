TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
#    main.cpp \
#    readPlyWithCn.cpp \
#    readPly.cpp \
#    info_insert_with_zip_iterator.cpp \
#    exportTri.cpp \
#    rayTriIntersection.cpp \
#    main.cpp
    cgal-introduction.cpp

QMAKE_CXXFLAGS += -std=c++11

# surpress wired warnings
QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-variable -Wno-unused-parameter

### UBUNTU

#INCLUDEPATH += usr/local/include

#LIBS += \
#-L/usr/local/lib/ -lCGAL \
#-L/usr/lib/ -lgmp \
#-L/usr/lib/ -lmpfr \


### MAC OS
macx: INCLUDEPATH +=  /usr/local/Cellar/cgal/4.13/include \ # CGAL
                /usr/local/Cellar/gmp/6.1.2_2/include \ # GMP
                /usr/local/Cellar/mpfr/4.0.1/include \ # MPFR
                /usr/local/Cellar/boost/1.68.0_1/include # BOOST

macx: LIBS += -L/usr/local/lib/ -lgmp
macx: LIBS += -L/usr/local/lib/ -lmpfr
macx: LIBS += -L/usr/local/lib/ -lCGAL




HEADERS += \
