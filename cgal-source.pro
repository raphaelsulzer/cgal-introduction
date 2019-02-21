TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.cpp \
#    readPlyWithCn.cpp \
    readPly.cpp \
#    info_insert_with_zip_iterator.cpp \
    exportTriWithCn.cpp \
#    exportTri.cpp \
    rayTriIntersection.cpp

QMAKE_CXXFLAGS += -std=c++11


INCLUDEPATH += usr/local/include

LIBS += \
-L/usr/local/lib/ -lCGAL \
-L/usr/lib/ -lgmp \
-L/usr/lib/ -lmpfr \


#LIBS +=     -L/usr/local/Cellar/boost/1.68.0_1/lib \ # BOOST
#            -L/usr/local/Cellar/mpfr/4.0.1/lib \ # MPFR
#            -L/usr/local/Cellar/gmp/6.1.2_2/lib \ # GMP
#            -L/usr/local/Cellar/cgal/4.13/lib # CGAL

#macx: INCLUDEPATH +=  /usr/local/Cellar/cgal/4.13/include \ # CGAL
#                /usr/local/Cellar/boost/1.68.0_1/include \ # BOOST
#                /usr/local/Cellar/gmp/6.1.2_2/include \ # GMP
#                /usr/local/Cellar/mpfr/4.0.1/include \ # MPFR


#macx: LIBS += -L/usr/local/lib/ -lgmp
#macx: LIBS += -L/usr/local/lib/ -lmpfr
#macx: LIBS += -L/usr/local/lib/ -lCGAL




HEADERS += \
