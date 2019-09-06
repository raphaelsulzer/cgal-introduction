#ifndef VECTORARITHMETIC_H
#define VECTORARITHMETIC_H

#include <cgal_typedefs.h>


double dot(const Vector& v1, const Vector& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Point& v1, const Point& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Vector& v1, const Point& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

double dot(const Point& v1, const Vector& v2){
    return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

Vector crossV(const Vector& v1, const Vector& v2){

    return Vector(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Vector crossV(const Point& v1, const Point& v2){

    return Vector(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point crossP(const Vector& v1, const Vector& v2){

    return Point(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point crossP(const Point& v1, const Point& v2){

    return Point(v1.y()*v2.z() - v1.z()*v2.y(),
                  v1.z()*v2.x() - v1.x()*v2.z(),
                  v1.x()*v2.y() - v1.y()*v2.x());
}

Point neg(const Point& p){

    return Point(-p.x(), -p.y(), -p.z());
}




#endif // VECTORARITHMETIC_H
