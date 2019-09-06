#ifndef COLLISIONDETECTION_H
#define COLLISIONDETECTION_H

// Copyright (C) 2014 Anders Logg and August Johansson
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2014-02-03
// Last changed: 2014-04-03

#include <cgal_typedefs.h>
#include <vectorArithmetic.h>


bool collides_tetrahedron_tetrahedron(const Tetrahedron& tetrahedron_0,
                                             const Tetrahedron& tetrahedron_1);


// Helper function for collides_tetrahedron_tetrahedron: checks if
// plane pv1 is a separating plane. Stores local coordinates bc
// and the mask bit mask_edges.
bool separating_plane_face_A_1(const std::vector<Point>& pv1,
                  const Point& n,
                  std::vector<double>& bc,
                  int& mask_edges);

// Helper function for collides_tetrahedron_tetrahedron: checks if
// plane v1, v2 is a separating plane. Stores local coordinates bc
// and the mask bit mask_edges.
bool separating_plane_face_A_2(const std::vector<Point>& v1,
                  const std::vector<Point>& v2,
                  const Point& n,
                  std::vector<double>& bc,
                  int& mask_edges);

// Helper function for collides_tetrahedron_tetrahedron: checks if
// plane pv2 is a separating plane.
bool separating_plane_face_B_1(const std::vector<Point>& P_V2,
                  const Point& n)
{
  return ((dot(P_V2[0],n) > 0) &&
      (dot(P_V2[1],n) > 0) &&
      (dot(P_V2[2],n) > 0) &&
      (dot(P_V2[3],n) > 0));
}

// Helper function for collides_tetrahedron_tetrahedron: checks if
// plane v1, v2 is a separating plane.
bool separating_plane_face_B_2(const std::vector<Point>& V1,
                  const std::vector<Point>& V2,
                  const Point& n)
{
  return ((dot((V1[0] - V2[1]),n) > 0) &&
      (dot((V1[1] - V2[1]),n) > 0) &&
      (dot((V1[2] - V2[1]),n) > 0) &&
      (dot((V1[3] - V2[1]),n) > 0));
}

// Helper function for collides_tetrahedron_tetrahedron: checks if
// edge is in the plane separating faces f0 and f1.
static bool separating_plane_edge_A(const std::vector<std::vector<double> >& coord_1,
                const std::vector<int>& masks,
                int f0,
                int f1);



#endif // COLLISIONDETECTION_H
