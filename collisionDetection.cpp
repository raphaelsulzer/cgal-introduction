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
// Modified by Chris Richardson, 2014.
//
// First added:  2014-02-03
// Last changed: 2014-04-03
//
//-----------------------------------------------------------------------------
// Special note regarding the function collides_tetrahedron_tetrahedron
//-----------------------------------------------------------------------------
//
// The source code for the tetrahedron-tetrahedron collision test is
// from Fabio Ganovelli, Federico Ponchio and Claudio Rocchini: Fast
// Tetrahedron-Tetrahedron Overlap Algorithm, Journal of Graphics
// Tools, 7(2), 2002, and is under the following copyright:
//
// Visual Computing Group
// IEI Institute, CNUCE Institute, CNR Pisa
//
// Copyright(C) 2002 by Fabio Ganovelli, Federico Ponchio and Claudio
// Rocchini
//
// All rights reserved.
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without
// fee, provided that the above copyright notice appear in all copies
// and that both that copyright notice and this permission notice
// appear in supporting documentation. the author makes no
// representations about the suitability of this software for any
// purpose. It is provided "as is" without express or implied
// warranty.
//
//-----------------------------------------------------------------------------
// Special note regarding the function collides_triangle_triangle
//-----------------------------------------------------------------------------
//
// The source code for the triangle-triangle collision test is from
// Tomas Moller: A Fast Triangle-Triangle Intersection Test, Journal
// of Graphics Tools, 2(2), 1997, and is in the public domain.
//
//-----------------------------------------------------------------------------

#include <collisionDetection.h>
#include <vectorArithmetic.h>


//-----------------------------------------------------------------------------
bool
collides_tetrahedron_tetrahedron
(const Tetrahedron& tetrahedron_0,
 const Tetrahedron& tetrahedron_1)
{
  // This algorithm checks whether two tetrahedra intersect.

  // Algorithm and source code from Fabio Ganovelli, Federico Ponchio
  // and Claudio Rocchini: Fast Tetrahedron-Tetrahedron Overlap
  // Algorithm, Journal of Graphics Tools, 7(2), 2002. DOI:
  // 10.1080/10867651.2002.10487557. Source code available at
  // http://web.archive.org/web/20031130075955/http://www.acm.org/jgt/papers/GanovelliPonchioRocchini02/tet_a_tet.html


  // Get the vertices as points
  std::vector<Point> V1(4), V2(4);
  for (std::size_t i = 0; i < 4; ++i)
  {
    V1[i] = tetrahedron_0.vertex(i);
    V2[i] = tetrahedron_1.vertex(i);;
  }

  // Get the vectors between V2 and V1[0]
  std::vector<Point> P_V1(4);
  for (std::size_t i = 0; i < 4; ++i){
      Vector temp = V2[i]-V1[0];
      P_V1[i] = Point(temp.x(), temp.y(), temp.z());

  }

  // Data structure for edges of V1 and V2
  std::vector<Vector> e_v1(5), e_v2(5);
  e_v1[0] = V1[1] - V1[0];
  e_v1[1] = V1[2] - V1[0];
  e_v1[2] = V1[3] - V1[0];
  Point n = crossP(e_v1[1], e_v1[0]);

  // Maybe flip normal. Normal should be outward.
  if (dot(n, e_v1[2]) > 0)
    n = neg(n);
  std::vector<int> masks(4);
  std::vector<std::vector<double>> Coord_1(4, std::vector<double>(4));
  if (separating_plane_face_A_1(P_V1, n, Coord_1[0], masks[0]))
    return false;
  n = crossP(e_v1[0],e_v1[2]);

  // Maybe flip normal
  if (dot(n, e_v1[1]) > 0)
    n = neg(n);
  if (separating_plane_face_A_1(P_V1, n, Coord_1[1], masks[1]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 1))
    return false;
  n = crossP(e_v1[2], e_v1[1]);

  // Maybe flip normal
  if (dot(n, e_v1[0]) > 0)
      n = neg(n);
  if (separating_plane_face_A_1(P_V1, n, Coord_1[2], masks[2]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 2))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 1,2))
    return false;
  e_v1[4] = V1[3] - V1[1];
  e_v1[3] = V1[2] - V1[1];
  n = crossP(e_v1[3], e_v1[4]);

  // Maybe flip normal. Note the < since e_v1[0]=v1-v0.
  if (dot(n, e_v1[0]) < 0)
      n = neg(n);
  if (separating_plane_face_A_2(V1, V2, n, Coord_1[3], masks[3]))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 0, 3))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 1, 3))
    return false;
  if (separating_plane_edge_A(Coord_1, masks, 2, 3))
    return false;
  if ((masks[0] | masks[1] | masks[2] | masks[3] )!= 15)
    return true;

  // From now on, if there is a separating plane, it is parallel to a
  // face of b.
  std::vector<Point> P_V2(4);
  for (std::size_t i = 0; i < 4; ++i){
      Vector temp = V1[i] - V2[0];
      P_V2[i] = Point(temp.x(), temp.y(), temp.z());
  }

  e_v2[0] = V2[1] - V2[0];
  e_v2[1] = V2[2] - V2[0];
  e_v2[2] = V2[3] - V2[0];
  n = crossP(e_v2[1], e_v2[0]);

  // Maybe flip normal
  if (dot(n, e_v2[2])>0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  n=crossP(e_v2[0], e_v2[2]);

  // Maybe flip normal
  if (dot(n, e_v2[1]) > 0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  n = crossP(e_v2[2], e_v2[1]);

  // Maybe flip normal
  if (dot(n, e_v2[0]) > 0)
    n = neg(n);
  if (separating_plane_face_B_1(P_V2, n))
    return false;
  e_v2[4] = V2[3] - V2[1];
  e_v2[3] = V2[2] - V2[1];
  n = crossP(e_v2[3], e_v2[4]);

  // Maybe flip normal. Note the < since e_v2[0] = V2[1] - V2[0].
  if (dot(n, e_v2[0]) < 0)
    n = neg(n);
  if (separating_plane_face_B_2(V1, V2, n))
    return false;

  return true;
}
//-----------------------------------------------------------------------------
bool
separating_plane_face_A_1(const std::vector<Point>& pv1,
                          const Point& n,
                          std::vector<double>& coord,
                          int&  mask_edges)
{
  // Helper function for tetrahedron-tetrahedron collision test:
  // checks if plane pv1 is a separating plane. Stores local
  // coordinates and the mask bit mask_edges.

  mask_edges = 0;
  const int shifts[4] = {1, 2, 4, 8};

  for (std::size_t i = 0; i < 4; ++i)
  {
    coord[i] = dot(pv1[i],n);
    if (coord[i] > 0)
      mask_edges |= shifts[i];
  }

  return (mask_edges == 15);
}
//-----------------------------------------------------------------------------
bool
separating_plane_face_A_2(const std::vector<Point>& V1,
                          const std::vector<Point>& V2,
                          const Point& n,
                          std::vector<double>& coord,
                          int&  mask_edges)
{
  // Helper function for tetrahedron-tetrahedron collision test:
  // checks if plane v1,v2 is a separating plane. Stores local
  // coordinates and the mask bit mask_edges.

  mask_edges = 0;
  const int shifts[4] = {1, 2, 4, 8};

  for (std::size_t i = 0; i < 4; ++i)
  {
    coord[i] = dot((V2[i] - V1[1]),n);
    if (coord[i] > 0)
      mask_edges |= shifts[i];
  }

  return (mask_edges == 15);
}
//-----------------------------------------------------------------------------
bool
separating_plane_edge_A(
  const std::vector<std::vector<double>>& coord_1,
  const std::vector<int>& masks, int f0, int f1)
{
  // Helper function for tetrahedron-tetrahedron collision: checks if
  // edge is in the plane separating faces f0 and f1.

  const std::vector<double>& coord_f0 = coord_1[f0];
  const std::vector<double>& coord_f1 = coord_1[f1];

  int maskf0 = masks[f0];
  int maskf1 = masks[f1];

  if ((maskf0 | maskf1) != 15) // if there is a vertex of b
    return false; // included in (-,-) return false

  maskf0 &= (maskf0 ^ maskf1); // exclude the vertices in (+,+)
  maskf1 &= (maskf0 ^ maskf1);

  // edge 0: 0--1
  if ((maskf0 & 1) && // the vertex 0 of b is in (-,+)
      (maskf1 & 2)) // the vertex 1 of b is in (+,-)
    if ((coord_f0[1]*coord_f1[0] - coord_f0[0]*coord_f1[1]) > 0)
      // the edge of b (0,1) intersect (-,-) (see the paper)
      return false;

  if ((maskf0 & 2) &&
      (maskf1 & 1))
    if ((coord_f0[1]*coord_f1[0] - coord_f0[0]*coord_f1[1]) < 0)
      return false;

  // edge 1: 0--2
  if ((maskf0 & 1) &&
      (maskf1 & 4))
    if ((coord_f0[2]*coord_f1[0] - coord_f0[0]*coord_f1[2]) > 0)
      return false;

  if ((maskf0 & 4) &&
      (maskf1 & 1))
    if ((coord_f0[2]*coord_f1[0] - coord_f0[0]*coord_f1[2]) < 0)
      return false;

  // edge 2: 0--3
  if ((maskf0 & 1) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[0] - coord_f0[0]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 1))
    if ((coord_f0[3]*coord_f1[0] - coord_f0[0]*coord_f1[3]) < 0)
      return false;

  // edge 3: 1--2
  if ((maskf0 & 2) &&
      (maskf1 & 4))
    if ((coord_f0[2]*coord_f1[1] - coord_f0[1]*coord_f1[2]) > 0)
      return false;

  if ((maskf0 & 4) &&
      (maskf1 & 2))
    if ((coord_f0[2]*coord_f1[1] - coord_f0[1]*coord_f1[2]) < 0)
      return false;

  // edge 4: 1--3
  if ((maskf0 & 2) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[1] - coord_f0[1]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 2))
    if ((coord_f0[3]*coord_f1[1] - coord_f0[1]*coord_f1[3]) < 0)
      return false;

  // edge 5: 2--3
  if ((maskf0 & 4) &&
      (maskf1 & 8))
    if ((coord_f0[3]*coord_f1[2] - coord_f0[2]*coord_f1[3]) > 0)
      return false;

  if ((maskf0 & 8) &&
      (maskf1 & 4))
    if ((coord_f0[3]*coord_f1[2] - coord_f0[2]*coord_f1[3]) < 0)
      return false;

  // Now there exists a separating plane supported by the edge shared
  // by f0 and f1.
  return true;
}
//-----------------------------------------------------------------------------
