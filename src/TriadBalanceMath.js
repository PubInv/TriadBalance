// Copyright 2019, Robert L. Read
// This file is part of TriadBalance.
//
// TriadBalance is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TriadBalance is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TriadBalance.  If not, see <https://www.gnu.org/licenses/>.

"use strict";


// TODO: This needs to be scaled!!!
function mean(wtc) {
  return vec.scale(wtc.length,vec.addAll(wtc));
}

// Test that the point is vertically oriented, origin-centered equilateral triangle.
function isCenteredEquilateral(wtc) {
  let d0 = vec.dist(wtc[0],wtc[1]);
  let d1 = vec.dist(wtc[1],wtc[2]);
  let d2 = vec.dist(wtc[2],wtc[0]);

  return vec.near(1e-5,[0,0],mean(wtc))
  // Third point vertical..
    && ((wtc[2][0] == 0) && (wtc[2][1] > 0))
  // equilateral
    && vec.scalarNear(1e-5,d0,d1) && vec.scalarNear(1e-5,d1,d2) && vec.scalarNear(1e-5,d2,d1);
}



// This is a tiny set of routines inspired by vec-la-fp, which
// does not currently handle 3d vectors...a real hero would
// fully extend that package to support 2d and 3d vectors.

let v3Add = function v3Add(a,b) {
  return [a[0]+b[0],a[1]+b[1],a[2]+b[2]];
}
let v3Sub = function v3Sub(a,b) {
  return [a[0]-b[0],a[1]-b[1],a[2]-b[2]];
}
let v3ManhattanDistance = function v3ManhattanDistance(a,b) {
  return Math.abs(a[0]-b[0]) + Math.abs(a[1]-b[1]) + Math.abs(a[2]-b[2]);
}
let v3dist = function v3dist(a,b) {
  return Math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2);
}
let v3mag = function v3mag(v) {
  return Math.sqrt((v[0])**2 + (v[1])**2 + (v[2])**2);
}
let v3scale = function v3scale(sc,v) {
  return [sc*v[0],sc*v[1],sc*v[2]];
}
let v3normalize = function normalize(v) {
  return v3scale( 1/v3c.mag(v),v);
}

var v3c = {
  add: v3Add,
  sub: v3Sub,
  manhattan: v3ManhattanDistance,
  dist: v3dist,
  mag: v3mag,
  scale: v3scale,
  normalize: v3normalize
}



const SQRT3 = Math.sqrt(3);

function GetRayToLineSegmentIntersection(rayOrigin,rayDirection,point1,point2)
{
  // This code from here: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin
  // Note this routine seems to depend on the chirality of the points; possibly it only counts an approach from one side.
  const rdn = vec.norm(rayDirection);
  const v1 = vec.sub(rayOrigin,point1);
  const v2 = vec.sub(point2,point1);

  const v3 = [-rdn[1],rdn[0]];
  const dot = vec.dot(v2,v3);

  if (vec.scalarNear(1e-5,dot,0))
    return null;

  const t1 = vec.perpdot(v2,v1) / dot;
  const t2 = vec.dot(v1,v3) / dot;

  if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)) {
    return [vec.add(rayOrigin,vec.scale(t1,rdn)),t1];
  }
  return null;
}

// This is the fundamental math behind a TriadBalance Diagram.
// I will probably put this in LaTeX soon.
// A balance diagram is a way of choosing a unit vector
// of n dimensions (most usefully 3) from a n-gon.
// If n > 3, it is not possible to complete map all points.
// The fundamental math of a TriadBalance diagram is to convert
// between a point on the 2-dimensional n-gon (the "representation") from
// the "balance" -- a n-dimensional vector.
// Call the function that produces the repesentation from a unit vector
// "r" and the function that produces the balance vector from the respresentation "b".
// Desiderator of these functions are:
// We want them to be inversions of each other.
// We want the center of the representation to map to a balanced vector.
// A representation at the vertex is a vector having a single 1 and the rest zeros.
// As we change the angle betwen the origin and the point in a representation toward a vertex,
// the value of that vertex should increase.
// As we move along such a line, we should not change the relative proportion of the
// other values. (this is vague).
// It should be easy to compute and explain (at least geometrically.)
// I now believe this should allow a "norm" to be passed in
// as a function. I think the L1 norm is probably better than L2
// norm for some functions, but it can be optional.
// It is essential that this function be invertible.

function L1NORM(v) {
  let s = v3c.manhattan([0,0,0],v);
  return v3c.scale(1/s,v);
}
function L2NORM(v) {
  return v3c.normalize(v);
}
function L1LENGTH(v) {
  return v3c.manhattan([0,0,0],v);
}
function L2LENGTH(v) {
  return v3c.mag(v);
}

var L1 = [L1NORM,L1LENGTH];
var L2 = [L2NORM,L2LENGTH];


// Under assumption of an upward facing
// equilateral triangle, return the edge
// intersected by the ray from the origin to tp,
// without dependence on vector libraries.
// The return value is an array:
// empty iff tp is the origin
// containing two values if it is a vertex
// containing one value otherwise.
// The edges are numbered anticlockwise,
// with the zeroeth edge at the bottom.
// If two edges are returned, this routine
// returns them in sorted order.
function eqEdgeAlgebraically(wtc,p) {
  if (vec.scalarNear(1e-5,p[0],0)) {
    if (vec.scalarNear(1e-5,p[1],0)) {
      return [];
    } else if (p[1] > 0) {
      return [1,2];
    } else {
      return [0];
    }
  } else {
    let m = p[1]/p[0];
    let m1 = -SQRT3/3;
    let m2 = SQRT3/3;
    if ((p[0] > 0) && (vec.scalarNear(1e-5,m,m1))) return [0,1];
    if ((p[0] > 0) && (m > m1)) return [1];
    if ((p[0] > 0) && (m < m1)) return [0];

    if ((p[0] < 0) && (vec.scalarNear(1e-5,m,m2))) return [0,2];
    if ((p[0] < 0) && (m < m2)) return [2];
    if ((p[0] < 0) && (m > m2)) return [0];
  }
}


// Here we return the point on the edge
// where the ray from the origin to tp intersects
// the triangle.
// The math here is created on the assumption
// of the triangle being perfectly equilateral with
// the centroid at the origin. This allows us to
// solve simultaneous equations to find the points.
// This is an alternative to vector-based methods
// that are of course in the end simlar, but we take
// advantage of known slopes to make it faster and simpler.
function eqPointOnEdgeAlgebraically(wtc,tp) {
  // we should probably check equilaterality and orientation here
  let es = eqEdgeAlgebraically(wtc,tp);
  if (es.length == 0) return null; // tp is the origin
  if (es.length == 2) { // we hit a vertex, but which one?
    if ((es[0] == 0) && (es[1] == 1)) {
      return wtc[1];
    } else
    if ((es[0] == 0) && (es[1] == 2)) {
      return wtc[0];
    } else
    if ((es[0] == 1) && (es[1] == 2)) {
      return wtc[2];
    }
  } else { // now we do a case split
    let xp = tp[0];
    let yp = tp[1];
    let a = vec.dist(wtc[0],wtc[1]);
    let B = a * SQRT3/6;
    if (vec.scalarNear(1e-5,xp,0)) {
      return (yp > 0) ? wtc[2] : [0,-B];
    }
    let m = yp/xp;
    if (es[0] == 0) {
      return [-B/m,-B];
    } else if (es[0] == 1) {
      let y = a / (3 *(1/SQRT3 + 1/m));
      let x = y / m;
      return [x,y];
    } else if (es[0] == 2) {
      let y = a / (3 *(1/SQRT3 - 1/m));
      let x = y / m;
      return [x,y];
    }
  }
}


function getEdgeAndPoint(wtc,p) {

  // If we are centered, vertical, pointing up, and equilateral,
  // we can use the more efficient algorithm.
  if (isCenteredEquilateral(wtc))
  {
    // this may return two, but we can just take the first
    return [eqEdgeAlgebraically(wtc,p)[0],
            eqPointOnEdgeAlgebraically(wtc,p)];
  }
  else
  {
    var point_on_edge;
    var fe_idx = -1; // index of the first edge we intersect
    for(var i = 0; i < 3 && fe_idx < 0; i++) {
      var r = GetRayToLineSegmentIntersection([0,0],p,wtc[i],wtc[(i+1) % 3]);
      if (r != null) { // if null, the ray did not intersect the edge
        fe_idx = i;
        point_on_edge = r[0]; // The first comp. of return value is intersection
      }
    }
    return [fe_idx,point_on_edge];
  }
}


// tp is a point in the 2-dimensional triangle space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 3-vector in the attribute space
function TriadBalance2to3(p,wtc,LXnorm_and_length = L2) {
  let LXnormalize = LXnorm_and_length[0];

  if (vec.scalarNear(1e-5,vec.mag(p),0)) {
    return LXnormalize([1,1,1]);
  }

  // Now we want to do a linear interpolation of how far we are from an edge,
  // but also how far the projection to the edge is between the vertices.
  // We must first decide which edges the line from the orign to p intersects.
  // If it intersects two segments, then it is aimed at a vertex.
  let [fe_idx,point_on_edge] = getEdgeAndPoint(wtc,p);

  // now point_on_edge is a point on edge fe_idx.
  const total_distance_to_edge = vec.dist([0,0],point_on_edge);

  // If the point is outside the triangle, we clamp (truncate if needed)
  // it's length so that it is precisely on the edge.
  const pc = vec.clampMag(0,total_distance_to_edge,p);
  const distance_to_p_o_e = vec.dist(pc,point_on_edge);
  var ratio_p_to_edge =  distance_to_p_o_e/total_distance_to_edge;

  let bal = v3c.scale(ratio_p_to_edge,
                      LXnormalize([1,1,1]));

  // Now the remainder of the contribution
  // to the unit vector should come from the two
  // points on the edge, in linear proportion.
  // These coordinates are fe_idx and (fe_idx+1) % 3.
  const d1 = vec.dist(wtc[fe_idx],point_on_edge);
  const d2 = vec.dist(wtc[(fe_idx+1) % 3],point_on_edge);

  let vs = [0,0,0];
  vs[fe_idx] = d2;
  vs[(fe_idx+1) % 3] = d1;

  let imb = v3c.scale(1 - ratio_p_to_edge,LXnormalize(vs));

  return v3c.add(imb,bal);
}

// vec is a 3-vector in the attribute space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 2-vector in the triangle space
function invertTriadBalance2to3(v,wtc,LXnorm_and_length = L2) {
  let length = LXnorm_and_length[1];

  let min = Math.min(Math.min(v[0],v[1]),v[2]);

  let imb = [v[0] - min,v[1] - min,v[2] - min];
  let bal = v3c.sub(v,imb);
  // Now that we have balance, we need to compute it's length,
  // which is dependent on the norm we chose!

  let imb_r = length(imb);
  let bal_r = length(bal);

  // Now we have the ratios. We need to determine the direction.
  // This is a function of the imbalance vector. We can determine
  // which side we are on, and then compute our position along that
  // to determine a point on the triangle, and then multiply by the imb_r
  // to obtain the actual point.
  // At least one value of imb will be zero.
  var from_v,to_v,ratio;
  // the points are OPPOSITE the zero
  // ratio will be the ratio along the triangle edge
  // it requires a little thought to understand which
  // of the other points should be the "from_v" and the "to_v"
  // for the interpolation which occurs later.
  let s = imb[0] + imb[1] + imb[2]; // one of these is always zero.
  if (imb[0] == 0) {
    from_v = wtc[2];
    to_v = wtc[1];
    ratio = imb[1]/s;
  } else if (imb[1] == 0) {
    from_v = wtc[0];
    to_v = wtc[2];
    ratio = imb[2]/s;
  } else if (imb[2] == 0) {
    from_v = wtc[1];
    to_v = wtc[0];
    ratio = imb[0]/s;
  }

  // The point on the triangle is by construction
  // on one edge of the triangle.
  const onTriangle = vec.lerp(from_v,to_v,ratio);
  // now onTriangle is a point on the triangle
  // now, having found that we interpolate a ray
  // to it of length imb_r...
  return vec.lerp([0,0],onTriangle,imb_r);
}
