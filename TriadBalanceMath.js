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


// BRANCH specific: This is my attempt to use vec-la-fp
// TODO: Possibly this (the perp dot product) could be added to vec-la-fp
function perpdot(v1,v2)
{
  return (v1[0]*v2[1]) - (v1[1]*v2[0]);
}
// This is a candidate for entry.
function collinear(a,b,c) {
  var ar = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (b[1] - c[1]);
  return near(ar,0);
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

function GetRayToLineSegmentIntersection(rayOrigin,rayDirection,point1,point2)
{
  // This code from here: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin
  // Note this routine seems to depend on the chirality of the points; possibly it only counts an approach from one side.

  const rdn = vec.norm(rayDirection);
  const v1 = vec.sub(rayOrigin,point1);
  const v2 = vec.sub(point2,point1);  

  const v3 = [-rdn[1],rdn[0]];        
  const dot = vec.dot(v2,v3);

  if (Math.abs(dot) < 0.000001)
    return null;

  const t1 = perpdot(v2,v1) / dot;
  const t2 = vec.dot(v1,v3) / dot

  if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)) {
    return [vec.add(rayOrigin,vec.scale(t1,rayDirection)),t1];
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
  //  return v.length();
  return v3c.mag(v);
}

var L1 = [L1NORM,L1LENGTH];
var L2 = [L2NORM,L2LENGTH];

function near(x,y,e = 1e-4) {
  return Math.abs(x - y) < e;
}

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
function eqEdgeAlgebraically(wtc,xp,yp) {
  if (near(xp,0)) {
    if (near(yp,0)) {
      return [];
    } else if (yp > 0) {
      return [1,2];
    } else {
      return [0];
    }
  } else {
    let m = yp/xp;
    let m1 = -Math.sqrt(3)/3;
    let m2 = Math.sqrt(3)/3;
    if ((xp > 0) && (near(m,m1))) return [0,1];
    if ((xp > 0) && (m > m1)) return [1];
    if ((xp > 0) && (m < m1)) return [0];
    
    if ((xp < 0) && (near(m,m2))) return [0,2];    
    if ((xp < 0) && (m < m2)) return [2];
    if ((xp < 0) && (m > m2)) return [0];
  }
}

function test_eqEdgeAlebraically(wtc) {
  var rof = eqEdgeAlgebraically(wtc,0,0);
  console.assert(rof.length == 0);

  var r0 = eqEdgeAlgebraically(wtc,wtc[0][0],wtc[0][1]);
  console.assert(r0.length == 2);
  
  var r1 = eqEdgeAlgebraically(wtc,wtc[1][0],wtc[1][1]);
  console.assert(r1.length == 2);
  
  var r2 = eqEdgeAlgebraically(wtc,wtc[2][0],wtc[2][1]);
  console.assert(r2.length == 2);

  // slope of 1
  var rone = eqEdgeAlgebraically(wtc,1,1);
  console.assert(rone.length == 1);
  console.assert(rone[0] == 1);

  var rneg_one = eqEdgeAlgebraically(wtc,-1,1);
  console.assert(rneg_one.length == 1);
  console.assert(rneg_one[0] == 2);

  var rdown = eqEdgeAlgebraically(wtc,0.01,-1);
  console.assert(rdown.length == 1);
  console.assert(rdown[0] == 0);
  return true;
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
  let es = eqEdgeAlgebraically(wtc,tp[0],tp[1]);
  if (es.length == 0) return null; // tp is the origin
  if (es.length == 2) { // we hit a vertex, but which one?
    if ((es[0] == 0) && (es[1] == 1)) {
      return [wtc[1][0],wtc[1][1]];
    } else
    if ((es[0] == 0) && (es[1] == 2)) {
      return [wtc[0][0],wtc[0][1]];
    } else
    if ((es[0] == 1) && (es[1] == 2)) {
      return [wtc[2][0],wtc[2][1]];
    } 
  } else { // now we do a case split
    let xp = tp[0];
    let yp = tp[1];
    let a = vec.dist(wtc[0],wtc[1]);
    let B = a * Math.sqrt(3)/6;    
    if (near(xp,0)) {
      return (yp > 0) ? [wtc[2][0],wtc[2][1]] : [0,-B];
    }
    let m = yp/xp;
    if (es[0] == 0) {
      return [-B/m,-B];
    } else if (es[0] == 1) {
      let y = a / (3 *(1/Math.sqrt(3) + 1/m));
      let x = y / m;
      return [x,y];      
    } else if (es[0] == 2) {
      let y = a / (3 *(1/Math.sqrt(3) - 1/m));
      let x = y / m;
      return [x,y];            
    }
  }
}

function pointsNear(a,b) {
  return near(a[0],b[0]) && near(a[1],b[1]);
}

function collinear(a,b,c) {
  var ar = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (b[1] - c[1]);
  return near(ar,0);
}

function test_eqPointOnAlgebraically(wtc) {
  let ro = [0,0];
  var rof = eqPointOnEdgeAlgebraically(wtc,ro);
  console.assert(rof == null);

  const halfvert = vec.scale(1/2,[wtc[0][0],wtc[0][1]]);
  var r0 = eqPointOnEdgeAlgebraically(wtc,halfvert);
  console.assert(pointsNear(r0,[wtc[0][0],wtc[0][1]]));
  
  var r1 = eqPointOnEdgeAlgebraically(wtc,[wtc[1][0],wtc[1][1]]);
  console.assert(pointsNear(r1,[wtc[1][0],wtc[1][1]]));
  
  var r2 = eqPointOnEdgeAlgebraically(wtc,[wtc[2][0],wtc[2][1]]);
  console.assert(pointsNear(r2,[wtc[2][0],wtc[2][1]]));

  // slope of 1
  var one = [1,1];
  var rone = eqPointOnEdgeAlgebraically(wtc,one);
  console.assert(collinear(rone,[wtc[1][0],wtc[1][1]],[wtc[2][0],wtc[2][1]]));

  var neg_one = [-1,1];
  var rneg_one = eqPointOnEdgeAlgebraically(wtc,neg_one);
  console.assert(collinear(rneg_one,[wtc[0][0],wtc[0][1]],[wtc[2][0],wtc[2][1]]));  


  var down = [0.01,-1];
  var rdown = eqPointOnEdgeAlgebraically(wtc,down);
  console.assert(collinear(rdown,[wtc[0][0],wtc[0][1]],[wtc[1][0],wtc[1][1]]));  

  return true;
}


// tp is a point in the 2-dimensional triangle space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 3-vector in the attribute space
function TriadBalance2to3(p,wtc,LXnorm_and_length = L2) {
  let LXnormalize = LXnorm_and_length[0];
  if (near(vec.mag(p),0,1e-5)) {
    return LXnormalize([1,1,1]);    
  }
  //  const p = [tp[0],tp[1]];
  
  const origin = [0,0];
  // Now we want to do a linear interpolation of how far we are from an edge,
  // but also how far the projection to the edge is between the vertices.
  // We must first decide which edges the line from the orign to p intersects.
  // If it intersects two segments, then it is aimed at a vertex.

  let USE_ALGEBRAIC_STRATEGY = 1;
  var point_on_edge; 
  var fe_idx = -1; // index of the first edge we intersect
  if (USE_ALGEBRAIC_STRATEGY)
  {
    // this may return two, but we can just take the first
    fe_idx = eqEdgeAlgebraically(wtc,p[0],p[1])[0];
    point_on_edge = eqPointOnEdgeAlgebraically(wtc,p);
  }
  else
  {
    fe_idx = -1;
    for(var i = 0; i < 3 && fe_idx < 0; i++) {
      var r = GetRayToLineSegmentIntersection(origin,p,wtc[i],wtc[(i +1) % 3]);
      if (r != null) { // if null, the ray did not intersect the edge
        fe_idx = i;
        point_on_edge = r[0]; // The first comp. of return value is intersection
      }
    }
  }

  // now point_on_edge is a point on edge fe_idx.     

  const total_distance_to_edge = vec.dist(origin,point_on_edge);
  
  // If the point is outside the triangle, we clamp (truncate if needed)
  // it's length so that it is precisely on the edge.
  //  p.clampLength(0,total_distance_to_edge);
  function clampLength(min,max,v) {
    const d = vec.mag(v);
    if (d < min)
      return vec.scale(min/d,v);
    else if (d > max)
      return vec.scale(max/d,v);
    else
      return v;
  }
  const pc = clampLength(0,total_distance_to_edge,p);

  const distance_to_p_o_e = vec.dist(pc,point_on_edge);
  
  var ratio_p_to_edge =  distance_to_p_o_e/total_distance_to_edge;
  
  let ubal = LXnormalize([1,1,1]);

  let bal =v3c.scale(ratio_p_to_edge,ubal);

  // Now the remainder of the contribution
  // to the unit vector should come from the two
  // points on the edge, in linear proportion.
  // These coordinates are fe_idx and (fe_idx+1) % 3.
  const d1 = vec.dist([wtc[fe_idx][0],wtc[fe_idx][1]],point_on_edge);
  const d2 = vec.dist([wtc[(fe_idx+1) % 3][0],wtc[(fe_idx+1) % 3][1]],point_on_edge);  
  
  let vs = [0,0,0];
  vs[fe_idx] = d2;
  vs[(fe_idx+1) % 3] = d1;
  
  let imb = v3c.scale(1 - ratio_p_to_edge,LXnormalize(vs));  
  
  // now construct a balanced vector proportional
  // to the length from the edge to the point p towards the axis
  // so that this be a unit vector if p is the origin.
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
//  v.clone();
//  bal.sub(imb);
  // Now that we have balance, we need to compute it's length,
  // which is dependent on the norm we chose!

  let imb_r = length(imb);
  let bal_r = length(bal);
  console.assert(near(bal_r+imb_r,1));
  if (!near(bal_r+imb_r,1)) {
    debugger;
  }

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
  var s = imb[0] + imb[1] + imb[2]; // one of these is always zero.
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

  const from_vv = [from_v[0],from_v[1]];
  const to_vv = [to_v[0],to_v[1]];  
  // The point on the triangle is by construction
  // on one edge of the triangle.
  const onTriangle = vec.lerp(from_vv,to_vv,ratio);
  // now onTriangle is a point on the triangle
  // now, having found that we interpolate a ray
  // to it of length imb_r...
  const origin = [0,0];
  return vec.lerp(origin,onTriangle,imb_r);
}

function testGetRayToLineSegmentIntersection(wtc) {
  let ro = [0,0];
  let rd = [1,1];
  let p1 = [0,10];
  let p2 = [10,0];
  var r = GetRayToLineSegmentIntersection(ro,rd,p1,p2)[0];
  console.assert(r[0] == r[1]);
  var r = GetRayToLineSegmentIntersection(ro,rd,p2,p1)[0];
  console.assert(r[0] == r[1]);

  var rd1 = [94.1015625,-36.36328125];
  let c0 = [wtc[0][0],wtc[0][1]];
  let c1 = [wtc[1][0],wtc[1][1]];
  let c2 = [wtc[2][0],wtc[2][1]];        

  var r12 = GetRayToLineSegmentIntersection(ro,rd1,c1,c2);
  console.assert(r12 != null);
  return true;
}


function testTriadBalance2to3(wtc) {
  let p = [30000,30];    
  let pv = TriadBalance2to3(p,wtc,L1);
  console.assert(near(L1LENGTH(pv),1));
  let py = [0,wtc[2][1]];
  let pyv = TriadBalance2to3(py,wtc,L1);
  console.assert(near(L1LENGTH(pyv),1));
  return true;  
}

function testOriginAndVertices(wtc) {
  // The origin should return a perfectly balanced vector
  let o = [0,0];
  let ov = TriadBalance2to3(o,wtc,L1);
  console.assert(near(ov[0],1/3));
  console.assert(near(ov[1],1/3));
  console.assert(near(ov[2],1/3));    

  // A vertex should return a vector with a 1 in exactly 1 position
  {
    let p = [wtc[0][0],wtc[0][1]];
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv[0],1));
  }

  {
    let p = [wtc[1][0],wtc[1][1]];
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv[1],1));
  }

  {
    let p = [wtc[2][0],wtc[2][1]];
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv[2],1));
  }
  return true;  
}

function testInversion(wtc) {
  // This insures we are within the triangle
  let d = vec.dist(wtc[0],wtc[1]);  
  let p = [d/8,d/8];    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);
  // test length here
  let vpcx = [vp_inv[0],vp_inv[1]];
  let vpc = vec.sub(vpcx,p);
//  vpc.sub(p);
  console.assert(near(vec.mag(vpc),0));
  
  let py = [0,wtc[2][1]];
  var vpy = TriadBalance2to3(py,wtc,L1);    
  var vpy_inv = invertTriadBalance2to3(vpy,wtc,L1);
  let vpycx = [vpy_inv[0],vpy_inv[1]];
  let vpyc = vec.sub(vpycx,py);
//  vpyc.sub(py);
  console.assert(near(vec.mag(vpyc),0));
  return true;  
}

function testInversionOutside(wtc) {
  // This insures we are outside the triangle
  let d = vec.dist(wtc[0],wtc[1]);  
  let p = [d,d];    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);

  // we now want to test that the invesion has a slope of 1
  console.assert(near(vp_inv[1]/vp_inv[0],1));
  console.assert(vp_inv[1] > 0);  
  return true;  
}

function testInversionNegativeY(wtc) {
  // This insures we are within the triangle  
  let d = vec.dist(wtc[0],wtc[1]);
  let p = [0,-d/4];    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);
  // test length here
  var vpc = vec.sub([vp_inv[0],vp_inv[1]],p);
  console.assert(near(vec.mag(vpc),0));
  return true;  
}

// Test via a circle completely within the triangle,
// thus exercising all angles (instigated by a bug.)
function testInversionWithACircle(wtc,NORM) {

  let wtcp = [wtc[1][0],wtc[1][1]];
  let radius = vec.mag(wtcp)/3;    
  let n = 10;
  let one_thirteenth = 2 * Math.PI / 13;
  for(var i = 0; i < n; i++) {
    let x = Math.sin(i*one_thirteenth);
    let y = Math.cos(i*one_thirteenth);
    // now x,y is a point on a circle within the trinagle
    // we will make sure the balance function inverts
    // to give the function back to us.
    {
      let p = [x,y];
      vec.scale(radius,p);
      var vp = TriadBalance2to3(p,wtc,NORM);
      var vp_inv = invertTriadBalance2to3(vp,wtc,NORM);
      var vpyc = vec.sub([vp_inv[0],vp_inv[1]],p);
      console.assert(near(vec.mag(vpyc),0));
      if (!near(vec.mag(vpyc),0)) {
        debugger;
      }
    }
  }
  return true;  
}

// wtc is the WORLD_TRIANGLE_COORDS, an array of three Vector2 objects.

function testEquilateralFunctions(wtc) {
  var allTrue = 1;
  allTrue &= test_eqEdgeAlebraically(wtc);
  allTrue &= test_eqPointOnAlgebraically(wtc);
  return allTrue;
}
function testAllTriadBalance(wtc) {
  var allTrue = 1;
  allTrue &= testEquilateralFunctions(wtc);
  allTrue &= testOriginAndVertices(wtc);    
  allTrue &= testGetRayToLineSegmentIntersection(wtc);
  allTrue &= testTriadBalance2to3(wtc);
  allTrue &= testInversion(wtc);
  allTrue &= testInversionNegativeY(wtc);
  allTrue &= testInversionOutside(wtc);
  let wtcp = wtc[1];
  

  let small_radius = vec.mag(wtcp)/3;    
  allTrue &= testInversionWithACircle(wtc,L1,small_radius);
  allTrue &= testInversionWithACircle(wtc,L2,small_radius);    
  let large_radius = vec.mag(wtcp)*3;    
  allTrue &= testInversionWithACircle(wtc,L1,large_radius);
  allTrue &= testInversionWithACircle(wtc,L2,large_radius);
  return allTrue;
}


