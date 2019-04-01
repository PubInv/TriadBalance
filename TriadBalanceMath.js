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
  let r = new THREE.Vector3(0,0,0).manhattanDistanceTo(v);
  return v.divideScalar(r);
}
function L2NORM(v) {
  return v.normalize();
}
function L1LENGTH(v) {
  let r = new THREE.Vector3(0,0,0).manhattanDistanceTo(v);
  return r;
}
function L2LENGTH(v) {
  return v.length();
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

  var r0 = eqEdgeAlgebraically(wtc,wtc[0].x,wtc[0].y);
  console.assert(r0.length == 2);
  
  var r1 = eqEdgeAlgebraically(wtc,wtc[1].x,wtc[1].y);
  console.assert(r1.length == 2);
  
  var r2 = eqEdgeAlgebraically(wtc,wtc[2].x,wtc[2].y);
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
      return [wtc[1].x,wtc[1].y];
    } else
    if ((es[0] == 0) && (es[1] == 2)) {
      return [wtc[0].x,wtc[0].y];
    } else
    if ((es[0] == 1) && (es[1] == 2)) {
      return [wtc[2].x,wtc[2].y];
    } 
  } else { // now we do a case split
    let xp = tp[0];
    let yp = tp[1];
    let a = wtc[0].distanceTo(wtc[1]);
    let B = a * Math.sqrt(3)/6;    
    if (near(xp,0)) {
      return (yp > 0) ? [wtc[2].x,wtc[2].y] : [0,-B];
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

  const halfvert = vec.scale(1/2,[wtc[0].x,wtc[0].y]);
  var r0 = eqPointOnEdgeAlgebraically(wtc,halfvert);
  console.assert(pointsNear(r0,[wtc[0].x,wtc[0].y]));
  
  var r1 = eqPointOnEdgeAlgebraically(wtc,[wtc[1].x,wtc[1].y]);
  console.assert(pointsNear(r1,[wtc[1].x,wtc[1].y]));
  
  var r2 = eqPointOnEdgeAlgebraically(wtc,[wtc[2].x,wtc[2].y]);
  console.assert(pointsNear(r2,[wtc[2].x,wtc[2].y]));

  // slope of 1
  var one = [1,1];
  var rone = eqPointOnEdgeAlgebraically(wtc,one);
  console.assert(collinear(rone,[wtc[1].x,wtc[1].y],[wtc[2].x,wtc[2].y]));

  var neg_one = [-1,1];
  var rneg_one = eqPointOnEdgeAlgebraically(wtc,neg_one);
  console.assert(collinear(rneg_one,[wtc[0].x,wtc[0].y],[wtc[2].x,wtc[2].y]));  


  var down = [0.01,-1];
  var rdown = eqPointOnEdgeAlgebraically(wtc,down);
  console.assert(collinear(rdown,[wtc[0].x,wtc[0].y],[wtc[1].x,wtc[1].y]));  

  return true;
}


// tp is a point in the 2-dimensional triangle space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 3-vector in the attribute space
function TriadBalance2to3(tp,wtc,LXnorm_and_length = L2) {
  let LXnormalize = LXnorm_and_length[0];
  if (near(tp.lengthSq(),0,1e-5)) {
    return LXnormalize(new THREE.Vector3(1,1,1));
  }
//  let p = new THREE.Vector2(tp.x,tp.y);
  const p = [tp.x,tp.y];
  
//  let origin = new THREE.Vector2(0,0);
  const origin = [0,0];
  // Now we want to do a linear interpolation of how far we are from an edge,
  // but also how far the projection to the edge is between the vertices.
  // We must first decide which edges the line from the orign to p intersects.
  // If it intersects two segments, then it is aimed at a vertex.

  let USE_ALGEBRAIC_STRATEGY = 1;
  var point_on_edge; 
  var fe_idx = -1; // index of the first edge we intersect
  console.log("XXXX",tp);
  if (USE_ALGEBRAIC_STRATEGY)
  {
    // this may return two, but we can just take the first
    fe_idx = eqEdgeAlgebraically(wtc,tp.x,tp.y)[0];
    point_on_edge = eqPointOnEdgeAlgebraically(wtc,[tp.x,tp.y]);
//    point_on_edge = new THREE.Vector2(poe[0],poe[1]);
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

  //  let total_distance_to_edge = origin.distanceTo(point_on_edge);
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

  //  let distance_to_p_o_e = p.distanceTo(point_on_edge);
  const distance_to_p_o_e = vec.dist(pc,point_on_edge);
  
  var ratio_p_to_edge =  distance_to_p_o_e/total_distance_to_edge;
  
  let bal = LXnormalize(new THREE.Vector3(1,1,1));

  bal.multiplyScalar(ratio_p_to_edge);

  // Now the remainder of the contribution
  // to the unit vector should come from the two
  // points on the edge, in linear proportion.
  // These coordinates are fe_idx and (fe_idx+1) % 3.
  //  var d1 = wtc[fe_idx].distanceTo(point_on_edge);
  const d1 = vec.dist([wtc[fe_idx].x,wtc[fe_idx].y],point_on_edge);
  //  var d2 = wtc[(fe_idx+1) % 3].distanceTo(point_on_edge);
  const d2 = vec.dist([wtc[(fe_idx+1) % 3].x,wtc[(fe_idx+1) % 3].y],point_on_edge);  
  
  let vs = [0,0,0];
  vs[fe_idx] = d2;
  vs[(fe_idx+1) % 3] = d1;
  
  let imb = LXnormalize(new THREE.Vector3(vs[0],vs[1],vs[2]));
  imb.multiplyScalar(1 - ratio_p_to_edge);
  
  // now construct a balanced vector proportional
  // to the length from the edge to the point p towards the axis
  // so that this be a unit vector if p is the origin.
  return new THREE.Vector3().add(imb).add(bal);
}

// vec is a 3-vector in the attribute space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 2-vector in the triangle space
function invertTriadBalance2to3(v,wtc,LXnorm_and_length = L2) {
  let length = LXnorm_and_length[1];
  let min = Math.min(Math.min(v.x,v.y),v.z);
  let imb = new THREE.Vector3(v.x - min, v.y - min, v.z - min);
  let bal = v.clone();
  bal.sub(imb);
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
  var s = imb.x + imb.y + imb.z; // one of these is always zero.
  if (imb.x == 0) {
    from_v = wtc[2];
    to_v = wtc[1];
    ratio = imb.y/s;
  } else if (imb.y == 0) {
    from_v = wtc[0];
    to_v = wtc[2];
    ratio = imb.z/s;        
  } else if (imb.z == 0) {
    from_v = wtc[1];
    to_v = wtc[0];
    ratio = imb.x/s;        
  }

  const from_vv = [from_v.x,from_v.y];
  const to_vv = [to_v.x,to_v.y];  
  // The point on the triangle is by construction
  // on one edge of the triangle.
//  var onTriangle = new THREE.Vector2();
  //  onTriangle.lerpVectors(from_v,to_v,ratio);
  const onTriangle = vec.lerp(from_vv,to_vv,ratio);
  // now onTriangle is a point on the triangle
  // now, having found that we interpolate a ray
  // to it of length imb_r...
//  let origin = new THREE.Vector2(0,0);
//  let inversion = new THREE.Vector2();
//  inversion.lerpVectors(origin,
//                        onTriangle,
  //                        imb_r);
  const origin = [0,0];
  const inversion = vec.lerp(origin,onTriangle,imb_r);
  
  return new THREE.Vector2(inversion[0],inversion[1]);
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
  let c0 = [wtc[0].x,wtc[0].y];
  let c1 = [wtc[1].x,wtc[1].y];
  let c2 = [wtc[2].x,wtc[2].y];        

  var r12 = GetRayToLineSegmentIntersection(ro,rd1,c1,c2);
  console.assert(r12 != null);
  return true;
}


function testTriadBalance2to3(wtc) {
  let p = new THREE.Vector2(30000,30);    
  let pv = TriadBalance2to3(p,wtc,L1);
  console.assert(near(L1LENGTH(pv),1));
  let py = new THREE.Vector2(0,wtc[2].y);
  let pyv = TriadBalance2to3(py,wtc,L1);
  console.assert(near(L1LENGTH(pyv),1));
  return true;  
}

function testOriginAndVertices(wtc) {
  // The origin should return a perfectly balanced vector
  let o = new THREE.Vector2(0,0);
  let ov = TriadBalance2to3(o,wtc,L1);
  console.assert(near(ov.x,1/3));
  console.assert(near(ov.y,1/3));
  console.assert(near(ov.z,1/3));    

  // A vertex should return a vector with a 1 in exactly 1 position
  {
    let p = new THREE.Vector2(wtc[0].x,wtc[0].y);
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv.x,1));
  }

  {
    let p = new THREE.Vector2(wtc[1].x,wtc[1].y);
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv.y,1));
  }

  {
    let p = new THREE.Vector2(wtc[2].x,wtc[2].y);
    let pv = TriadBalance2to3(p,wtc,L1);
    console.assert(near(pv.z,1));
  }
  return true;  
}

function testInversion(wtc) {
  // This insures we are within the triangle
  let d = wtc[0].distanceTo(wtc[1]);  
  let p = new THREE.Vector2(d/8,d/8);    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);
  // test length here
  var vpc = vp_inv.clone();
  vpc.sub(p);
  console.assert(near(vpc.length(),0));
  
  let py = new THREE.Vector2(0,wtc[2].y);
  var vpy = TriadBalance2to3(py,wtc,L1);    
  var vpy_inv = invertTriadBalance2to3(vpy,wtc,L1);
  
  var vpyc = vpy_inv.clone();
  vpyc.sub(py);
  console.assert(near(vpyc.length(),0));
  return true;  
}

function testInversionOutside(wtc) {
  // This insures we are outside the triangle
  let d = wtc[0].distanceTo(wtc[1]);  
  let p = new THREE.Vector2(d,d);    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);

  // we now want to test that the invesion has a slope of 1
  console.assert(near(vp_inv.y/vp_inv.x,1));
  console.assert(vp_inv.y > 0);  
  return true;  
}

function testInversionNegativeY(wtc) {
  // This insures we are within the triangle  
  let d = wtc[0].distanceTo(wtc[1]);
  let p = new THREE.Vector2(0,-d/4);    
  var vp = TriadBalance2to3(p,wtc,L1);
  var vp_inv = invertTriadBalance2to3(vp,wtc,L1);
  // test length here
  var vpc = vp_inv.clone();
  vpc.sub(p);
  console.assert(near(vpc.length(),0));
  return true;  
}

// Test via a circle completely within the triangle,
// thus exercising all angles (instigated by a bug.)
function testInversionWithACircle(wtc,NORM) {

  let wtcp = new THREE.Vector2(wtc[1].x,
                               wtc[1].y);
  let radius = wtcp.length()/3;    
  let n = 10;
  let one_thirteenth = 2 * Math.PI / 13;
  for(var i = 0; i < n; i++) {
    let x = Math.sin(i*one_thirteenth);
    let y = Math.cos(i*one_thirteenth);
    // now x,y is a point on a circle within the trinagle
    // we will make sure the balance function inverts
    // to give the function back to us.
    {
      let p = new THREE.Vector2(x,y);
      p.multiplyScalar(radius);        
      var vp = TriadBalance2to3(p,wtc,NORM);
      var vp_inv = invertTriadBalance2to3(vp,wtc,NORM);
      var vpyc = vp_inv.clone();
      vpyc.sub(p);
      console.assert(near(vpyc.length(),0));
      if (!near(vpyc.length(),0)) {
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

  let wtcp = new THREE.Vector2(wtc[1].x,
                               wtc[1].y);
  
  let small_radius = wtcp.length()/3;    
  allTrue &= testInversionWithACircle(wtc,L1,small_radius);
  allTrue &= testInversionWithACircle(wtc,L2,small_radius);    
  let large_radius = wtcp.length()*3;    
  allTrue &= testInversionWithACircle(wtc,L1,large_radius);
  allTrue &= testInversionWithACircle(wtc,L2,large_radius);
  return allTrue;
}


