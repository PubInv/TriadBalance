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

// It is extremely valuable to be able to compute the
// inverse of the TriadBalance algorithm for testing,
// although this is not quite a true invesion because
// points outside the triangle are brought to the triangle
// edge.

"use strict";

function computeFinverseF(wtc,norm,p) {
   return invertTriadBalance2to3(
     TriadBalance2to3(p,wtc,norm),
     wtc,
     norm);
}

function test_eqPointOnAlgebraically(wtc) {
  var rof = eqPointOnEdgeAlgebraically(wtc,[0,0]);
  console.assert(rof == null);

  const halfvert = vec.scale(1/2,wtc[0]);
  var r0 = eqPointOnEdgeAlgebraically(wtc,halfvert);
  console.assert(vec.near(1e-5,r0,wtc[0]));

  var r1 = eqPointOnEdgeAlgebraically(wtc,wtc[1]);
  console.assert(vec.near(1e-5,r1,wtc[1]));

  var r2 = eqPointOnEdgeAlgebraically(wtc,wtc[2]);
  console.assert(vec.near(1e-5,r2,wtc[2]));

  // slope of 1
  var rone = eqPointOnEdgeAlgebraically(wtc,[1,1]);
  console.assert(vec.colinear(wtc[1],wtc[2],rone));

  // slope of -1
  var rneg_one = eqPointOnEdgeAlgebraically(wtc,[-1,1]);
  console.assert(vec.colinear(wtc[0],wtc[2],rneg_one));

  // almost straight down
  var rdown = eqPointOnEdgeAlgebraically(wtc,[0.01,-1]);
  console.assert(vec.colinear(rdown,wtc[0],wtc[1]));

  return true;
}

function test_eqEdgeAlebraically(wtc) {
  var rof = eqEdgeAlgebraically(wtc,[0,0]);
  console.assert(rof.length == 0);

  var r0 = eqEdgeAlgebraically(wtc,wtc[0]);
  console.assert(r0.length == 2);

  var r1 = eqEdgeAlgebraically(wtc,wtc[1]);
  console.assert(r1.length == 2);

  var r2 = eqEdgeAlgebraically(wtc,wtc[2]);
  console.assert(r2.length == 2);

  // slope of 1
  var rone = eqEdgeAlgebraically(wtc,[1,1]);
  console.assert(rone.length == 1);
  console.assert(rone[0] == 1);

  var rneg_one = eqEdgeAlgebraically(wtc,[-1,1]);
  console.assert(rneg_one.length == 1);
  console.assert(rneg_one[0] == 2);

  var rdown = eqEdgeAlgebraically(wtc,[0.01,-1]);
  console.assert(rdown.length == 1);
  console.assert(rdown[0] == 0);
  return true;
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
  let c0 = wtc[0];
  let c1 = wtc[1];
  let c2 = wtc[2];

  var r12 = GetRayToLineSegmentIntersection(ro,rd1,c1,c2);
  console.assert(r12 != null);

  {
    var down = GetRayToLineSegmentIntersection(ro,[0,-4999],c0,c1);
    console.assert(vec.scalarNear(1e-5,down[0][1],wtc[0][1]));
  }
  return true;
}


function testTriadBalance2to3(wtc) {
  let pv = TriadBalance2to3([30000,30],wtc,L1);
  console.assert(vec.scalarNear(1e-5,L1LENGTH(pv),1));

  let pyv = TriadBalance2to3([0,wtc[2][1]],wtc,L1);
  console.assert(vec.scalarNear(1e-5,L1LENGTH(pyv),1));
  return true;
}

function testOriginAndVertices(wtc) {
  // The origin should return a perfectly balanced vector
  let ov = TriadBalance2to3([0,0],wtc,L1);
  console.assert(vec.scalarNear(1e-5,ov[0],1/3));
  console.assert(vec.scalarNear(1e-5,ov[1],1/3));
  console.assert(vec.scalarNear(1e-5,ov[2],1/3));

  // A vertex should return a vector with a 1 in exactly 1 position
  console.assert(vec.scalarNear(1e-5,TriadBalance2to3(wtc[0],wtc,L1)[0],1));
  console.assert(vec.scalarNear(1e-5,TriadBalance2to3(wtc[1],wtc,L1)[1],1));
  console.assert(vec.scalarNear(1e-5,TriadBalance2to3(wtc[2],wtc,L1)[2],1));

  return true;
}

// It is extremely valuable to be able to compute the
// inverse of the TriadBalance algorithm for testing,
// although this is not quite a true invesion because
// points outside the triangle are brought to the triangle
// edge.
function computeFinverseF(wtc,norm,p) {
   return invertTriadBalance2to3(
     TriadBalance2to3(p,wtc,norm),
     wtc,
     norm);
}

function testInversion(wtc) {
  // This insures we are within the triangle
  let d = vec.dist(wtc[0],wtc[1]);
  let p = [d/8,d/8];

  let vpc = vec.sub(
    computeFinverseF(wtc,L1,p),
    p);
  console.assert(vec.scalarNear(1e-5,vec.mag(vpc),0));

  let py = [0,wtc[2][1]];
  let vpyc = vec.sub(
    computeFinverseF(wtc,L1,py),
    py);
  console.assert(vec.scalarNear(1e-5,vec.mag(vpyc),0));
  return true;
}

function testInversionOutside(wtc) {
  // This insures we are outside the triangle
  let d = vec.dist(wtc[0],wtc[1]);

  var vp_inv = computeFinverseF(wtc,L1,[d,d]);
  // we now want to test that the invesion has a slope of 1
  console.assert(vec.scalarNear(1e-5,vp_inv[1]/vp_inv[0],1));
  console.assert(vp_inv[1] > 0);
  return true;
}

function testInversionNegativeY(wtc) {
  // This insures we are within the triangle
  let d = vec.dist(wtc[0],wtc[1]);
  let p = [0,-d/4];
  var vp_inv = computeFinverseF(wtc,L1,p);
  // test length here
  var vpc = vec.sub(vp_inv,p);
  console.assert(vec.scalarNear(1e-5,vec.mag(vpc),0));
  return true;
}

// Test via a circle completely within the triangle,
// thus exercising all angles (instigated by a bug.)
function testInversionWithACircle(wtc,NORM) {

  let wtcp = wtc[1];
  let radius = vec.mag(wtcp)/3;
  let n = 13;
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
      var vp_inv =  computeFinverseF(wtc,NORM,p);
      var vpyc = vec.sub(vp_inv,p);
      console.assert(vec.scalarNear(1e-5,vec.mag(vpyc),0));
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
function testAllTriadBalance(upward,wtc) {
  var allTrue = 1;
  if (upward) {
    allTrue &= testEquilateralFunctions(wtc);
    allTrue &= testGetRayToLineSegmentIntersection(wtc);
  }
  allTrue &= testOriginAndVertices(wtc);
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
