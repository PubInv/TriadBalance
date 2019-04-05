(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.triadbtest = f()}})(function(){var define,module,exports;return (function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});

function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

// curry :: (a -> b -> ... -> n) -> (a -> b) -> (b -> ...) -> (... -> n)
var curry = function curry(fn) {
  var curried = function curried() {
    for (var _len = arguments.length, args = Array(_len), _key = 0; _key < _len; _key++) {
      args[_key] = arguments[_key];
    }

    if (args.length >= fn.length) {
      return fn.apply(undefined, args);
    }
    return function () {
      for (var _len2 = arguments.length, argsNext = Array(_len2), _key2 = 0; _key2 < _len2; _key2++) {
        argsNext[_key2] = arguments[_key2];
      }

      return curried.apply(undefined, args.concat(argsNext));
    };
  };
  return curried;
};

// pipe :: (a -> b) -> (b -> ...) -> (... -> n)
var pipe = function pipe(fn1) {
  for (var _len3 = arguments.length, functions = Array(_len3 > 1 ? _len3 - 1 : 0), _key3 = 1; _key3 < _len3; _key3++) {
    functions[_key3 - 1] = arguments[_key3];
  }

  return function () {
    return functions.reduce(function (acc, fn) {
      return fn(acc);
    }, fn1.apply(undefined, arguments));
  };
};

// compose :: (... -> n) -> (b -> ...) -> (a -> b)
var compose = function compose() {
  for (var _len4 = arguments.length, functions = Array(_len4), _key4 = 0; _key4 < _len4; _key4++) {
    functions[_key4] = arguments[_key4];
  }

  return pipe.apply(undefined, _toConsumableArray(functions.reverse()));
};

// vAdd :: Vector -> Vector -> Vector
var vAdd = curry(function (v, v2) {
  return [v[0] + v2[0], v[1] + v2[1]];
});

// vAdd3 :: Vector -> Vector -> Vector -> Vector
var vAdd3 = curry(function (v, v2, v3) {
  return [v[0] + v2[0] + v3[0], v[1] + v2[1] + v3[1]];
});

// vAddAll :: [Vector] -> Vector
var vAddAll = function vAddAll(vs) {
  return vs.reduce(vAdd, [0, 0]);
};

// vSub :: Vector -> Vector -> Vector
var vSub = curry(function (v, v2) {
  return [v[0] - v2[0], v[1] - v2[1]];
});

// vSub3 :: Vector -> Vector -> Vector -> Vector
var vSub3 = curry(function (v, v2, v3) {
  return [v[0] - v2[0] - v3[0], v[1] - v2[1] - v3[1]];
});

// vSubAll :: [Vector] -> Vector
var vSubAll = function vSubAll(vs) {
  return vs.slice(1).reduce(vSub, vs.slice(0, 1)[0]);
};

// vMag :: Vector -> Number
var vMag = function vMag(v) {
  return Math.sqrt(v[0] * v[0] + v[1] * v[1]);
};

// vNormal :: Vector -> Vector
var vNormal = function vNormal(v) {
  return [-v[1], v[0]];
};

// vScale :: Number -> Vector
var vScale = curry(function (sc, v) {
  return [v[0] * sc, v[1] * sc];
});

// vTowards :: Number -> Vector -> Vector -> Vector
var vTowards = curry(function (t, v1, v2) {
  var d = vSub(v2, v1);
  var sc = vMag(d) * t;
  return vAdd(v1, vScale(sc, vNorm(d)));
});

// vLerp :: Vector -> Vector -> Number -> Vector
var vLerp = curry(function (v1, v2, t) {
  return vTowards(t, v1, v2);
});

// vScalarNear :: Number -> Number -> Number -> bool
var vScalarNear = curry(function (e, a, b) {
  return Math.abs(a - b) < e;
});

// vNear :: Number -> Vector -> Vector -> bool
var vNear = curry(function (e, a, b) {
  return vScalarNear(e, a[0], b[0]) && vScalarNear(e, a[1], b[1]);
});

// vClampMag :: Number -> Number -> Vector -> Vector
var vClampMag = curry(function (min, max, v) {
  var d = vec.mag(v);
  if (d < min) return vec.scale(min / d, v);else if (d > max) return vec.scale(max / d, v);else return v;
});

// vNorm :: Vector -> Vector
var vNorm = function vNorm(v) {
  var mag = vMag(v);
  return [v[0] / mag, v[1] / mag];
};

// mId :: Matrix
var mId = Object.freeze([1, 0, 0, 0, 1, 0, 0, 0, 1]);

// vCreateMatrix :: Number -> Number -> Number -> Number -> Number -> Number -> Matrix
var vCreateMatrix = function vCreateMatrix() {
  var a = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 1;
  var b = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : 0;
  var c = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 0;
  var d = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 1;
  var tx = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : 0;
  var ty = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : 0;
  return [a, c, tx, b, d, ty, 0, 0, 1];
};

// vTransform :: Matrix -> Vector -> Vector
var vTransform = curry(function (m, v) {
  return [v[0] * m[0] + v[1] * m[1] + m[2], v[0] * m[3] + v[1] * m[4] + m[5]];
});

// mCompose :: Matrix -> Matrix -> Matrix
var mCompose = curry(function (m, m2) {
  return [m[0] * m2[0] + m[1] * m2[3] + m[2] * m2[6], m[0] * m2[1] + m[1] * m2[4] + m[2] * m2[7], m[0] * m2[2] + m[1] * m2[5] + m[2] * m2[8], m[3] * m2[0] + m[4] * m2[3] + m[5] * m2[6], m[3] * m2[1] + m[4] * m2[4] + m[5] * m2[7], m[3] * m2[2] + m[4] * m2[5] + m[5] * m2[8], m[6] * m2[0] + m[7] * m2[3] + m[8] * m2[6], m[6] * m2[1] + m[7] * m2[4] + m[8] * m2[7], m[6] * m2[2] + m[7] * m2[5] + m[8] * m2[8]];
});

// mRotate :: Number -> Matrix -> Matrix
var mRotate = function mRotate(a) {
  return mCompose([Math.cos(a), -Math.sin(a), 0, Math.sin(a), Math.cos(a), 0, 0, 0, 1]);
};

// mTranslate :: Vector -> Matrix -> Matrix
var mTranslate = function mTranslate(v) {
  return mCompose([1, 0, v[0], 0, 1, v[1], 0, 0, 1]);
};

// mScale :: Vector -> Matrix -> Matrix
var mScale = function mScale(v) {
  return mCompose([v[0], 0, 0, 0, v[1], 0, 0, 0, 1]);
};

// mShear :: Vector -> Matrix -> Matrix
var mShear = function mShear(v) {
  return mCompose([1, v[0], 0, v[1], 1, 0, 0, 0, 1]);
};

// vRotate :: Number -> Vector -> Vector
var vRotate = curry(function (a, v) {
  return [v[0] * Math.cos(a) - v[1] * Math.sin(a), v[0] * Math.sin(a) + v[1] * Math.cos(a)];
});

// vRotatePointAround :: Number -> Vector -> Vector -> Vector
var vRotatePointAround = curry(function (a, cp, v) {
  var v2 = vSub(v, cp);
  return vAdd(cp, [v2[0] * Math.cos(a) - v2[1] * Math.sin(a), v2[0] * Math.sin(a) + v2[1] * Math.cos(a)]);
});

// vMidpoint :: Vector -> Vector -> Vector
var vMidpoint = curry(function (v, v2) {
  return vScale(0.5, vAdd(v, v2));
});

// vAngle :: Number -> Vector
var vAngle = function vAngle(a) {
  return [Math.cos(a), Math.sin(a)];
};

// vAlongAngle :: Number -> Number -> Vector
var vAlongAngle = curry(function (a, r, v) {
  return compose(vAdd(v), vScale(r), vAngle)(a);
});

// vFastDist :: Vector -> Vector -> Number
var vFastDist = curry(function (v, v2) {
  return Math.pow(v2[0] - v[0], 2) + Math.pow(v2[1] - v[1], 2);
});

// vDist :: Vector -> Vector -> Number
var vDist = curry(function (v, v2) {
  return Math.hypot(v2[0] - v[0], v2[1] - v[1]);
});

// vDot :: Vector -> Vector -> Number
var vDot = curry(function (v, v2) {
  return v[0] * v2[0] + v[1] * v2[1];
});

// vPerpDot :: Vector -> Vector -> Number
var vPerpDot = curry(function (v, v2) {
  return v[0] * v2[1] - v[1] * v2[0];
});

// vTriangleArea :: Vector -> Vector -> Vector -> Number
var vTriangleArea = curry(function (a, b, c) {
  return ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2;
});

// vColinear :: Vector -> Vector -> Vector -> bool
var vColinear = curry(function (v0, v1, v2) {
  return vScalarNear(1e-4, vTriangleArea(v0, v1, v2), 0);
});

// vDet :: Matrix -> Number
var vDet = function vDet(m) {
  return m[0] * m[4] - m[3] * m[1];
};

var vec = {
  add: vAdd,
  add3: vAdd3,
  addAll: vAddAll,
  sub: vSub,
  sub3: vSub3,
  subAll: vSubAll,
  mag: vMag,
  normal: vNormal,
  scale: vScale,
  towards: vTowards,
  lerp: vLerp,
  scalarNear: vScalarNear,
  near: vNear,
  clampMag: vClampMag,
  norm: vNorm,
  mId: mId,
  createMatrix: vCreateMatrix,
  transform: vTransform,
  mCompose: mCompose,
  mRotate: mRotate,
  mTranslate: mTranslate,
  mScale: mScale,
  mShear: mShear,
  rotate: vRotate,
  rotatePointAround: vRotatePointAround,
  midpoint: vMidpoint,
  angle: vAngle,
  alongAngle: vAlongAngle,
  fastDist: vFastDist,
  dist: vDist,
  dot: vDot,
  perpdot: vPerpDot,
  triangleArea: vTriangleArea,
  colinear: vColinear,
  det: vDet
};

/* start exports */
exports.default = vec;
exports.vec = vec;
exports.vAdd = vAdd;
exports.vAdd3 = vAdd3;
exports.vAddAll = vAddAll;
exports.vSub = vSub;
exports.vSub3 = vSub3;
exports.vSubAll = vSubAll;
exports.vMag = vMag;
exports.vNormal = vNormal;
exports.vScale = vScale;
exports.vTowards = vTowards;
exports.vLerp = vLerp;
exports.vScalarNear = vScalarNear;
exports.vNear = vNear;
exports.vClampMag = vClampMag;
exports.vNorm = vNorm;
exports.mId = mId;
exports.vCreateMatrix = vCreateMatrix;
exports.vTransform = vTransform;
exports.mCompose = mCompose;
exports.mRotate = mRotate;
exports.mTranslate = mTranslate;
exports.mScale = mScale;
exports.mShear = mShear;
exports.vRotate = vRotate;
exports.vRotatePointAround = vRotatePointAround;
exports.vMidpoint = vMidpoint;
exports.vAngle = vAngle;
exports.vAlongAngle = vAlongAngle;
exports.vFastDist = vFastDist;
exports.vDist = vDist;
exports.vDot = vDot;
exports.vPerpDot = vPerpDot;
exports.vTriangleArea = vTriangleArea;
exports.vColinear = vColinear;
exports.vDet = vDet;
/* end exports */
},{}],2:[function(require,module,exports){
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

var vecModule = require("../js/vec.module.js");

// vec-la-fp places nice names in a member named "vec"
var vec = vecModule.vec;

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

module.exports = {
  TriadBalance2to3: TriadBalance2to3,
  invertTriadBalance2to3: invertTriadBalance2to3,
  eqPointOnEdgeAlgebraically: eqPointOnEdgeAlgebraically,
  eqEdgeAlgebraically: eqEdgeAlgebraically,
  GetRayToLineSegmentIntersection : GetRayToLineSegmentIntersection,
  L1LENGTH: L1LENGTH,
  L2LENGTH: L2LENGTH,
  L1: L1,
  L2: L2};

},{"../js/vec.module.js":1}],3:[function(require,module,exports){
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

var m = require("../src/TriadBalanceMath.js");
var vecModule = require("../js/vec.module.js");

// vec-la-fp places nice names in a member named "vec"
var vec = vecModule.vec;

function computeFinverseF(wtc,norm,p) {
   return m.invertTriadBalance2to3(
     m.TriadBalance2to3(p,wtc,norm),
     wtc,
     norm);
}

function test_eqPointOnAlgebraically(wtc) {
  var rof = m.eqPointOnEdgeAlgebraically(wtc,[0,0]);
  console.assert(rof == null);

  const halfvert = vec.scale(1/2,wtc[0]);
  var r0 = m.eqPointOnEdgeAlgebraically(wtc,halfvert);
  console.assert(vec.near(1e-5,r0,wtc[0]));

  var r1 = m.eqPointOnEdgeAlgebraically(wtc,wtc[1]);
  console.assert(vec.near(1e-5,r1,wtc[1]));

  var r2 = m.eqPointOnEdgeAlgebraically(wtc,wtc[2]);
  console.assert(vec.near(1e-5,r2,wtc[2]));

  // slope of 1
  var rone = m.eqPointOnEdgeAlgebraically(wtc,[1,1]);
  console.assert(vec.colinear(wtc[1],wtc[2],rone));

  // slope of -1
  var rneg_one = m.eqPointOnEdgeAlgebraically(wtc,[-1,1]);
  console.assert(vec.colinear(wtc[0],wtc[2],rneg_one));

  // almost straight down
  var rdown = m.eqPointOnEdgeAlgebraically(wtc,[0.01,-1]);
  console.assert(vec.colinear(rdown,wtc[0],wtc[1]));

  return true;
}

function test_eqEdgeAlebraically(wtc) {
  var rof = m.eqEdgeAlgebraically(wtc,[0,0]);
  console.assert(rof.length == 0);

  var r0 = m.eqEdgeAlgebraically(wtc,wtc[0]);
  console.assert(r0.length == 2);

  var r1 = m.eqEdgeAlgebraically(wtc,wtc[1]);
  console.assert(r1.length == 2);

  var r2 = m.eqEdgeAlgebraically(wtc,wtc[2]);
  console.assert(r2.length == 2);

  // slope of 1
  var rone = m.eqEdgeAlgebraically(wtc,[1,1]);
  console.assert(rone.length == 1);
  console.assert(rone[0] == 1);

  var rneg_one = m.eqEdgeAlgebraically(wtc,[-1,1]);
  console.assert(rneg_one.length == 1);
  console.assert(rneg_one[0] == 2);

  var rdown = m.eqEdgeAlgebraically(wtc,[0.01,-1]);
  console.assert(rdown.length == 1);
  console.assert(rdown[0] == 0);
  return true;
}


function testGetRayToLineSegmentIntersection(wtc) {
  let ro = [0,0];
  let rd = [1,1];
  let p1 = [0,10];
  let p2 = [10,0];
  var r = m.GetRayToLineSegmentIntersection(ro,rd,p1,p2)[0];
  console.assert(r[0] == r[1]);
  var r = m.GetRayToLineSegmentIntersection(ro,rd,p2,p1)[0];
  console.assert(r[0] == r[1]);

  var rd1 = [94.1015625,-36.36328125];
  let c0 = wtc[0];
  let c1 = wtc[1];
  let c2 = wtc[2];

  var r12 = m.GetRayToLineSegmentIntersection(ro,rd1,c1,c2);
  console.assert(r12 != null);

  {
    var down = m.GetRayToLineSegmentIntersection(ro,[0,-4999],c0,c1);
    console.assert(vec.scalarNear(1e-5,down[0][1],wtc[0][1]));
  }
  return true;
}


function testTriadBalance2to3(wtc) {
  let pv = m.TriadBalance2to3([30000,30],wtc,m.L1);
  console.assert(vec.scalarNear(1e-5,m.L1LENGTH(pv),1));

  let pyv = m.TriadBalance2to3([0,wtc[2][1]],wtc,m.L1);
  console.assert(vec.scalarNear(1e-5,m.L1LENGTH(pyv),1));
  return true;
}

function testOriginAndVertices(wtc) {
  // The origin should return a perfectly balanced vector
  let ov = m.TriadBalance2to3([0,0],wtc,m.L1);
  console.assert(vec.scalarNear(1e-5,ov[0],1/3));
  console.assert(vec.scalarNear(1e-5,ov[1],1/3));
  console.assert(vec.scalarNear(1e-5,ov[2],1/3));

  // A vertex should return a vector with a 1 in exactly 1 position
  console.assert(vec.scalarNear(1e-5,m.TriadBalance2to3(wtc[0],wtc,m.L1)[0],1));
  console.assert(vec.scalarNear(1e-5,m.TriadBalance2to3(wtc[1],wtc,m.L1)[1],1));
  console.assert(vec.scalarNear(1e-5,m.TriadBalance2to3(wtc[2],wtc,m.L1)[2],1));

  return true;
}

// It is extremely valuable to be able to compute the
// inverse of the TriadBalance algorithm for testing,
// although this is not quite a true invesion because
// points outside the triangle are brought to the triangle
// edge.
function computeFinverseF(wtc,norm,p) {
   return m.invertTriadBalance2to3(
     m.TriadBalance2to3(p,wtc,norm),
     wtc,
     norm);
}

function testInversion(wtc) {
  // This insures we are within the triangle
  let d = vec.dist(wtc[0],wtc[1]);
  let p = [d/8,d/8];

  let vpc = vec.sub(
    computeFinverseF(wtc,m.L1,p),
    p);
  console.assert(vec.scalarNear(1e-5,vec.mag(vpc),0));

  let py = [0,wtc[2][1]];
  let vpyc = vec.sub(
    computeFinverseF(wtc,m.L1,py),
    py);
  console.assert(vec.scalarNear(1e-5,vec.mag(vpyc),0));
  return true;
}

function testInversionOutside(wtc) {
  // This insures we are outside the triangle
  let d = vec.dist(wtc[0],wtc[1]);

  var vp_inv = computeFinverseF(wtc,m.L1,[d,d]);
  // we now want to test that the invesion has a slope of 1
  console.assert(vec.scalarNear(1e-5,vp_inv[1]/vp_inv[0],1));
  console.assert(vp_inv[1] > 0);
  return true;
}

function testInversionNegativeY(wtc) {
  // This insures we are within the triangle
  let d = vec.dist(wtc[0],wtc[1]);
  let p = [0,-d/4];
  var vp_inv = computeFinverseF(wtc,m.L1,p);
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
  allTrue &= testInversionWithACircle(wtc,m.L1,small_radius);
  allTrue &= testInversionWithACircle(wtc,m.L2,small_radius);
  let large_radius = vec.mag(wtcp)*3;
  allTrue &= testInversionWithACircle(wtc,m.L1,large_radius);
  allTrue &= testInversionWithACircle(wtc,m.L2,large_radius);
  return allTrue;
}

module.exports = {
  testAllTriadBalance: testAllTriadBalance
};

},{"../js/vec.module.js":1,"../src/TriadBalanceMath.js":2}]},{},[3])(3)
});

//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyaWZ5L25vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJqcy92ZWMubW9kdWxlLmpzIiwic3JjL1RyaWFkQmFsYW5jZU1hdGguanMiLCJ0ZXN0L1RyaWFkQmFsYW5jZVRlc3QuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDdlRBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDL1dBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0EiLCJmaWxlIjoiZ2VuZXJhdGVkLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXNDb250ZW50IjpbIihmdW5jdGlvbigpe2Z1bmN0aW9uIHIoZSxuLHQpe2Z1bmN0aW9uIG8oaSxmKXtpZighbltpXSl7aWYoIWVbaV0pe3ZhciBjPVwiZnVuY3Rpb25cIj09dHlwZW9mIHJlcXVpcmUmJnJlcXVpcmU7aWYoIWYmJmMpcmV0dXJuIGMoaSwhMCk7aWYodSlyZXR1cm4gdShpLCEwKTt2YXIgYT1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK2krXCInXCIpO3Rocm93IGEuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixhfXZhciBwPW5baV09e2V4cG9ydHM6e319O2VbaV1bMF0uY2FsbChwLmV4cG9ydHMsZnVuY3Rpb24ocil7dmFyIG49ZVtpXVsxXVtyXTtyZXR1cm4gbyhufHxyKX0scCxwLmV4cG9ydHMscixlLG4sdCl9cmV0dXJuIG5baV0uZXhwb3J0c31mb3IodmFyIHU9XCJmdW5jdGlvblwiPT10eXBlb2YgcmVxdWlyZSYmcmVxdWlyZSxpPTA7aTx0Lmxlbmd0aDtpKyspbyh0W2ldKTtyZXR1cm4gb31yZXR1cm4gcn0pKCkiLCJcInVzZSBzdHJpY3RcIjtcblxuT2JqZWN0LmRlZmluZVByb3BlcnR5KGV4cG9ydHMsIFwiX19lc01vZHVsZVwiLCB7XG4gIHZhbHVlOiB0cnVlXG59KTtcblxuZnVuY3Rpb24gX3RvQ29uc3VtYWJsZUFycmF5KGFycikgeyBpZiAoQXJyYXkuaXNBcnJheShhcnIpKSB7IGZvciAodmFyIGkgPSAwLCBhcnIyID0gQXJyYXkoYXJyLmxlbmd0aCk7IGkgPCBhcnIubGVuZ3RoOyBpKyspIHsgYXJyMltpXSA9IGFycltpXTsgfSByZXR1cm4gYXJyMjsgfSBlbHNlIHsgcmV0dXJuIEFycmF5LmZyb20oYXJyKTsgfSB9XG5cbi8vIGN1cnJ5IDo6IChhIC0+IGIgLT4gLi4uIC0+IG4pIC0+IChhIC0+IGIpIC0+IChiIC0+IC4uLikgLT4gKC4uLiAtPiBuKVxudmFyIGN1cnJ5ID0gZnVuY3Rpb24gY3VycnkoZm4pIHtcbiAgdmFyIGN1cnJpZWQgPSBmdW5jdGlvbiBjdXJyaWVkKCkge1xuICAgIGZvciAodmFyIF9sZW4gPSBhcmd1bWVudHMubGVuZ3RoLCBhcmdzID0gQXJyYXkoX2xlbiksIF9rZXkgPSAwOyBfa2V5IDwgX2xlbjsgX2tleSsrKSB7XG4gICAgICBhcmdzW19rZXldID0gYXJndW1lbnRzW19rZXldO1xuICAgIH1cblxuICAgIGlmIChhcmdzLmxlbmd0aCA+PSBmbi5sZW5ndGgpIHtcbiAgICAgIHJldHVybiBmbi5hcHBseSh1bmRlZmluZWQsIGFyZ3MpO1xuICAgIH1cbiAgICByZXR1cm4gZnVuY3Rpb24gKCkge1xuICAgICAgZm9yICh2YXIgX2xlbjIgPSBhcmd1bWVudHMubGVuZ3RoLCBhcmdzTmV4dCA9IEFycmF5KF9sZW4yKSwgX2tleTIgPSAwOyBfa2V5MiA8IF9sZW4yOyBfa2V5MisrKSB7XG4gICAgICAgIGFyZ3NOZXh0W19rZXkyXSA9IGFyZ3VtZW50c1tfa2V5Ml07XG4gICAgICB9XG5cbiAgICAgIHJldHVybiBjdXJyaWVkLmFwcGx5KHVuZGVmaW5lZCwgYXJncy5jb25jYXQoYXJnc05leHQpKTtcbiAgICB9O1xuICB9O1xuICByZXR1cm4gY3VycmllZDtcbn07XG5cbi8vIHBpcGUgOjogKGEgLT4gYikgLT4gKGIgLT4gLi4uKSAtPiAoLi4uIC0+IG4pXG52YXIgcGlwZSA9IGZ1bmN0aW9uIHBpcGUoZm4xKSB7XG4gIGZvciAodmFyIF9sZW4zID0gYXJndW1lbnRzLmxlbmd0aCwgZnVuY3Rpb25zID0gQXJyYXkoX2xlbjMgPiAxID8gX2xlbjMgLSAxIDogMCksIF9rZXkzID0gMTsgX2tleTMgPCBfbGVuMzsgX2tleTMrKykge1xuICAgIGZ1bmN0aW9uc1tfa2V5MyAtIDFdID0gYXJndW1lbnRzW19rZXkzXTtcbiAgfVxuXG4gIHJldHVybiBmdW5jdGlvbiAoKSB7XG4gICAgcmV0dXJuIGZ1bmN0aW9ucy5yZWR1Y2UoZnVuY3Rpb24gKGFjYywgZm4pIHtcbiAgICAgIHJldHVybiBmbihhY2MpO1xuICAgIH0sIGZuMS5hcHBseSh1bmRlZmluZWQsIGFyZ3VtZW50cykpO1xuICB9O1xufTtcblxuLy8gY29tcG9zZSA6OiAoLi4uIC0+IG4pIC0+IChiIC0+IC4uLikgLT4gKGEgLT4gYilcbnZhciBjb21wb3NlID0gZnVuY3Rpb24gY29tcG9zZSgpIHtcbiAgZm9yICh2YXIgX2xlbjQgPSBhcmd1bWVudHMubGVuZ3RoLCBmdW5jdGlvbnMgPSBBcnJheShfbGVuNCksIF9rZXk0ID0gMDsgX2tleTQgPCBfbGVuNDsgX2tleTQrKykge1xuICAgIGZ1bmN0aW9uc1tfa2V5NF0gPSBhcmd1bWVudHNbX2tleTRdO1xuICB9XG5cbiAgcmV0dXJuIHBpcGUuYXBwbHkodW5kZWZpbmVkLCBfdG9Db25zdW1hYmxlQXJyYXkoZnVuY3Rpb25zLnJldmVyc2UoKSkpO1xufTtcblxuLy8gdkFkZCA6OiBWZWN0b3IgLT4gVmVjdG9yIC0+IFZlY3RvclxudmFyIHZBZGQgPSBjdXJyeShmdW5jdGlvbiAodiwgdjIpIHtcbiAgcmV0dXJuIFt2WzBdICsgdjJbMF0sIHZbMV0gKyB2MlsxXV07XG59KTtcblxuLy8gdkFkZDMgOjogVmVjdG9yIC0+IFZlY3RvciAtPiBWZWN0b3IgLT4gVmVjdG9yXG52YXIgdkFkZDMgPSBjdXJyeShmdW5jdGlvbiAodiwgdjIsIHYzKSB7XG4gIHJldHVybiBbdlswXSArIHYyWzBdICsgdjNbMF0sIHZbMV0gKyB2MlsxXSArIHYzWzFdXTtcbn0pO1xuXG4vLyB2QWRkQWxsIDo6IFtWZWN0b3JdIC0+IFZlY3RvclxudmFyIHZBZGRBbGwgPSBmdW5jdGlvbiB2QWRkQWxsKHZzKSB7XG4gIHJldHVybiB2cy5yZWR1Y2UodkFkZCwgWzAsIDBdKTtcbn07XG5cbi8vIHZTdWIgOjogVmVjdG9yIC0+IFZlY3RvciAtPiBWZWN0b3JcbnZhciB2U3ViID0gY3VycnkoZnVuY3Rpb24gKHYsIHYyKSB7XG4gIHJldHVybiBbdlswXSAtIHYyWzBdLCB2WzFdIC0gdjJbMV1dO1xufSk7XG5cbi8vIHZTdWIzIDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gVmVjdG9yIC0+IFZlY3RvclxudmFyIHZTdWIzID0gY3VycnkoZnVuY3Rpb24gKHYsIHYyLCB2Mykge1xuICByZXR1cm4gW3ZbMF0gLSB2MlswXSAtIHYzWzBdLCB2WzFdIC0gdjJbMV0gLSB2M1sxXV07XG59KTtcblxuLy8gdlN1YkFsbCA6OiBbVmVjdG9yXSAtPiBWZWN0b3JcbnZhciB2U3ViQWxsID0gZnVuY3Rpb24gdlN1YkFsbCh2cykge1xuICByZXR1cm4gdnMuc2xpY2UoMSkucmVkdWNlKHZTdWIsIHZzLnNsaWNlKDAsIDEpWzBdKTtcbn07XG5cbi8vIHZNYWcgOjogVmVjdG9yIC0+IE51bWJlclxudmFyIHZNYWcgPSBmdW5jdGlvbiB2TWFnKHYpIHtcbiAgcmV0dXJuIE1hdGguc3FydCh2WzBdICogdlswXSArIHZbMV0gKiB2WzFdKTtcbn07XG5cbi8vIHZOb3JtYWwgOjogVmVjdG9yIC0+IFZlY3RvclxudmFyIHZOb3JtYWwgPSBmdW5jdGlvbiB2Tm9ybWFsKHYpIHtcbiAgcmV0dXJuIFstdlsxXSwgdlswXV07XG59O1xuXG4vLyB2U2NhbGUgOjogTnVtYmVyIC0+IFZlY3RvclxudmFyIHZTY2FsZSA9IGN1cnJ5KGZ1bmN0aW9uIChzYywgdikge1xuICByZXR1cm4gW3ZbMF0gKiBzYywgdlsxXSAqIHNjXTtcbn0pO1xuXG4vLyB2VG93YXJkcyA6OiBOdW1iZXIgLT4gVmVjdG9yIC0+IFZlY3RvciAtPiBWZWN0b3JcbnZhciB2VG93YXJkcyA9IGN1cnJ5KGZ1bmN0aW9uICh0LCB2MSwgdjIpIHtcbiAgdmFyIGQgPSB2U3ViKHYyLCB2MSk7XG4gIHZhciBzYyA9IHZNYWcoZCkgKiB0O1xuICByZXR1cm4gdkFkZCh2MSwgdlNjYWxlKHNjLCB2Tm9ybShkKSkpO1xufSk7XG5cbi8vIHZMZXJwIDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gTnVtYmVyIC0+IFZlY3RvclxudmFyIHZMZXJwID0gY3VycnkoZnVuY3Rpb24gKHYxLCB2MiwgdCkge1xuICByZXR1cm4gdlRvd2FyZHModCwgdjEsIHYyKTtcbn0pO1xuXG4vLyB2U2NhbGFyTmVhciA6OiBOdW1iZXIgLT4gTnVtYmVyIC0+IE51bWJlciAtPiBib29sXG52YXIgdlNjYWxhck5lYXIgPSBjdXJyeShmdW5jdGlvbiAoZSwgYSwgYikge1xuICByZXR1cm4gTWF0aC5hYnMoYSAtIGIpIDwgZTtcbn0pO1xuXG4vLyB2TmVhciA6OiBOdW1iZXIgLT4gVmVjdG9yIC0+IFZlY3RvciAtPiBib29sXG52YXIgdk5lYXIgPSBjdXJyeShmdW5jdGlvbiAoZSwgYSwgYikge1xuICByZXR1cm4gdlNjYWxhck5lYXIoZSwgYVswXSwgYlswXSkgJiYgdlNjYWxhck5lYXIoZSwgYVsxXSwgYlsxXSk7XG59KTtcblxuLy8gdkNsYW1wTWFnIDo6IE51bWJlciAtPiBOdW1iZXIgLT4gVmVjdG9yIC0+IFZlY3RvclxudmFyIHZDbGFtcE1hZyA9IGN1cnJ5KGZ1bmN0aW9uIChtaW4sIG1heCwgdikge1xuICB2YXIgZCA9IHZlYy5tYWcodik7XG4gIGlmIChkIDwgbWluKSByZXR1cm4gdmVjLnNjYWxlKG1pbiAvIGQsIHYpO2Vsc2UgaWYgKGQgPiBtYXgpIHJldHVybiB2ZWMuc2NhbGUobWF4IC8gZCwgdik7ZWxzZSByZXR1cm4gdjtcbn0pO1xuXG4vLyB2Tm9ybSA6OiBWZWN0b3IgLT4gVmVjdG9yXG52YXIgdk5vcm0gPSBmdW5jdGlvbiB2Tm9ybSh2KSB7XG4gIHZhciBtYWcgPSB2TWFnKHYpO1xuICByZXR1cm4gW3ZbMF0gLyBtYWcsIHZbMV0gLyBtYWddO1xufTtcblxuLy8gbUlkIDo6IE1hdHJpeFxudmFyIG1JZCA9IE9iamVjdC5mcmVlemUoWzEsIDAsIDAsIDAsIDEsIDAsIDAsIDAsIDFdKTtcblxuLy8gdkNyZWF0ZU1hdHJpeCA6OiBOdW1iZXIgLT4gTnVtYmVyIC0+IE51bWJlciAtPiBOdW1iZXIgLT4gTnVtYmVyIC0+IE51bWJlciAtPiBNYXRyaXhcbnZhciB2Q3JlYXRlTWF0cml4ID0gZnVuY3Rpb24gdkNyZWF0ZU1hdHJpeCgpIHtcbiAgdmFyIGEgPSBhcmd1bWVudHMubGVuZ3RoID4gMCAmJiBhcmd1bWVudHNbMF0gIT09IHVuZGVmaW5lZCA/IGFyZ3VtZW50c1swXSA6IDE7XG4gIHZhciBiID0gYXJndW1lbnRzLmxlbmd0aCA+IDEgJiYgYXJndW1lbnRzWzFdICE9PSB1bmRlZmluZWQgPyBhcmd1bWVudHNbMV0gOiAwO1xuICB2YXIgYyA9IGFyZ3VtZW50cy5sZW5ndGggPiAyICYmIGFyZ3VtZW50c1syXSAhPT0gdW5kZWZpbmVkID8gYXJndW1lbnRzWzJdIDogMDtcbiAgdmFyIGQgPSBhcmd1bWVudHMubGVuZ3RoID4gMyAmJiBhcmd1bWVudHNbM10gIT09IHVuZGVmaW5lZCA/IGFyZ3VtZW50c1szXSA6IDE7XG4gIHZhciB0eCA9IGFyZ3VtZW50cy5sZW5ndGggPiA0ICYmIGFyZ3VtZW50c1s0XSAhPT0gdW5kZWZpbmVkID8gYXJndW1lbnRzWzRdIDogMDtcbiAgdmFyIHR5ID0gYXJndW1lbnRzLmxlbmd0aCA+IDUgJiYgYXJndW1lbnRzWzVdICE9PSB1bmRlZmluZWQgPyBhcmd1bWVudHNbNV0gOiAwO1xuICByZXR1cm4gW2EsIGMsIHR4LCBiLCBkLCB0eSwgMCwgMCwgMV07XG59O1xuXG4vLyB2VHJhbnNmb3JtIDo6IE1hdHJpeCAtPiBWZWN0b3IgLT4gVmVjdG9yXG52YXIgdlRyYW5zZm9ybSA9IGN1cnJ5KGZ1bmN0aW9uIChtLCB2KSB7XG4gIHJldHVybiBbdlswXSAqIG1bMF0gKyB2WzFdICogbVsxXSArIG1bMl0sIHZbMF0gKiBtWzNdICsgdlsxXSAqIG1bNF0gKyBtWzVdXTtcbn0pO1xuXG4vLyBtQ29tcG9zZSA6OiBNYXRyaXggLT4gTWF0cml4IC0+IE1hdHJpeFxudmFyIG1Db21wb3NlID0gY3VycnkoZnVuY3Rpb24gKG0sIG0yKSB7XG4gIHJldHVybiBbbVswXSAqIG0yWzBdICsgbVsxXSAqIG0yWzNdICsgbVsyXSAqIG0yWzZdLCBtWzBdICogbTJbMV0gKyBtWzFdICogbTJbNF0gKyBtWzJdICogbTJbN10sIG1bMF0gKiBtMlsyXSArIG1bMV0gKiBtMls1XSArIG1bMl0gKiBtMls4XSwgbVszXSAqIG0yWzBdICsgbVs0XSAqIG0yWzNdICsgbVs1XSAqIG0yWzZdLCBtWzNdICogbTJbMV0gKyBtWzRdICogbTJbNF0gKyBtWzVdICogbTJbN10sIG1bM10gKiBtMlsyXSArIG1bNF0gKiBtMls1XSArIG1bNV0gKiBtMls4XSwgbVs2XSAqIG0yWzBdICsgbVs3XSAqIG0yWzNdICsgbVs4XSAqIG0yWzZdLCBtWzZdICogbTJbMV0gKyBtWzddICogbTJbNF0gKyBtWzhdICogbTJbN10sIG1bNl0gKiBtMlsyXSArIG1bN10gKiBtMls1XSArIG1bOF0gKiBtMls4XV07XG59KTtcblxuLy8gbVJvdGF0ZSA6OiBOdW1iZXIgLT4gTWF0cml4IC0+IE1hdHJpeFxudmFyIG1Sb3RhdGUgPSBmdW5jdGlvbiBtUm90YXRlKGEpIHtcbiAgcmV0dXJuIG1Db21wb3NlKFtNYXRoLmNvcyhhKSwgLU1hdGguc2luKGEpLCAwLCBNYXRoLnNpbihhKSwgTWF0aC5jb3MoYSksIDAsIDAsIDAsIDFdKTtcbn07XG5cbi8vIG1UcmFuc2xhdGUgOjogVmVjdG9yIC0+IE1hdHJpeCAtPiBNYXRyaXhcbnZhciBtVHJhbnNsYXRlID0gZnVuY3Rpb24gbVRyYW5zbGF0ZSh2KSB7XG4gIHJldHVybiBtQ29tcG9zZShbMSwgMCwgdlswXSwgMCwgMSwgdlsxXSwgMCwgMCwgMV0pO1xufTtcblxuLy8gbVNjYWxlIDo6IFZlY3RvciAtPiBNYXRyaXggLT4gTWF0cml4XG52YXIgbVNjYWxlID0gZnVuY3Rpb24gbVNjYWxlKHYpIHtcbiAgcmV0dXJuIG1Db21wb3NlKFt2WzBdLCAwLCAwLCAwLCB2WzFdLCAwLCAwLCAwLCAxXSk7XG59O1xuXG4vLyBtU2hlYXIgOjogVmVjdG9yIC0+IE1hdHJpeCAtPiBNYXRyaXhcbnZhciBtU2hlYXIgPSBmdW5jdGlvbiBtU2hlYXIodikge1xuICByZXR1cm4gbUNvbXBvc2UoWzEsIHZbMF0sIDAsIHZbMV0sIDEsIDAsIDAsIDAsIDFdKTtcbn07XG5cbi8vIHZSb3RhdGUgOjogTnVtYmVyIC0+IFZlY3RvciAtPiBWZWN0b3JcbnZhciB2Um90YXRlID0gY3VycnkoZnVuY3Rpb24gKGEsIHYpIHtcbiAgcmV0dXJuIFt2WzBdICogTWF0aC5jb3MoYSkgLSB2WzFdICogTWF0aC5zaW4oYSksIHZbMF0gKiBNYXRoLnNpbihhKSArIHZbMV0gKiBNYXRoLmNvcyhhKV07XG59KTtcblxuLy8gdlJvdGF0ZVBvaW50QXJvdW5kIDo6IE51bWJlciAtPiBWZWN0b3IgLT4gVmVjdG9yIC0+IFZlY3RvclxudmFyIHZSb3RhdGVQb2ludEFyb3VuZCA9IGN1cnJ5KGZ1bmN0aW9uIChhLCBjcCwgdikge1xuICB2YXIgdjIgPSB2U3ViKHYsIGNwKTtcbiAgcmV0dXJuIHZBZGQoY3AsIFt2MlswXSAqIE1hdGguY29zKGEpIC0gdjJbMV0gKiBNYXRoLnNpbihhKSwgdjJbMF0gKiBNYXRoLnNpbihhKSArIHYyWzFdICogTWF0aC5jb3MoYSldKTtcbn0pO1xuXG4vLyB2TWlkcG9pbnQgOjogVmVjdG9yIC0+IFZlY3RvciAtPiBWZWN0b3JcbnZhciB2TWlkcG9pbnQgPSBjdXJyeShmdW5jdGlvbiAodiwgdjIpIHtcbiAgcmV0dXJuIHZTY2FsZSgwLjUsIHZBZGQodiwgdjIpKTtcbn0pO1xuXG4vLyB2QW5nbGUgOjogTnVtYmVyIC0+IFZlY3RvclxudmFyIHZBbmdsZSA9IGZ1bmN0aW9uIHZBbmdsZShhKSB7XG4gIHJldHVybiBbTWF0aC5jb3MoYSksIE1hdGguc2luKGEpXTtcbn07XG5cbi8vIHZBbG9uZ0FuZ2xlIDo6IE51bWJlciAtPiBOdW1iZXIgLT4gVmVjdG9yXG52YXIgdkFsb25nQW5nbGUgPSBjdXJyeShmdW5jdGlvbiAoYSwgciwgdikge1xuICByZXR1cm4gY29tcG9zZSh2QWRkKHYpLCB2U2NhbGUociksIHZBbmdsZSkoYSk7XG59KTtcblxuLy8gdkZhc3REaXN0IDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gTnVtYmVyXG52YXIgdkZhc3REaXN0ID0gY3VycnkoZnVuY3Rpb24gKHYsIHYyKSB7XG4gIHJldHVybiBNYXRoLnBvdyh2MlswXSAtIHZbMF0sIDIpICsgTWF0aC5wb3codjJbMV0gLSB2WzFdLCAyKTtcbn0pO1xuXG4vLyB2RGlzdCA6OiBWZWN0b3IgLT4gVmVjdG9yIC0+IE51bWJlclxudmFyIHZEaXN0ID0gY3VycnkoZnVuY3Rpb24gKHYsIHYyKSB7XG4gIHJldHVybiBNYXRoLmh5cG90KHYyWzBdIC0gdlswXSwgdjJbMV0gLSB2WzFdKTtcbn0pO1xuXG4vLyB2RG90IDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gTnVtYmVyXG52YXIgdkRvdCA9IGN1cnJ5KGZ1bmN0aW9uICh2LCB2Mikge1xuICByZXR1cm4gdlswXSAqIHYyWzBdICsgdlsxXSAqIHYyWzFdO1xufSk7XG5cbi8vIHZQZXJwRG90IDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gTnVtYmVyXG52YXIgdlBlcnBEb3QgPSBjdXJyeShmdW5jdGlvbiAodiwgdjIpIHtcbiAgcmV0dXJuIHZbMF0gKiB2MlsxXSAtIHZbMV0gKiB2MlswXTtcbn0pO1xuXG4vLyB2VHJpYW5nbGVBcmVhIDo6IFZlY3RvciAtPiBWZWN0b3IgLT4gVmVjdG9yIC0+IE51bWJlclxudmFyIHZUcmlhbmdsZUFyZWEgPSBjdXJyeShmdW5jdGlvbiAoYSwgYiwgYykge1xuICByZXR1cm4gKChiWzBdIC0gYVswXSkgKiAoY1sxXSAtIGFbMV0pIC0gKGNbMF0gLSBhWzBdKSAqIChiWzFdIC0gYVsxXSkpIC8gMjtcbn0pO1xuXG4vLyB2Q29saW5lYXIgOjogVmVjdG9yIC0+IFZlY3RvciAtPiBWZWN0b3IgLT4gYm9vbFxudmFyIHZDb2xpbmVhciA9IGN1cnJ5KGZ1bmN0aW9uICh2MCwgdjEsIHYyKSB7XG4gIHJldHVybiB2U2NhbGFyTmVhcigxZS00LCB2VHJpYW5nbGVBcmVhKHYwLCB2MSwgdjIpLCAwKTtcbn0pO1xuXG4vLyB2RGV0IDo6IE1hdHJpeCAtPiBOdW1iZXJcbnZhciB2RGV0ID0gZnVuY3Rpb24gdkRldChtKSB7XG4gIHJldHVybiBtWzBdICogbVs0XSAtIG1bM10gKiBtWzFdO1xufTtcblxudmFyIHZlYyA9IHtcbiAgYWRkOiB2QWRkLFxuICBhZGQzOiB2QWRkMyxcbiAgYWRkQWxsOiB2QWRkQWxsLFxuICBzdWI6IHZTdWIsXG4gIHN1YjM6IHZTdWIzLFxuICBzdWJBbGw6IHZTdWJBbGwsXG4gIG1hZzogdk1hZyxcbiAgbm9ybWFsOiB2Tm9ybWFsLFxuICBzY2FsZTogdlNjYWxlLFxuICB0b3dhcmRzOiB2VG93YXJkcyxcbiAgbGVycDogdkxlcnAsXG4gIHNjYWxhck5lYXI6IHZTY2FsYXJOZWFyLFxuICBuZWFyOiB2TmVhcixcbiAgY2xhbXBNYWc6IHZDbGFtcE1hZyxcbiAgbm9ybTogdk5vcm0sXG4gIG1JZDogbUlkLFxuICBjcmVhdGVNYXRyaXg6IHZDcmVhdGVNYXRyaXgsXG4gIHRyYW5zZm9ybTogdlRyYW5zZm9ybSxcbiAgbUNvbXBvc2U6IG1Db21wb3NlLFxuICBtUm90YXRlOiBtUm90YXRlLFxuICBtVHJhbnNsYXRlOiBtVHJhbnNsYXRlLFxuICBtU2NhbGU6IG1TY2FsZSxcbiAgbVNoZWFyOiBtU2hlYXIsXG4gIHJvdGF0ZTogdlJvdGF0ZSxcbiAgcm90YXRlUG9pbnRBcm91bmQ6IHZSb3RhdGVQb2ludEFyb3VuZCxcbiAgbWlkcG9pbnQ6IHZNaWRwb2ludCxcbiAgYW5nbGU6IHZBbmdsZSxcbiAgYWxvbmdBbmdsZTogdkFsb25nQW5nbGUsXG4gIGZhc3REaXN0OiB2RmFzdERpc3QsXG4gIGRpc3Q6IHZEaXN0LFxuICBkb3Q6IHZEb3QsXG4gIHBlcnBkb3Q6IHZQZXJwRG90LFxuICB0cmlhbmdsZUFyZWE6IHZUcmlhbmdsZUFyZWEsXG4gIGNvbGluZWFyOiB2Q29saW5lYXIsXG4gIGRldDogdkRldFxufTtcblxuLyogc3RhcnQgZXhwb3J0cyAqL1xuZXhwb3J0cy5kZWZhdWx0ID0gdmVjO1xuZXhwb3J0cy52ZWMgPSB2ZWM7XG5leHBvcnRzLnZBZGQgPSB2QWRkO1xuZXhwb3J0cy52QWRkMyA9IHZBZGQzO1xuZXhwb3J0cy52QWRkQWxsID0gdkFkZEFsbDtcbmV4cG9ydHMudlN1YiA9IHZTdWI7XG5leHBvcnRzLnZTdWIzID0gdlN1YjM7XG5leHBvcnRzLnZTdWJBbGwgPSB2U3ViQWxsO1xuZXhwb3J0cy52TWFnID0gdk1hZztcbmV4cG9ydHMudk5vcm1hbCA9IHZOb3JtYWw7XG5leHBvcnRzLnZTY2FsZSA9IHZTY2FsZTtcbmV4cG9ydHMudlRvd2FyZHMgPSB2VG93YXJkcztcbmV4cG9ydHMudkxlcnAgPSB2TGVycDtcbmV4cG9ydHMudlNjYWxhck5lYXIgPSB2U2NhbGFyTmVhcjtcbmV4cG9ydHMudk5lYXIgPSB2TmVhcjtcbmV4cG9ydHMudkNsYW1wTWFnID0gdkNsYW1wTWFnO1xuZXhwb3J0cy52Tm9ybSA9IHZOb3JtO1xuZXhwb3J0cy5tSWQgPSBtSWQ7XG5leHBvcnRzLnZDcmVhdGVNYXRyaXggPSB2Q3JlYXRlTWF0cml4O1xuZXhwb3J0cy52VHJhbnNmb3JtID0gdlRyYW5zZm9ybTtcbmV4cG9ydHMubUNvbXBvc2UgPSBtQ29tcG9zZTtcbmV4cG9ydHMubVJvdGF0ZSA9IG1Sb3RhdGU7XG5leHBvcnRzLm1UcmFuc2xhdGUgPSBtVHJhbnNsYXRlO1xuZXhwb3J0cy5tU2NhbGUgPSBtU2NhbGU7XG5leHBvcnRzLm1TaGVhciA9IG1TaGVhcjtcbmV4cG9ydHMudlJvdGF0ZSA9IHZSb3RhdGU7XG5leHBvcnRzLnZSb3RhdGVQb2ludEFyb3VuZCA9IHZSb3RhdGVQb2ludEFyb3VuZDtcbmV4cG9ydHMudk1pZHBvaW50ID0gdk1pZHBvaW50O1xuZXhwb3J0cy52QW5nbGUgPSB2QW5nbGU7XG5leHBvcnRzLnZBbG9uZ0FuZ2xlID0gdkFsb25nQW5nbGU7XG5leHBvcnRzLnZGYXN0RGlzdCA9IHZGYXN0RGlzdDtcbmV4cG9ydHMudkRpc3QgPSB2RGlzdDtcbmV4cG9ydHMudkRvdCA9IHZEb3Q7XG5leHBvcnRzLnZQZXJwRG90ID0gdlBlcnBEb3Q7XG5leHBvcnRzLnZUcmlhbmdsZUFyZWEgPSB2VHJpYW5nbGVBcmVhO1xuZXhwb3J0cy52Q29saW5lYXIgPSB2Q29saW5lYXI7XG5leHBvcnRzLnZEZXQgPSB2RGV0O1xuLyogZW5kIGV4cG9ydHMgKi8iLCIvLyBDb3B5cmlnaHQgMjAxOSwgUm9iZXJ0IEwuIFJlYWRcbi8vIFRoaXMgZmlsZSBpcyBwYXJ0IG9mIFRyaWFkQmFsYW5jZS5cbi8vXG4vLyBUcmlhZEJhbGFuY2UgaXMgZnJlZSBzb2Z0d2FyZTogeW91IGNhbiByZWRpc3RyaWJ1dGUgaXQgYW5kL29yIG1vZGlmeVxuLy8gaXQgdW5kZXIgdGhlIHRlcm1zIG9mIHRoZSBHTlUgR2VuZXJhbCBQdWJsaWMgTGljZW5zZSBhcyBwdWJsaXNoZWQgYnlcbi8vIHRoZSBGcmVlIFNvZnR3YXJlIEZvdW5kYXRpb24sIGVpdGhlciB2ZXJzaW9uIDMgb2YgdGhlIExpY2Vuc2UsIG9yXG4vLyAoYXQgeW91ciBvcHRpb24pIGFueSBsYXRlciB2ZXJzaW9uLlxuLy9cbi8vIFRyaWFkQmFsYW5jZSBpcyBkaXN0cmlidXRlZCBpbiB0aGUgaG9wZSB0aGF0IGl0IHdpbGwgYmUgdXNlZnVsLFxuLy8gYnV0IFdJVEhPVVQgQU5ZIFdBUlJBTlRZOyB3aXRob3V0IGV2ZW4gdGhlIGltcGxpZWQgd2FycmFudHkgb2Zcbi8vIE1FUkNIQU5UQUJJTElUWSBvciBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRS4gIFNlZSB0aGVcbi8vIEdOVSBHZW5lcmFsIFB1YmxpYyBMaWNlbnNlIGZvciBtb3JlIGRldGFpbHMuXG4vL1xuLy8gWW91IHNob3VsZCBoYXZlIHJlY2VpdmVkIGEgY29weSBvZiB0aGUgR05VIEdlbmVyYWwgUHVibGljIExpY2Vuc2Vcbi8vIGFsb25nIHdpdGggVHJpYWRCYWxhbmNlLiAgSWYgbm90LCBzZWUgPGh0dHBzOi8vd3d3LmdudS5vcmcvbGljZW5zZXMvPi5cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbnZhciB2ZWNNb2R1bGUgPSByZXF1aXJlKFwiLi4vanMvdmVjLm1vZHVsZS5qc1wiKTtcblxuLy8gdmVjLWxhLWZwIHBsYWNlcyBuaWNlIG5hbWVzIGluIGEgbWVtYmVyIG5hbWVkIFwidmVjXCJcbnZhciB2ZWMgPSB2ZWNNb2R1bGUudmVjO1xuXG4vLyBUT0RPOiBUaGlzIG5lZWRzIHRvIGJlIHNjYWxlZCEhIVxuZnVuY3Rpb24gbWVhbih3dGMpIHtcbiAgcmV0dXJuIHZlYy5zY2FsZSh3dGMubGVuZ3RoLHZlYy5hZGRBbGwod3RjKSk7XG59XG5cbi8vIFRlc3QgdGhhdCB0aGUgcG9pbnQgaXMgdmVydGljYWxseSBvcmllbnRlZCwgb3JpZ2luLWNlbnRlcmVkIGVxdWlsYXRlcmFsIHRyaWFuZ2xlLlxuZnVuY3Rpb24gaXNDZW50ZXJlZEVxdWlsYXRlcmFsKHd0Yykge1xuICBsZXQgZDAgPSB2ZWMuZGlzdCh3dGNbMF0sd3RjWzFdKTtcbiAgbGV0IGQxID0gdmVjLmRpc3Qod3RjWzFdLHd0Y1syXSk7XG4gIGxldCBkMiA9IHZlYy5kaXN0KHd0Y1syXSx3dGNbMF0pO1xuXG4gIHJldHVybiB2ZWMubmVhcigxZS01LFswLDBdLG1lYW4od3RjKSlcbiAgLy8gVGhpcmQgcG9pbnQgdmVydGljYWwuLlxuICAgICYmICgod3RjWzJdWzBdID09IDApICYmICh3dGNbMl1bMV0gPiAwKSlcbiAgLy8gZXF1aWxhdGVyYWxcbiAgICAmJiB2ZWMuc2NhbGFyTmVhcigxZS01LGQwLGQxKSAmJiB2ZWMuc2NhbGFyTmVhcigxZS01LGQxLGQyKSAmJiB2ZWMuc2NhbGFyTmVhcigxZS01LGQyLGQxKTtcbn1cblxuXG5cbi8vIFRoaXMgaXMgYSB0aW55IHNldCBvZiByb3V0aW5lcyBpbnNwaXJlZCBieSB2ZWMtbGEtZnAsIHdoaWNoXG4vLyBkb2VzIG5vdCBjdXJyZW50bHkgaGFuZGxlIDNkIHZlY3RvcnMuLi5hIHJlYWwgaGVybyB3b3VsZFxuLy8gZnVsbHkgZXh0ZW5kIHRoYXQgcGFja2FnZSB0byBzdXBwb3J0IDJkIGFuZCAzZCB2ZWN0b3JzLlxuXG5sZXQgdjNBZGQgPSBmdW5jdGlvbiB2M0FkZChhLGIpIHtcbiAgcmV0dXJuIFthWzBdK2JbMF0sYVsxXStiWzFdLGFbMl0rYlsyXV07XG59XG5sZXQgdjNTdWIgPSBmdW5jdGlvbiB2M1N1YihhLGIpIHtcbiAgcmV0dXJuIFthWzBdLWJbMF0sYVsxXS1iWzFdLGFbMl0tYlsyXV07XG59XG5sZXQgdjNNYW5oYXR0YW5EaXN0YW5jZSA9IGZ1bmN0aW9uIHYzTWFuaGF0dGFuRGlzdGFuY2UoYSxiKSB7XG4gIHJldHVybiBNYXRoLmFicyhhWzBdLWJbMF0pICsgTWF0aC5hYnMoYVsxXS1iWzFdKSArIE1hdGguYWJzKGFbMl0tYlsyXSk7XG59XG5sZXQgdjNkaXN0ID0gZnVuY3Rpb24gdjNkaXN0KGEsYikge1xuICByZXR1cm4gTWF0aC5zcXJ0KChhWzBdLWJbMF0pKioyICsgKGFbMV0tYlsxXSkqKjIgKyAoYVsyXS1iWzJdKSoqMik7XG59XG5sZXQgdjNtYWcgPSBmdW5jdGlvbiB2M21hZyh2KSB7XG4gIHJldHVybiBNYXRoLnNxcnQoKHZbMF0pKioyICsgKHZbMV0pKioyICsgKHZbMl0pKioyKTtcbn1cbmxldCB2M3NjYWxlID0gZnVuY3Rpb24gdjNzY2FsZShzYyx2KSB7XG4gIHJldHVybiBbc2MqdlswXSxzYyp2WzFdLHNjKnZbMl1dO1xufVxubGV0IHYzbm9ybWFsaXplID0gZnVuY3Rpb24gbm9ybWFsaXplKHYpIHtcbiAgcmV0dXJuIHYzc2NhbGUoIDEvdjNjLm1hZyh2KSx2KTtcbn1cblxudmFyIHYzYyA9IHtcbiAgYWRkOiB2M0FkZCxcbiAgc3ViOiB2M1N1YixcbiAgbWFuaGF0dGFuOiB2M01hbmhhdHRhbkRpc3RhbmNlLFxuICBkaXN0OiB2M2Rpc3QsXG4gIG1hZzogdjNtYWcsXG4gIHNjYWxlOiB2M3NjYWxlLFxuICBub3JtYWxpemU6IHYzbm9ybWFsaXplXG59XG5cblxuXG5jb25zdCBTUVJUMyA9IE1hdGguc3FydCgzKTtcblxuZnVuY3Rpb24gR2V0UmF5VG9MaW5lU2VnbWVudEludGVyc2VjdGlvbihyYXlPcmlnaW4scmF5RGlyZWN0aW9uLHBvaW50MSxwb2ludDIpXG57XG4gIC8vIFRoaXMgY29kZSBmcm9tIGhlcmU6IGh0dHBzOi8vc3RhY2tvdmVyZmxvdy5jb20vcXVlc3Rpb25zLzE0MzA3MTU4L2hvdy1kby15b3UtY2hlY2stZm9yLWludGVyc2VjdGlvbi1iZXR3ZWVuLWEtbGluZS1zZWdtZW50LWFuZC1hLWxpbmUtcmF5LWVtYW5hdGluXG4gIC8vIE5vdGUgdGhpcyByb3V0aW5lIHNlZW1zIHRvIGRlcGVuZCBvbiB0aGUgY2hpcmFsaXR5IG9mIHRoZSBwb2ludHM7IHBvc3NpYmx5IGl0IG9ubHkgY291bnRzIGFuIGFwcHJvYWNoIGZyb20gb25lIHNpZGUuXG4gIGNvbnN0IHJkbiA9IHZlYy5ub3JtKHJheURpcmVjdGlvbik7XG4gIGNvbnN0IHYxID0gdmVjLnN1YihyYXlPcmlnaW4scG9pbnQxKTtcbiAgY29uc3QgdjIgPSB2ZWMuc3ViKHBvaW50Mixwb2ludDEpO1xuXG4gIGNvbnN0IHYzID0gWy1yZG5bMV0scmRuWzBdXTtcbiAgY29uc3QgZG90ID0gdmVjLmRvdCh2Mix2Myk7XG5cbiAgaWYgKHZlYy5zY2FsYXJOZWFyKDFlLTUsZG90LDApKVxuICAgIHJldHVybiBudWxsO1xuXG4gIGNvbnN0IHQxID0gdmVjLnBlcnBkb3QodjIsdjEpIC8gZG90O1xuICBjb25zdCB0MiA9IHZlYy5kb3QodjEsdjMpIC8gZG90O1xuXG4gIGlmICh0MSA+PSAwLjAgJiYgKHQyID49IDAuMCAmJiB0MiA8PSAxLjApKSB7XG4gICAgcmV0dXJuIFt2ZWMuYWRkKHJheU9yaWdpbix2ZWMuc2NhbGUodDEscmRuKSksdDFdO1xuICB9XG4gIHJldHVybiBudWxsO1xufVxuXG4vLyBUaGlzIGlzIHRoZSBmdW5kYW1lbnRhbCBtYXRoIGJlaGluZCBhIFRyaWFkQmFsYW5jZSBEaWFncmFtLlxuLy8gSSB3aWxsIHByb2JhYmx5IHB1dCB0aGlzIGluIExhVGVYIHNvb24uXG4vLyBBIGJhbGFuY2UgZGlhZ3JhbSBpcyBhIHdheSBvZiBjaG9vc2luZyBhIHVuaXQgdmVjdG9yXG4vLyBvZiBuIGRpbWVuc2lvbnMgKG1vc3QgdXNlZnVsbHkgMykgZnJvbSBhIG4tZ29uLlxuLy8gSWYgbiA+IDMsIGl0IGlzIG5vdCBwb3NzaWJsZSB0byBjb21wbGV0ZSBtYXAgYWxsIHBvaW50cy5cbi8vIFRoZSBmdW5kYW1lbnRhbCBtYXRoIG9mIGEgVHJpYWRCYWxhbmNlIGRpYWdyYW0gaXMgdG8gY29udmVydFxuLy8gYmV0d2VlbiBhIHBvaW50IG9uIHRoZSAyLWRpbWVuc2lvbmFsIG4tZ29uICh0aGUgXCJyZXByZXNlbnRhdGlvblwiKSBmcm9tXG4vLyB0aGUgXCJiYWxhbmNlXCIgLS0gYSBuLWRpbWVuc2lvbmFsIHZlY3Rvci5cbi8vIENhbGwgdGhlIGZ1bmN0aW9uIHRoYXQgcHJvZHVjZXMgdGhlIHJlcGVzZW50YXRpb24gZnJvbSBhIHVuaXQgdmVjdG9yXG4vLyBcInJcIiBhbmQgdGhlIGZ1bmN0aW9uIHRoYXQgcHJvZHVjZXMgdGhlIGJhbGFuY2UgdmVjdG9yIGZyb20gdGhlIHJlc3ByZXNlbnRhdGlvbiBcImJcIi5cbi8vIERlc2lkZXJhdG9yIG9mIHRoZXNlIGZ1bmN0aW9ucyBhcmU6XG4vLyBXZSB3YW50IHRoZW0gdG8gYmUgaW52ZXJzaW9ucyBvZiBlYWNoIG90aGVyLlxuLy8gV2Ugd2FudCB0aGUgY2VudGVyIG9mIHRoZSByZXByZXNlbnRhdGlvbiB0byBtYXAgdG8gYSBiYWxhbmNlZCB2ZWN0b3IuXG4vLyBBIHJlcHJlc2VudGF0aW9uIGF0IHRoZSB2ZXJ0ZXggaXMgYSB2ZWN0b3IgaGF2aW5nIGEgc2luZ2xlIDEgYW5kIHRoZSByZXN0IHplcm9zLlxuLy8gQXMgd2UgY2hhbmdlIHRoZSBhbmdsZSBiZXR3ZW4gdGhlIG9yaWdpbiBhbmQgdGhlIHBvaW50IGluIGEgcmVwcmVzZW50YXRpb24gdG93YXJkIGEgdmVydGV4LFxuLy8gdGhlIHZhbHVlIG9mIHRoYXQgdmVydGV4IHNob3VsZCBpbmNyZWFzZS5cbi8vIEFzIHdlIG1vdmUgYWxvbmcgc3VjaCBhIGxpbmUsIHdlIHNob3VsZCBub3QgY2hhbmdlIHRoZSByZWxhdGl2ZSBwcm9wb3J0aW9uIG9mIHRoZVxuLy8gb3RoZXIgdmFsdWVzLiAodGhpcyBpcyB2YWd1ZSkuXG4vLyBJdCBzaG91bGQgYmUgZWFzeSB0byBjb21wdXRlIGFuZCBleHBsYWluIChhdCBsZWFzdCBnZW9tZXRyaWNhbGx5Lilcbi8vIEkgbm93IGJlbGlldmUgdGhpcyBzaG91bGQgYWxsb3cgYSBcIm5vcm1cIiB0byBiZSBwYXNzZWQgaW5cbi8vIGFzIGEgZnVuY3Rpb24uIEkgdGhpbmsgdGhlIEwxIG5vcm0gaXMgcHJvYmFibHkgYmV0dGVyIHRoYW4gTDJcbi8vIG5vcm0gZm9yIHNvbWUgZnVuY3Rpb25zLCBidXQgaXQgY2FuIGJlIG9wdGlvbmFsLlxuLy8gSXQgaXMgZXNzZW50aWFsIHRoYXQgdGhpcyBmdW5jdGlvbiBiZSBpbnZlcnRpYmxlLlxuXG5mdW5jdGlvbiBMMU5PUk0odikge1xuICBsZXQgcyA9IHYzYy5tYW5oYXR0YW4oWzAsMCwwXSx2KTtcbiAgcmV0dXJuIHYzYy5zY2FsZSgxL3Msdik7XG59XG5mdW5jdGlvbiBMMk5PUk0odikge1xuICByZXR1cm4gdjNjLm5vcm1hbGl6ZSh2KTtcbn1cbmZ1bmN0aW9uIEwxTEVOR1RIKHYpIHtcbiAgcmV0dXJuIHYzYy5tYW5oYXR0YW4oWzAsMCwwXSx2KTtcbn1cbmZ1bmN0aW9uIEwyTEVOR1RIKHYpIHtcbiAgcmV0dXJuIHYzYy5tYWcodik7XG59XG5cbnZhciBMMSA9IFtMMU5PUk0sTDFMRU5HVEhdO1xudmFyIEwyID0gW0wyTk9STSxMMkxFTkdUSF07XG5cblxuLy8gVW5kZXIgYXNzdW1wdGlvbiBvZiBhbiB1cHdhcmQgZmFjaW5nXG4vLyBlcXVpbGF0ZXJhbCB0cmlhbmdsZSwgcmV0dXJuIHRoZSBlZGdlXG4vLyBpbnRlcnNlY3RlZCBieSB0aGUgcmF5IGZyb20gdGhlIG9yaWdpbiB0byB0cCxcbi8vIHdpdGhvdXQgZGVwZW5kZW5jZSBvbiB2ZWN0b3IgbGlicmFyaWVzLlxuLy8gVGhlIHJldHVybiB2YWx1ZSBpcyBhbiBhcnJheTpcbi8vIGVtcHR5IGlmZiB0cCBpcyB0aGUgb3JpZ2luXG4vLyBjb250YWluaW5nIHR3byB2YWx1ZXMgaWYgaXQgaXMgYSB2ZXJ0ZXhcbi8vIGNvbnRhaW5pbmcgb25lIHZhbHVlIG90aGVyd2lzZS5cbi8vIFRoZSBlZGdlcyBhcmUgbnVtYmVyZWQgYW50aWNsb2Nrd2lzZSxcbi8vIHdpdGggdGhlIHplcm9ldGggZWRnZSBhdCB0aGUgYm90dG9tLlxuLy8gSWYgdHdvIGVkZ2VzIGFyZSByZXR1cm5lZCwgdGhpcyByb3V0aW5lXG4vLyByZXR1cm5zIHRoZW0gaW4gc29ydGVkIG9yZGVyLlxuZnVuY3Rpb24gZXFFZGdlQWxnZWJyYWljYWxseSh3dGMscCkge1xuICBpZiAodmVjLnNjYWxhck5lYXIoMWUtNSxwWzBdLDApKSB7XG4gICAgaWYgKHZlYy5zY2FsYXJOZWFyKDFlLTUscFsxXSwwKSkge1xuICAgICAgcmV0dXJuIFtdO1xuICAgIH0gZWxzZSBpZiAocFsxXSA+IDApIHtcbiAgICAgIHJldHVybiBbMSwyXTtcbiAgICB9IGVsc2Uge1xuICAgICAgcmV0dXJuIFswXTtcbiAgICB9XG4gIH0gZWxzZSB7XG4gICAgbGV0IG0gPSBwWzFdL3BbMF07XG4gICAgbGV0IG0xID0gLVNRUlQzLzM7XG4gICAgbGV0IG0yID0gU1FSVDMvMztcbiAgICBpZiAoKHBbMF0gPiAwKSAmJiAodmVjLnNjYWxhck5lYXIoMWUtNSxtLG0xKSkpIHJldHVybiBbMCwxXTtcbiAgICBpZiAoKHBbMF0gPiAwKSAmJiAobSA+IG0xKSkgcmV0dXJuIFsxXTtcbiAgICBpZiAoKHBbMF0gPiAwKSAmJiAobSA8IG0xKSkgcmV0dXJuIFswXTtcblxuICAgIGlmICgocFswXSA8IDApICYmICh2ZWMuc2NhbGFyTmVhcigxZS01LG0sbTIpKSkgcmV0dXJuIFswLDJdO1xuICAgIGlmICgocFswXSA8IDApICYmIChtIDwgbTIpKSByZXR1cm4gWzJdO1xuICAgIGlmICgocFswXSA8IDApICYmIChtID4gbTIpKSByZXR1cm4gWzBdO1xuICB9XG59XG5cblxuLy8gSGVyZSB3ZSByZXR1cm4gdGhlIHBvaW50IG9uIHRoZSBlZGdlXG4vLyB3aGVyZSB0aGUgcmF5IGZyb20gdGhlIG9yaWdpbiB0byB0cCBpbnRlcnNlY3RzXG4vLyB0aGUgdHJpYW5nbGUuXG4vLyBUaGUgbWF0aCBoZXJlIGlzIGNyZWF0ZWQgb24gdGhlIGFzc3VtcHRpb25cbi8vIG9mIHRoZSB0cmlhbmdsZSBiZWluZyBwZXJmZWN0bHkgZXF1aWxhdGVyYWwgd2l0aFxuLy8gdGhlIGNlbnRyb2lkIGF0IHRoZSBvcmlnaW4uIFRoaXMgYWxsb3dzIHVzIHRvXG4vLyBzb2x2ZSBzaW11bHRhbmVvdXMgZXF1YXRpb25zIHRvIGZpbmQgdGhlIHBvaW50cy5cbi8vIFRoaXMgaXMgYW4gYWx0ZXJuYXRpdmUgdG8gdmVjdG9yLWJhc2VkIG1ldGhvZHNcbi8vIHRoYXQgYXJlIG9mIGNvdXJzZSBpbiB0aGUgZW5kIHNpbWxhciwgYnV0IHdlIHRha2Vcbi8vIGFkdmFudGFnZSBvZiBrbm93biBzbG9wZXMgdG8gbWFrZSBpdCBmYXN0ZXIgYW5kIHNpbXBsZXIuXG5mdW5jdGlvbiBlcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSh3dGMsdHApIHtcbiAgLy8gd2Ugc2hvdWxkIHByb2JhYmx5IGNoZWNrIGVxdWlsYXRlcmFsaXR5IGFuZCBvcmllbnRhdGlvbiBoZXJlXG4gIGxldCBlcyA9IGVxRWRnZUFsZ2VicmFpY2FsbHkod3RjLHRwKTtcbiAgaWYgKGVzLmxlbmd0aCA9PSAwKSByZXR1cm4gbnVsbDsgLy8gdHAgaXMgdGhlIG9yaWdpblxuICBpZiAoZXMubGVuZ3RoID09IDIpIHsgLy8gd2UgaGl0IGEgdmVydGV4LCBidXQgd2hpY2ggb25lP1xuICAgIGlmICgoZXNbMF0gPT0gMCkgJiYgKGVzWzFdID09IDEpKSB7XG4gICAgICByZXR1cm4gd3RjWzFdO1xuICAgIH0gZWxzZVxuICAgIGlmICgoZXNbMF0gPT0gMCkgJiYgKGVzWzFdID09IDIpKSB7XG4gICAgICByZXR1cm4gd3RjWzBdO1xuICAgIH0gZWxzZVxuICAgIGlmICgoZXNbMF0gPT0gMSkgJiYgKGVzWzFdID09IDIpKSB7XG4gICAgICByZXR1cm4gd3RjWzJdO1xuICAgIH1cbiAgfSBlbHNlIHsgLy8gbm93IHdlIGRvIGEgY2FzZSBzcGxpdFxuICAgIGxldCB4cCA9IHRwWzBdO1xuICAgIGxldCB5cCA9IHRwWzFdO1xuICAgIGxldCBhID0gdmVjLmRpc3Qod3RjWzBdLHd0Y1sxXSk7XG4gICAgbGV0IEIgPSBhICogU1FSVDMvNjtcbiAgICBpZiAodmVjLnNjYWxhck5lYXIoMWUtNSx4cCwwKSkge1xuICAgICAgcmV0dXJuICh5cCA+IDApID8gd3RjWzJdIDogWzAsLUJdO1xuICAgIH1cbiAgICBsZXQgbSA9IHlwL3hwO1xuICAgIGlmIChlc1swXSA9PSAwKSB7XG4gICAgICByZXR1cm4gWy1CL20sLUJdO1xuICAgIH0gZWxzZSBpZiAoZXNbMF0gPT0gMSkge1xuICAgICAgbGV0IHkgPSBhIC8gKDMgKigxL1NRUlQzICsgMS9tKSk7XG4gICAgICBsZXQgeCA9IHkgLyBtO1xuICAgICAgcmV0dXJuIFt4LHldO1xuICAgIH0gZWxzZSBpZiAoZXNbMF0gPT0gMikge1xuICAgICAgbGV0IHkgPSBhIC8gKDMgKigxL1NRUlQzIC0gMS9tKSk7XG4gICAgICBsZXQgeCA9IHkgLyBtO1xuICAgICAgcmV0dXJuIFt4LHldO1xuICAgIH1cbiAgfVxufVxuXG5cbmZ1bmN0aW9uIGdldEVkZ2VBbmRQb2ludCh3dGMscCkge1xuXG4gIC8vIElmIHdlIGFyZSBjZW50ZXJlZCwgdmVydGljYWwsIHBvaW50aW5nIHVwLCBhbmQgZXF1aWxhdGVyYWwsXG4gIC8vIHdlIGNhbiB1c2UgdGhlIG1vcmUgZWZmaWNpZW50IGFsZ29yaXRobS5cbiAgaWYgKGlzQ2VudGVyZWRFcXVpbGF0ZXJhbCh3dGMpKVxuICB7XG4gICAgLy8gdGhpcyBtYXkgcmV0dXJuIHR3bywgYnV0IHdlIGNhbiBqdXN0IHRha2UgdGhlIGZpcnN0XG4gICAgcmV0dXJuIFtlcUVkZ2VBbGdlYnJhaWNhbGx5KHd0YyxwKVswXSxcbiAgICAgICAgICAgIGVxUG9pbnRPbkVkZ2VBbGdlYnJhaWNhbGx5KHd0YyxwKV07XG4gIH1cbiAgZWxzZVxuICB7XG4gICAgdmFyIHBvaW50X29uX2VkZ2U7XG4gICAgdmFyIGZlX2lkeCA9IC0xOyAvLyBpbmRleCBvZiB0aGUgZmlyc3QgZWRnZSB3ZSBpbnRlcnNlY3RcbiAgICBmb3IodmFyIGkgPSAwOyBpIDwgMyAmJiBmZV9pZHggPCAwOyBpKyspIHtcbiAgICAgIHZhciByID0gR2V0UmF5VG9MaW5lU2VnbWVudEludGVyc2VjdGlvbihbMCwwXSxwLHd0Y1tpXSx3dGNbKGkrMSkgJSAzXSk7XG4gICAgICBpZiAociAhPSBudWxsKSB7IC8vIGlmIG51bGwsIHRoZSByYXkgZGlkIG5vdCBpbnRlcnNlY3QgdGhlIGVkZ2VcbiAgICAgICAgZmVfaWR4ID0gaTtcbiAgICAgICAgcG9pbnRfb25fZWRnZSA9IHJbMF07IC8vIFRoZSBmaXJzdCBjb21wLiBvZiByZXR1cm4gdmFsdWUgaXMgaW50ZXJzZWN0aW9uXG4gICAgICB9XG4gICAgfVxuICAgIHJldHVybiBbZmVfaWR4LHBvaW50X29uX2VkZ2VdO1xuICB9XG59XG5cblxuLy8gdHAgaXMgYSBwb2ludCBpbiB0aGUgMi1kaW1lbnNpb25hbCB0cmlhbmdsZSBzcGFjZVxuLy8gd3RjIGFyZSB0aGUgdGhyZWUgdmVydGljZXMgb2YgYW4gZXFpbGF0ZXJhbCB0cmlhbmdsZSB3aG9zZSBjZW50cm9pZCBpcyB0aGUgb3JpZ2luXG4vLyBMWG5vcm1fYW5kX2xlbmd0aCBpcyBhIHBhaXIgb2YgZnVuY3Rpb25zIHRvIHRvIG5vcm1hbGl6ZSBhIHZlY3RvciBhbmQgY29tcHV0ZSB0aGUgbGVuZ3RoXG4vLyByZXR1cm4gdGhlIGNvcnJlc3BvbmRpbmcgMy12ZWN0b3IgaW4gdGhlIGF0dHJpYnV0ZSBzcGFjZVxuZnVuY3Rpb24gVHJpYWRCYWxhbmNlMnRvMyhwLHd0YyxMWG5vcm1fYW5kX2xlbmd0aCA9IEwyKSB7XG4gIGxldCBMWG5vcm1hbGl6ZSA9IExYbm9ybV9hbmRfbGVuZ3RoWzBdO1xuXG4gIGlmICh2ZWMuc2NhbGFyTmVhcigxZS01LHZlYy5tYWcocCksMCkpIHtcbiAgICByZXR1cm4gTFhub3JtYWxpemUoWzEsMSwxXSk7XG4gIH1cblxuICAvLyBOb3cgd2Ugd2FudCB0byBkbyBhIGxpbmVhciBpbnRlcnBvbGF0aW9uIG9mIGhvdyBmYXIgd2UgYXJlIGZyb20gYW4gZWRnZSxcbiAgLy8gYnV0IGFsc28gaG93IGZhciB0aGUgcHJvamVjdGlvbiB0byB0aGUgZWRnZSBpcyBiZXR3ZWVuIHRoZSB2ZXJ0aWNlcy5cbiAgLy8gV2UgbXVzdCBmaXJzdCBkZWNpZGUgd2hpY2ggZWRnZXMgdGhlIGxpbmUgZnJvbSB0aGUgb3JpZ24gdG8gcCBpbnRlcnNlY3RzLlxuICAvLyBJZiBpdCBpbnRlcnNlY3RzIHR3byBzZWdtZW50cywgdGhlbiBpdCBpcyBhaW1lZCBhdCBhIHZlcnRleC5cbiAgbGV0IFtmZV9pZHgscG9pbnRfb25fZWRnZV0gPSBnZXRFZGdlQW5kUG9pbnQod3RjLHApO1xuXG4gIC8vIG5vdyBwb2ludF9vbl9lZGdlIGlzIGEgcG9pbnQgb24gZWRnZSBmZV9pZHguXG4gIGNvbnN0IHRvdGFsX2Rpc3RhbmNlX3RvX2VkZ2UgPSB2ZWMuZGlzdChbMCwwXSxwb2ludF9vbl9lZGdlKTtcblxuICAvLyBJZiB0aGUgcG9pbnQgaXMgb3V0c2lkZSB0aGUgdHJpYW5nbGUsIHdlIGNsYW1wICh0cnVuY2F0ZSBpZiBuZWVkZWQpXG4gIC8vIGl0J3MgbGVuZ3RoIHNvIHRoYXQgaXQgaXMgcHJlY2lzZWx5IG9uIHRoZSBlZGdlLlxuICBjb25zdCBwYyA9IHZlYy5jbGFtcE1hZygwLHRvdGFsX2Rpc3RhbmNlX3RvX2VkZ2UscCk7XG5cbiAgY29uc3QgZGlzdGFuY2VfdG9fcF9vX2UgPSB2ZWMuZGlzdChwYyxwb2ludF9vbl9lZGdlKTtcbiAgdmFyIHJhdGlvX3BfdG9fZWRnZSA9ICBkaXN0YW5jZV90b19wX29fZS90b3RhbF9kaXN0YW5jZV90b19lZGdlO1xuXG4gIGxldCBiYWwgPSB2M2Muc2NhbGUocmF0aW9fcF90b19lZGdlLFxuICAgICAgICAgICAgICAgICAgICAgIExYbm9ybWFsaXplKFsxLDEsMV0pKTtcblxuICAvLyBOb3cgdGhlIHJlbWFpbmRlciBvZiB0aGUgY29udHJpYnV0aW9uXG4gIC8vIHRvIHRoZSB1bml0IHZlY3RvciBzaG91bGQgY29tZSBmcm9tIHRoZSB0d29cbiAgLy8gcG9pbnRzIG9uIHRoZSBlZGdlLCBpbiBsaW5lYXIgcHJvcG9ydGlvbi5cbiAgLy8gVGhlc2UgY29vcmRpbmF0ZXMgYXJlIGZlX2lkeCBhbmQgKGZlX2lkeCsxKSAlIDMuXG4gIGNvbnN0IGQxID0gdmVjLmRpc3Qod3RjW2ZlX2lkeF0scG9pbnRfb25fZWRnZSk7XG4gIGNvbnN0IGQyID0gdmVjLmRpc3Qod3RjWyhmZV9pZHgrMSkgJSAzXSxwb2ludF9vbl9lZGdlKTtcblxuICBsZXQgdnMgPSBbMCwwLDBdO1xuICB2c1tmZV9pZHhdID0gZDI7XG4gIHZzWyhmZV9pZHgrMSkgJSAzXSA9IGQxO1xuXG4gIGxldCBpbWIgPSB2M2Muc2NhbGUoMSAtIHJhdGlvX3BfdG9fZWRnZSxMWG5vcm1hbGl6ZSh2cykpO1xuXG4gIHJldHVybiB2M2MuYWRkKGltYixiYWwpO1xufVxuXG4vLyB2ZWMgaXMgYSAzLXZlY3RvciBpbiB0aGUgYXR0cmlidXRlIHNwYWNlXG4vLyB3dGMgYXJlIHRoZSB0aHJlZSB2ZXJ0aWNlcyBvZiBhbiBlcWlsYXRlcmFsIHRyaWFuZ2xlIHdob3NlIGNlbnRyb2lkIGlzIHRoZSBvcmlnaW5cbi8vIExYbm9ybV9hbmRfbGVuZ3RoIGlzIGEgcGFpciBvZiBmdW5jdGlvbnMgdG8gdG8gbm9ybWFsaXplIGEgdmVjdG9yIGFuZCBjb21wdXRlIHRoZSBsZW5ndGhcbi8vIHJldHVybiB0aGUgY29ycmVzcG9uZGluZyAyLXZlY3RvciBpbiB0aGUgdHJpYW5nbGUgc3BhY2VcbmZ1bmN0aW9uIGludmVydFRyaWFkQmFsYW5jZTJ0bzModix3dGMsTFhub3JtX2FuZF9sZW5ndGggPSBMMikge1xuICBsZXQgbGVuZ3RoID0gTFhub3JtX2FuZF9sZW5ndGhbMV07XG5cbiAgbGV0IG1pbiA9IE1hdGgubWluKE1hdGgubWluKHZbMF0sdlsxXSksdlsyXSk7XG5cbiAgbGV0IGltYiA9IFt2WzBdIC0gbWluLHZbMV0gLSBtaW4sdlsyXSAtIG1pbl07XG4gIGxldCBiYWwgPSB2M2Muc3ViKHYsaW1iKTtcbiAgLy8gTm93IHRoYXQgd2UgaGF2ZSBiYWxhbmNlLCB3ZSBuZWVkIHRvIGNvbXB1dGUgaXQncyBsZW5ndGgsXG4gIC8vIHdoaWNoIGlzIGRlcGVuZGVudCBvbiB0aGUgbm9ybSB3ZSBjaG9zZSFcblxuICBsZXQgaW1iX3IgPSBsZW5ndGgoaW1iKTtcbiAgbGV0IGJhbF9yID0gbGVuZ3RoKGJhbCk7XG5cbiAgLy8gTm93IHdlIGhhdmUgdGhlIHJhdGlvcy4gV2UgbmVlZCB0byBkZXRlcm1pbmUgdGhlIGRpcmVjdGlvbi5cbiAgLy8gVGhpcyBpcyBhIGZ1bmN0aW9uIG9mIHRoZSBpbWJhbGFuY2UgdmVjdG9yLiBXZSBjYW4gZGV0ZXJtaW5lXG4gIC8vIHdoaWNoIHNpZGUgd2UgYXJlIG9uLCBhbmQgdGhlbiBjb21wdXRlIG91ciBwb3NpdGlvbiBhbG9uZyB0aGF0XG4gIC8vIHRvIGRldGVybWluZSBhIHBvaW50IG9uIHRoZSB0cmlhbmdsZSwgYW5kIHRoZW4gbXVsdGlwbHkgYnkgdGhlIGltYl9yXG4gIC8vIHRvIG9idGFpbiB0aGUgYWN0dWFsIHBvaW50LlxuICAvLyBBdCBsZWFzdCBvbmUgdmFsdWUgb2YgaW1iIHdpbGwgYmUgemVyby5cbiAgdmFyIGZyb21fdix0b192LHJhdGlvO1xuICAvLyB0aGUgcG9pbnRzIGFyZSBPUFBPU0lURSB0aGUgemVyb1xuICAvLyByYXRpbyB3aWxsIGJlIHRoZSByYXRpbyBhbG9uZyB0aGUgdHJpYW5nbGUgZWRnZVxuICAvLyBpdCByZXF1aXJlcyBhIGxpdHRsZSB0aG91Z2h0IHRvIHVuZGVyc3RhbmQgd2hpY2hcbiAgLy8gb2YgdGhlIG90aGVyIHBvaW50cyBzaG91bGQgYmUgdGhlIFwiZnJvbV92XCIgYW5kIHRoZSBcInRvX3ZcIlxuICAvLyBmb3IgdGhlIGludGVycG9sYXRpb24gd2hpY2ggb2NjdXJzIGxhdGVyLlxuICBsZXQgcyA9IGltYlswXSArIGltYlsxXSArIGltYlsyXTsgLy8gb25lIG9mIHRoZXNlIGlzIGFsd2F5cyB6ZXJvLlxuICBpZiAoaW1iWzBdID09IDApIHtcbiAgICBmcm9tX3YgPSB3dGNbMl07XG4gICAgdG9fdiA9IHd0Y1sxXTtcbiAgICByYXRpbyA9IGltYlsxXS9zO1xuICB9IGVsc2UgaWYgKGltYlsxXSA9PSAwKSB7XG4gICAgZnJvbV92ID0gd3RjWzBdO1xuICAgIHRvX3YgPSB3dGNbMl07XG4gICAgcmF0aW8gPSBpbWJbMl0vcztcbiAgfSBlbHNlIGlmIChpbWJbMl0gPT0gMCkge1xuICAgIGZyb21fdiA9IHd0Y1sxXTtcbiAgICB0b192ID0gd3RjWzBdO1xuICAgIHJhdGlvID0gaW1iWzBdL3M7XG4gIH1cblxuICAvLyBUaGUgcG9pbnQgb24gdGhlIHRyaWFuZ2xlIGlzIGJ5IGNvbnN0cnVjdGlvblxuICAvLyBvbiBvbmUgZWRnZSBvZiB0aGUgdHJpYW5nbGUuXG4gIGNvbnN0IG9uVHJpYW5nbGUgPSB2ZWMubGVycChmcm9tX3YsdG9fdixyYXRpbyk7XG4gIC8vIG5vdyBvblRyaWFuZ2xlIGlzIGEgcG9pbnQgb24gdGhlIHRyaWFuZ2xlXG4gIC8vIG5vdywgaGF2aW5nIGZvdW5kIHRoYXQgd2UgaW50ZXJwb2xhdGUgYSByYXlcbiAgLy8gdG8gaXQgb2YgbGVuZ3RoIGltYl9yLi4uXG4gIHJldHVybiB2ZWMubGVycChbMCwwXSxvblRyaWFuZ2xlLGltYl9yKTtcbn1cblxubW9kdWxlLmV4cG9ydHMgPSB7XG4gIFRyaWFkQmFsYW5jZTJ0bzM6IFRyaWFkQmFsYW5jZTJ0bzMsXG4gIGludmVydFRyaWFkQmFsYW5jZTJ0bzM6IGludmVydFRyaWFkQmFsYW5jZTJ0bzMsXG4gIGVxUG9pbnRPbkVkZ2VBbGdlYnJhaWNhbGx5OiBlcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSxcbiAgZXFFZGdlQWxnZWJyYWljYWxseTogZXFFZGdlQWxnZWJyYWljYWxseSxcbiAgR2V0UmF5VG9MaW5lU2VnbWVudEludGVyc2VjdGlvbiA6IEdldFJheVRvTGluZVNlZ21lbnRJbnRlcnNlY3Rpb24sXG4gIEwxTEVOR1RIOiBMMUxFTkdUSCxcbiAgTDJMRU5HVEg6IEwyTEVOR1RILFxuICBMMTogTDEsXG4gIEwyOiBMMn07XG4iLCIvLyBDb3B5cmlnaHQgMjAxOSwgUm9iZXJ0IEwuIFJlYWRcbi8vIFRoaXMgZmlsZSBpcyBwYXJ0IG9mIFRyaWFkQmFsYW5jZS5cbi8vXG4vLyBUcmlhZEJhbGFuY2UgaXMgZnJlZSBzb2Z0d2FyZTogeW91IGNhbiByZWRpc3RyaWJ1dGUgaXQgYW5kL29yIG1vZGlmeVxuLy8gaXQgdW5kZXIgdGhlIHRlcm1zIG9mIHRoZSBHTlUgR2VuZXJhbCBQdWJsaWMgTGljZW5zZSBhcyBwdWJsaXNoZWQgYnlcbi8vIHRoZSBGcmVlIFNvZnR3YXJlIEZvdW5kYXRpb24sIGVpdGhlciB2ZXJzaW9uIDMgb2YgdGhlIExpY2Vuc2UsIG9yXG4vLyAoYXQgeW91ciBvcHRpb24pIGFueSBsYXRlciB2ZXJzaW9uLlxuLy9cbi8vIFRyaWFkQmFsYW5jZSBpcyBkaXN0cmlidXRlZCBpbiB0aGUgaG9wZSB0aGF0IGl0IHdpbGwgYmUgdXNlZnVsLFxuLy8gYnV0IFdJVEhPVVQgQU5ZIFdBUlJBTlRZOyB3aXRob3V0IGV2ZW4gdGhlIGltcGxpZWQgd2FycmFudHkgb2Zcbi8vIE1FUkNIQU5UQUJJTElUWSBvciBGSVRORVNTIEZPUiBBIFBBUlRJQ1VMQVIgUFVSUE9TRS4gIFNlZSB0aGVcbi8vIEdOVSBHZW5lcmFsIFB1YmxpYyBMaWNlbnNlIGZvciBtb3JlIGRldGFpbHMuXG4vL1xuLy8gWW91IHNob3VsZCBoYXZlIHJlY2VpdmVkIGEgY29weSBvZiB0aGUgR05VIEdlbmVyYWwgUHVibGljIExpY2Vuc2Vcbi8vIGFsb25nIHdpdGggVHJpYWRCYWxhbmNlLiAgSWYgbm90LCBzZWUgPGh0dHBzOi8vd3d3LmdudS5vcmcvbGljZW5zZXMvPi5cblxuLy8gSXQgaXMgZXh0cmVtZWx5IHZhbHVhYmxlIHRvIGJlIGFibGUgdG8gY29tcHV0ZSB0aGVcbi8vIGludmVyc2Ugb2YgdGhlIFRyaWFkQmFsYW5jZSBhbGdvcml0aG0gZm9yIHRlc3RpbmcsXG4vLyBhbHRob3VnaCB0aGlzIGlzIG5vdCBxdWl0ZSBhIHRydWUgaW52ZXNpb24gYmVjYXVzZVxuLy8gcG9pbnRzIG91dHNpZGUgdGhlIHRyaWFuZ2xlIGFyZSBicm91Z2h0IHRvIHRoZSB0cmlhbmdsZVxuLy8gZWRnZS5cblxuXCJ1c2Ugc3RyaWN0XCI7XG5cbnZhciBtID0gcmVxdWlyZShcIi4uL3NyYy9UcmlhZEJhbGFuY2VNYXRoLmpzXCIpO1xudmFyIHZlY01vZHVsZSA9IHJlcXVpcmUoXCIuLi9qcy92ZWMubW9kdWxlLmpzXCIpO1xuXG4vLyB2ZWMtbGEtZnAgcGxhY2VzIG5pY2UgbmFtZXMgaW4gYSBtZW1iZXIgbmFtZWQgXCJ2ZWNcIlxudmFyIHZlYyA9IHZlY01vZHVsZS52ZWM7XG5cbmZ1bmN0aW9uIGNvbXB1dGVGaW52ZXJzZUYod3RjLG5vcm0scCkge1xuICAgcmV0dXJuIG0uaW52ZXJ0VHJpYWRCYWxhbmNlMnRvMyhcbiAgICAgbS5UcmlhZEJhbGFuY2UydG8zKHAsd3RjLG5vcm0pLFxuICAgICB3dGMsXG4gICAgIG5vcm0pO1xufVxuXG5mdW5jdGlvbiB0ZXN0X2VxUG9pbnRPbkFsZ2VicmFpY2FsbHkod3RjKSB7XG4gIHZhciByb2YgPSBtLmVxUG9pbnRPbkVkZ2VBbGdlYnJhaWNhbGx5KHd0YyxbMCwwXSk7XG4gIGNvbnNvbGUuYXNzZXJ0KHJvZiA9PSBudWxsKTtcblxuICBjb25zdCBoYWxmdmVydCA9IHZlYy5zY2FsZSgxLzIsd3RjWzBdKTtcbiAgdmFyIHIwID0gbS5lcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSh3dGMsaGFsZnZlcnQpO1xuICBjb25zb2xlLmFzc2VydCh2ZWMubmVhcigxZS01LHIwLHd0Y1swXSkpO1xuXG4gIHZhciByMSA9IG0uZXFQb2ludE9uRWRnZUFsZ2VicmFpY2FsbHkod3RjLHd0Y1sxXSk7XG4gIGNvbnNvbGUuYXNzZXJ0KHZlYy5uZWFyKDFlLTUscjEsd3RjWzFdKSk7XG5cbiAgdmFyIHIyID0gbS5lcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSh3dGMsd3RjWzJdKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLm5lYXIoMWUtNSxyMix3dGNbMl0pKTtcblxuICAvLyBzbG9wZSBvZiAxXG4gIHZhciByb25lID0gbS5lcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSh3dGMsWzEsMV0pO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuY29saW5lYXIod3RjWzFdLHd0Y1syXSxyb25lKSk7XG5cbiAgLy8gc2xvcGUgb2YgLTFcbiAgdmFyIHJuZWdfb25lID0gbS5lcVBvaW50T25FZGdlQWxnZWJyYWljYWxseSh3dGMsWy0xLDFdKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLmNvbGluZWFyKHd0Y1swXSx3dGNbMl0scm5lZ19vbmUpKTtcblxuICAvLyBhbG1vc3Qgc3RyYWlnaHQgZG93blxuICB2YXIgcmRvd24gPSBtLmVxUG9pbnRPbkVkZ2VBbGdlYnJhaWNhbGx5KHd0YyxbMC4wMSwtMV0pO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuY29saW5lYXIocmRvd24sd3RjWzBdLHd0Y1sxXSkpO1xuXG4gIHJldHVybiB0cnVlO1xufVxuXG5mdW5jdGlvbiB0ZXN0X2VxRWRnZUFsZWJyYWljYWxseSh3dGMpIHtcbiAgdmFyIHJvZiA9IG0uZXFFZGdlQWxnZWJyYWljYWxseSh3dGMsWzAsMF0pO1xuICBjb25zb2xlLmFzc2VydChyb2YubGVuZ3RoID09IDApO1xuXG4gIHZhciByMCA9IG0uZXFFZGdlQWxnZWJyYWljYWxseSh3dGMsd3RjWzBdKTtcbiAgY29uc29sZS5hc3NlcnQocjAubGVuZ3RoID09IDIpO1xuXG4gIHZhciByMSA9IG0uZXFFZGdlQWxnZWJyYWljYWxseSh3dGMsd3RjWzFdKTtcbiAgY29uc29sZS5hc3NlcnQocjEubGVuZ3RoID09IDIpO1xuXG4gIHZhciByMiA9IG0uZXFFZGdlQWxnZWJyYWljYWxseSh3dGMsd3RjWzJdKTtcbiAgY29uc29sZS5hc3NlcnQocjIubGVuZ3RoID09IDIpO1xuXG4gIC8vIHNsb3BlIG9mIDFcbiAgdmFyIHJvbmUgPSBtLmVxRWRnZUFsZ2VicmFpY2FsbHkod3RjLFsxLDFdKTtcbiAgY29uc29sZS5hc3NlcnQocm9uZS5sZW5ndGggPT0gMSk7XG4gIGNvbnNvbGUuYXNzZXJ0KHJvbmVbMF0gPT0gMSk7XG5cbiAgdmFyIHJuZWdfb25lID0gbS5lcUVkZ2VBbGdlYnJhaWNhbGx5KHd0YyxbLTEsMV0pO1xuICBjb25zb2xlLmFzc2VydChybmVnX29uZS5sZW5ndGggPT0gMSk7XG4gIGNvbnNvbGUuYXNzZXJ0KHJuZWdfb25lWzBdID09IDIpO1xuXG4gIHZhciByZG93biA9IG0uZXFFZGdlQWxnZWJyYWljYWxseSh3dGMsWzAuMDEsLTFdKTtcbiAgY29uc29sZS5hc3NlcnQocmRvd24ubGVuZ3RoID09IDEpO1xuICBjb25zb2xlLmFzc2VydChyZG93blswXSA9PSAwKTtcbiAgcmV0dXJuIHRydWU7XG59XG5cblxuZnVuY3Rpb24gdGVzdEdldFJheVRvTGluZVNlZ21lbnRJbnRlcnNlY3Rpb24od3RjKSB7XG4gIGxldCBybyA9IFswLDBdO1xuICBsZXQgcmQgPSBbMSwxXTtcbiAgbGV0IHAxID0gWzAsMTBdO1xuICBsZXQgcDIgPSBbMTAsMF07XG4gIHZhciByID0gbS5HZXRSYXlUb0xpbmVTZWdtZW50SW50ZXJzZWN0aW9uKHJvLHJkLHAxLHAyKVswXTtcbiAgY29uc29sZS5hc3NlcnQoclswXSA9PSByWzFdKTtcbiAgdmFyIHIgPSBtLkdldFJheVRvTGluZVNlZ21lbnRJbnRlcnNlY3Rpb24ocm8scmQscDIscDEpWzBdO1xuICBjb25zb2xlLmFzc2VydChyWzBdID09IHJbMV0pO1xuXG4gIHZhciByZDEgPSBbOTQuMTAxNTYyNSwtMzYuMzYzMjgxMjVdO1xuICBsZXQgYzAgPSB3dGNbMF07XG4gIGxldCBjMSA9IHd0Y1sxXTtcbiAgbGV0IGMyID0gd3RjWzJdO1xuXG4gIHZhciByMTIgPSBtLkdldFJheVRvTGluZVNlZ21lbnRJbnRlcnNlY3Rpb24ocm8scmQxLGMxLGMyKTtcbiAgY29uc29sZS5hc3NlcnQocjEyICE9IG51bGwpO1xuXG4gIHtcbiAgICB2YXIgZG93biA9IG0uR2V0UmF5VG9MaW5lU2VnbWVudEludGVyc2VjdGlvbihybyxbMCwtNDk5OV0sYzAsYzEpO1xuICAgIGNvbnNvbGUuYXNzZXJ0KHZlYy5zY2FsYXJOZWFyKDFlLTUsZG93blswXVsxXSx3dGNbMF1bMV0pKTtcbiAgfVxuICByZXR1cm4gdHJ1ZTtcbn1cblxuXG5mdW5jdGlvbiB0ZXN0VHJpYWRCYWxhbmNlMnRvMyh3dGMpIHtcbiAgbGV0IHB2ID0gbS5UcmlhZEJhbGFuY2UydG8zKFszMDAwMCwzMF0sd3RjLG0uTDEpO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuc2NhbGFyTmVhcigxZS01LG0uTDFMRU5HVEgocHYpLDEpKTtcblxuICBsZXQgcHl2ID0gbS5UcmlhZEJhbGFuY2UydG8zKFswLHd0Y1syXVsxXV0sd3RjLG0uTDEpO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuc2NhbGFyTmVhcigxZS01LG0uTDFMRU5HVEgocHl2KSwxKSk7XG4gIHJldHVybiB0cnVlO1xufVxuXG5mdW5jdGlvbiB0ZXN0T3JpZ2luQW5kVmVydGljZXMod3RjKSB7XG4gIC8vIFRoZSBvcmlnaW4gc2hvdWxkIHJldHVybiBhIHBlcmZlY3RseSBiYWxhbmNlZCB2ZWN0b3JcbiAgbGV0IG92ID0gbS5UcmlhZEJhbGFuY2UydG8zKFswLDBdLHd0YyxtLkwxKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSxvdlswXSwxLzMpKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSxvdlsxXSwxLzMpKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSxvdlsyXSwxLzMpKTtcblxuICAvLyBBIHZlcnRleCBzaG91bGQgcmV0dXJuIGEgdmVjdG9yIHdpdGggYSAxIGluIGV4YWN0bHkgMSBwb3NpdGlvblxuICBjb25zb2xlLmFzc2VydCh2ZWMuc2NhbGFyTmVhcigxZS01LG0uVHJpYWRCYWxhbmNlMnRvMyh3dGNbMF0sd3RjLG0uTDEpWzBdLDEpKTtcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSxtLlRyaWFkQmFsYW5jZTJ0bzMod3RjWzFdLHd0YyxtLkwxKVsxXSwxKSk7XG4gIGNvbnNvbGUuYXNzZXJ0KHZlYy5zY2FsYXJOZWFyKDFlLTUsbS5UcmlhZEJhbGFuY2UydG8zKHd0Y1syXSx3dGMsbS5MMSlbMl0sMSkpO1xuXG4gIHJldHVybiB0cnVlO1xufVxuXG4vLyBJdCBpcyBleHRyZW1lbHkgdmFsdWFibGUgdG8gYmUgYWJsZSB0byBjb21wdXRlIHRoZVxuLy8gaW52ZXJzZSBvZiB0aGUgVHJpYWRCYWxhbmNlIGFsZ29yaXRobSBmb3IgdGVzdGluZyxcbi8vIGFsdGhvdWdoIHRoaXMgaXMgbm90IHF1aXRlIGEgdHJ1ZSBpbnZlc2lvbiBiZWNhdXNlXG4vLyBwb2ludHMgb3V0c2lkZSB0aGUgdHJpYW5nbGUgYXJlIGJyb3VnaHQgdG8gdGhlIHRyaWFuZ2xlXG4vLyBlZGdlLlxuZnVuY3Rpb24gY29tcHV0ZUZpbnZlcnNlRih3dGMsbm9ybSxwKSB7XG4gICByZXR1cm4gbS5pbnZlcnRUcmlhZEJhbGFuY2UydG8zKFxuICAgICBtLlRyaWFkQmFsYW5jZTJ0bzMocCx3dGMsbm9ybSksXG4gICAgIHd0YyxcbiAgICAgbm9ybSk7XG59XG5cbmZ1bmN0aW9uIHRlc3RJbnZlcnNpb24od3RjKSB7XG4gIC8vIFRoaXMgaW5zdXJlcyB3ZSBhcmUgd2l0aGluIHRoZSB0cmlhbmdsZVxuICBsZXQgZCA9IHZlYy5kaXN0KHd0Y1swXSx3dGNbMV0pO1xuICBsZXQgcCA9IFtkLzgsZC84XTtcblxuICBsZXQgdnBjID0gdmVjLnN1YihcbiAgICBjb21wdXRlRmludmVyc2VGKHd0YyxtLkwxLHApLFxuICAgIHApO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuc2NhbGFyTmVhcigxZS01LHZlYy5tYWcodnBjKSwwKSk7XG5cbiAgbGV0IHB5ID0gWzAsd3RjWzJdWzFdXTtcbiAgbGV0IHZweWMgPSB2ZWMuc3ViKFxuICAgIGNvbXB1dGVGaW52ZXJzZUYod3RjLG0uTDEscHkpLFxuICAgIHB5KTtcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSx2ZWMubWFnKHZweWMpLDApKTtcbiAgcmV0dXJuIHRydWU7XG59XG5cbmZ1bmN0aW9uIHRlc3RJbnZlcnNpb25PdXRzaWRlKHd0Yykge1xuICAvLyBUaGlzIGluc3VyZXMgd2UgYXJlIG91dHNpZGUgdGhlIHRyaWFuZ2xlXG4gIGxldCBkID0gdmVjLmRpc3Qod3RjWzBdLHd0Y1sxXSk7XG5cbiAgdmFyIHZwX2ludiA9IGNvbXB1dGVGaW52ZXJzZUYod3RjLG0uTDEsW2QsZF0pO1xuICAvLyB3ZSBub3cgd2FudCB0byB0ZXN0IHRoYXQgdGhlIGludmVzaW9uIGhhcyBhIHNsb3BlIG9mIDFcbiAgY29uc29sZS5hc3NlcnQodmVjLnNjYWxhck5lYXIoMWUtNSx2cF9pbnZbMV0vdnBfaW52WzBdLDEpKTtcbiAgY29uc29sZS5hc3NlcnQodnBfaW52WzFdID4gMCk7XG4gIHJldHVybiB0cnVlO1xufVxuXG5mdW5jdGlvbiB0ZXN0SW52ZXJzaW9uTmVnYXRpdmVZKHd0Yykge1xuICAvLyBUaGlzIGluc3VyZXMgd2UgYXJlIHdpdGhpbiB0aGUgdHJpYW5nbGVcbiAgbGV0IGQgPSB2ZWMuZGlzdCh3dGNbMF0sd3RjWzFdKTtcbiAgbGV0IHAgPSBbMCwtZC80XTtcbiAgdmFyIHZwX2ludiA9IGNvbXB1dGVGaW52ZXJzZUYod3RjLG0uTDEscCk7XG4gIC8vIHRlc3QgbGVuZ3RoIGhlcmVcbiAgdmFyIHZwYyA9IHZlYy5zdWIodnBfaW52LHApO1xuICBjb25zb2xlLmFzc2VydCh2ZWMuc2NhbGFyTmVhcigxZS01LHZlYy5tYWcodnBjKSwwKSk7XG4gIHJldHVybiB0cnVlO1xufVxuXG4vLyBUZXN0IHZpYSBhIGNpcmNsZSBjb21wbGV0ZWx5IHdpdGhpbiB0aGUgdHJpYW5nbGUsXG4vLyB0aHVzIGV4ZXJjaXNpbmcgYWxsIGFuZ2xlcyAoaW5zdGlnYXRlZCBieSBhIGJ1Zy4pXG5mdW5jdGlvbiB0ZXN0SW52ZXJzaW9uV2l0aEFDaXJjbGUod3RjLE5PUk0pIHtcblxuICBsZXQgd3RjcCA9IHd0Y1sxXTtcbiAgbGV0IHJhZGl1cyA9IHZlYy5tYWcod3RjcCkvMztcbiAgbGV0IG4gPSAxMztcbiAgbGV0IG9uZV90aGlydGVlbnRoID0gMiAqIE1hdGguUEkgLyAxMztcbiAgZm9yKHZhciBpID0gMDsgaSA8IG47IGkrKykge1xuICAgIGxldCB4ID0gTWF0aC5zaW4oaSpvbmVfdGhpcnRlZW50aCk7XG4gICAgbGV0IHkgPSBNYXRoLmNvcyhpKm9uZV90aGlydGVlbnRoKTtcbiAgICAvLyBub3cgeCx5IGlzIGEgcG9pbnQgb24gYSBjaXJjbGUgd2l0aGluIHRoZSB0cmluYWdsZVxuICAgIC8vIHdlIHdpbGwgbWFrZSBzdXJlIHRoZSBiYWxhbmNlIGZ1bmN0aW9uIGludmVydHNcbiAgICAvLyB0byBnaXZlIHRoZSBmdW5jdGlvbiBiYWNrIHRvIHVzLlxuICAgIHtcbiAgICAgIGxldCBwID0gW3gseV07XG4gICAgICB2ZWMuc2NhbGUocmFkaXVzLHApO1xuICAgICAgdmFyIHZwX2ludiA9ICBjb21wdXRlRmludmVyc2VGKHd0YyxOT1JNLHApO1xuICAgICAgdmFyIHZweWMgPSB2ZWMuc3ViKHZwX2ludixwKTtcbiAgICAgIGNvbnNvbGUuYXNzZXJ0KHZlYy5zY2FsYXJOZWFyKDFlLTUsdmVjLm1hZyh2cHljKSwwKSk7XG4gICAgfVxuICB9XG4gIHJldHVybiB0cnVlO1xufVxuXG4vLyB3dGMgaXMgdGhlIFdPUkxEX1RSSUFOR0xFX0NPT1JEUywgYW4gYXJyYXkgb2YgdGhyZWUgVmVjdG9yMiBvYmplY3RzLlxuZnVuY3Rpb24gdGVzdEVxdWlsYXRlcmFsRnVuY3Rpb25zKHd0Yykge1xuICB2YXIgYWxsVHJ1ZSA9IDE7XG4gIGFsbFRydWUgJj0gdGVzdF9lcUVkZ2VBbGVicmFpY2FsbHkod3RjKTtcbiAgYWxsVHJ1ZSAmPSB0ZXN0X2VxUG9pbnRPbkFsZ2VicmFpY2FsbHkod3RjKTtcbiAgcmV0dXJuIGFsbFRydWU7XG59XG5cbmZ1bmN0aW9uIHRlc3RBbGxUcmlhZEJhbGFuY2UodXB3YXJkLHd0Yykge1xuICB2YXIgYWxsVHJ1ZSA9IDE7XG4gIGlmICh1cHdhcmQpIHtcbiAgICBhbGxUcnVlICY9IHRlc3RFcXVpbGF0ZXJhbEZ1bmN0aW9ucyh3dGMpO1xuICAgIGFsbFRydWUgJj0gdGVzdEdldFJheVRvTGluZVNlZ21lbnRJbnRlcnNlY3Rpb24od3RjKTtcbiAgfVxuICBhbGxUcnVlICY9IHRlc3RPcmlnaW5BbmRWZXJ0aWNlcyh3dGMpO1xuICBhbGxUcnVlICY9IHRlc3RUcmlhZEJhbGFuY2UydG8zKHd0Yyk7XG4gIGFsbFRydWUgJj0gdGVzdEludmVyc2lvbih3dGMpO1xuICBhbGxUcnVlICY9IHRlc3RJbnZlcnNpb25OZWdhdGl2ZVkod3RjKTtcbiAgYWxsVHJ1ZSAmPSB0ZXN0SW52ZXJzaW9uT3V0c2lkZSh3dGMpO1xuICBsZXQgd3RjcCA9IHd0Y1sxXTtcblxuXG4gIGxldCBzbWFsbF9yYWRpdXMgPSB2ZWMubWFnKHd0Y3ApLzM7XG4gIGFsbFRydWUgJj0gdGVzdEludmVyc2lvbldpdGhBQ2lyY2xlKHd0YyxtLkwxLHNtYWxsX3JhZGl1cyk7XG4gIGFsbFRydWUgJj0gdGVzdEludmVyc2lvbldpdGhBQ2lyY2xlKHd0YyxtLkwyLHNtYWxsX3JhZGl1cyk7XG4gIGxldCBsYXJnZV9yYWRpdXMgPSB2ZWMubWFnKHd0Y3ApKjM7XG4gIGFsbFRydWUgJj0gdGVzdEludmVyc2lvbldpdGhBQ2lyY2xlKHd0YyxtLkwxLGxhcmdlX3JhZGl1cyk7XG4gIGFsbFRydWUgJj0gdGVzdEludmVyc2lvbldpdGhBQ2lyY2xlKHd0YyxtLkwyLGxhcmdlX3JhZGl1cyk7XG4gIHJldHVybiBhbGxUcnVlO1xufVxuXG5tb2R1bGUuZXhwb3J0cyA9IHtcbiAgdGVzdEFsbFRyaWFkQmFsYW5jZTogdGVzdEFsbFRyaWFkQmFsYW5jZVxufTtcbiJdfQ==
