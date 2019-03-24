// Copyright 2019, Robert L. Read
// This file is part of TriBalance.
//
// TriBalance is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TriBalance is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TriBalance.  If not, see <https://www.gnu.org/licenses/>.




function GetRayToLineSegmentIntersection(rayOrigin,rayDirection,point1,point2)
{
// This code from here: https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin
// Note this routine seems to depend on the chirality of the points; possibly it only counts an approach from one side.

    var ro = new THREE.Vector2(rayOrigin.x,rayOrigin.y);
    var rd = new THREE.Vector2(rayDirection.x,rayDirection.y);
    rd.normalize();
    var p1 = new THREE.Vector2(point1.x,point1.y);
    var p2 = new THREE.Vector2(point2.x,point2.y);
    var v1 = new THREE.Vector2().subVectors(ro,p1);
    var v2 = new THREE.Vector2().subVectors(p2,p1);    

    // This is construct a perpendicular.
    var v3 = new THREE.Vector2(-rd.y,rd.x);
    var dot = v2.dot(v3);

    if (Math.abs(dot) < 0.000001)
            return null;

    var t1 = v2.cross(v1) / dot;
    
    var t2 = v1.dot(v3) /dot;    

    if (t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)) {
        rd.multiplyScalar(t1);
        return [ro.add(rd),t1];
    }
    return null;
}

// This is the fundamental math behind a TriBalance Diagram.
// I will probably put this in LaTeX soon.
// A balance diagram is a way of choosing a unit vector
// of n dimensions (most usefully 3) from a n-gon.
// If n > 3, it is not possible to complete map all points.
// The fundamental math of a TriBalance diagram is to convert
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

// tp is the a point in the 2-dimensional triangle space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 3-vector in the attribute space
function TriBalance2to3(tp,wtc,LXnorm_and_length = L2) {
    let LXnormalize = LXnorm_and_length[0];
    if (near(tp.lengthSq(),0,1e-5)) {
        return LXnormalize(new THREE.Vector3(1,1,1));
    }
    let p = new THREE.Vector2(tp.x,tp.y);
    
    let origin = new THREE.Vector2(0,0);    
    // Now we want to do a linear interpolation of how far we are from an edge,
    // but also how far the projection to the edge is between the vertices.
    // We must first decide which edges the line form the orign to p intersects.
    // If it intersects two segments, then it is aimed at a vertex.
    var point_on_edge; 
    var fe_idx = -1; // index of the first edge we intersect
    for(var i = 0; i < 3 && fe_idx < 0; i++) {
        var r = GetRayToLineSegmentIntersection(origin,p,wtc[i],wtc[(i +1) % 3]);
        if (r != null) { // if null, the ray did not intersect the edge
            fe_idx = i;
            point_on_edge = r[0]; // The first comp. of return value is intersection
        }
    }
    // now point_on_edge is a point on edge fe_idx.     

    let total_distance_to_edge = origin.distanceTo(point_on_edge);
    // If the point is outside the triangle, we clamp (truncate if needed)
    // it's length so that it is precisely on the edge.
    p.clampLength(0,total_distance_to_edge);

    let distance_to_p_o_e = p.distanceTo(point_on_edge);
    
    var ratio_p_to_edge =  distance_to_p_o_e/total_distance_to_edge;
    
    let bal = LXnormalize(new THREE.Vector3(1,1,1));
    bal.multiplyScalar(ratio_p_to_edge);

    // Now the remainder of the contribution
    // to the unit vector should come from the two
    // points on the edge, in linear proportion.
    // These coordinates are fe_idx and (fe_idx+1) % 3.
    var d1 = wtc[fe_idx].distanceTo(point_on_edge);
    var d2 = wtc[(fe_idx+1) % 3].distanceTo(point_on_edge);
    
    let vs = [0,0,0];
    vs[fe_idx] = d2;
    vs[(fe_idx+1) % 3] = d1;
    
    let imb = LXnormalize(new THREE.Vector3(vs[0],vs[1],vs[2]));
    imb.multiplyScalar(1 - ratio_p_to_edge);
    
    // now construct a balanced vector proportional
    // to the length from the edge to the p towards the axis
    // so that this be a unit vector if p is the origin.
    return new THREE.Vector3().add(imb).add(bal);
}

// vec is a 3-vector in the attribute space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 2-vector in the triangle space
function invertTriBalance2to3(vec,wtc,LXnorm_and_length = L2) {
    let length = LXnorm_and_length[1];
    let min = Math.min(Math.min(vec.x,vec.y),vec.z);
    let imb = new THREE.Vector3(vec.x - min, vec.y - min, vec.z - min);
    let bal = vec.clone();
    bal.sub(imb);
    // Now that we have balance, we need to compute it's length,
    // which is dependent on the norm we chose!

    let imb_r = length(imb);
    let bal_r = length(bal);
    console.assert(Math.abs((bal_r+imb_r) - 1) <   1e-5);

    // Now we have the ratios. We need to determine the direction.
    // this is a function of the imbalance vector. We could determine
    // which side we are on, and then compute our position along that
    // to determine a point on the triangle, and then multiply by the imb_r
    // to obtain the actual point.
    // At least one value of imb will be zero.
    var from,to,ratio;
    // the points are OPPOSITE the zero
    // ratio will be the ratio along the triangle edge
    // it requires a little thought to understand which
    // of the other points should be the "from" and the "to"
    // for the interpolation which occurs later.
    var s = imb.x + imb.y + imb.z; // one of these is always zero.
    if (imb.x == 0) {
        from = wtc[2];
        to = wtc[1];
        ratio = imb.y/s;
    } else if (imb.y == 0) {
        from = wtc[0];
        to = wtc[2];
        ratio = imb.z/s;        
    } else if (imb.z == 0) {
        from = wtc[1];
        to = wtc[0];
        ratio = imb.x/s;        
    }

    // The point on the triangle is by construction
    // on one edge of the triangle.
    var onTriangle = new THREE.Vector2();
    onTriangle.lerpVectors(from,to,ratio);
    // now onTriangle is a point on the triangle
    // now, having found that we interpolate a ray
    // to it of length imb_r...
    let origin = new THREE.Vector2(0,0);
    let inversion = new THREE.Vector2();
    inversion.lerpVectors(origin,
                          onTriangle,
                          imb_r);
    return inversion;
}

function testGetRayToLineSegmentIntersection(wtc) {
    let ro = new THREE.Vector2(0,0);
    let rd = new THREE.Vector2(1,1);
    let p1 = new THREE.Vector2(0,10);
    let p2 = new THREE.Vector2(10,0);
    var r = GetRayToLineSegmentIntersection(ro,rd,p1,p2)[0];
    console.assert(r.x == r.y);
    var r = GetRayToLineSegmentIntersection(ro,rd,p2,p1)[0];
    console.assert(r.x == r.y);

    var rd1 = new THREE.Vector2(94.1015625,-36.36328125);
    let c0 = new THREE.Vector2(wtc[0].x,wtc[0].y);
    let c1 = new THREE.Vector2(wtc[1].x,wtc[1].y);
    let c2 = new THREE.Vector2(wtc[2].x,wtc[2].y);        

    var r12 = GetRayToLineSegmentIntersection(ro,rd1,c1,c2);
    console.assert(r12 != null);
}


function testTriBalance2to3(wtc) {
    let p = new THREE.Vector2(30000,30);    
    let pv = TriBalance2to3(p,wtc,L1);
    console.assert(near(L1LENGTH(pv),1));
    let py = new THREE.Vector2(0,wtc[2].y);
    let pyv = TriBalance2to3(py,wtc,L1);
    console.assert(near(L1LENGTH(pyv),1));    
}

function testOriginAndVertices(wtc) {
    // The origin should return a perfectly balanced vector
    let o = new THREE.Vector2(0,0);
    let ov = TriBalance2to3(o,wtc,L1);
    console.assert(near(ov.x,1/3));
    console.assert(near(ov.y,1/3));
    console.assert(near(ov.z,1/3));    

    // A vertex should return a vector with a 1 in exactly 1 position
    {
        let p = new THREE.Vector2(wtc[0].x,wtc[0].y);
        let pv = TriBalance2to3(p,wtc,L1);
        console.assert(near(pv.x,1));
    }

    {
        let p = new THREE.Vector2(wtc[1].x,wtc[1].y);
        let pv = TriBalance2to3(p,wtc,L1);
        console.assert(near(pv.y,1));
    }

    {
        let p = new THREE.Vector2(wtc[2].x,wtc[2].y);
        let pv = TriBalance2to3(p,wtc,L1);
        console.assert(near(pv.z,1));
    }
}

function testInversion(wtc) {
    let p = new THREE.Vector2(30,30);    
    var vp = TriBalance2to3(p,wtc,L1);
    var vp_inv = invertTriBalance2to3(vp,wtc,L1);
    // test length here
    var vpc = vp_inv.clone();
    vpc.sub(p);
    console.assert(vpc.length() < 1e-4);
    
    let py = new THREE.Vector2(0,wtc[2].y);
    var vpy = TriBalance2to3(py,wtc,L1);    
    var vpy_inv = invertTriBalance2to3(vpy,wtc,L1);
    
    var vpyc = vpy_inv.clone();
    vpyc.sub(py);
    console.assert(vpyc.length() < 1e-4);
}

function testInversionNegativeY(wtc) {
    let p = new THREE.Vector2(0,-30);    
    var vp = TriBalance2to3(p,wtc,L1);
    var vp_inv = invertTriBalance2to3(vp,wtc,L1);
    // test length here
    var vpc = vp_inv.clone();
    vpc.sub(p);
    console.assert(vpc.length() < 1e-4);
}

// Test via a circle completely within the triangle,
// thus exercising all angles (instigated by a bug.)
function testInversionWithACircle(wtc,NORM) {

    let wtcp = new THREE.Vector2(wtc[1].x,
                                 wtc[1].y);
    let radius = wtcp.length()/3;    
    let n = 10;
    let thirty_six = 2 * Math.PI / 10;
    for(var i = 0; i < n; i++) {
        let x = Math.sin(i*thirty_six);
        let y = Math.cos(i*thirty_six);
        // now x,y is a point on a circle within the trinagle
        // we will make sure the balance function inverts
        // to give the function back to us.
        {
            let p = new THREE.Vector2(x,y);
            p.multiplyScalar(radius);        
            var vp = TriBalance2to3(p,wtc,NORM);
            var vp_inv = invertTriBalance2to3(vp,wtc,NORM);
            var vpyc = vp_inv.clone();
            vpyc.sub(p);
            console.assert(vpyc.length() < 1e-4);
        }
    }
}

// wtc is the WORLD_TRIANGLE_COORDS, an array of three Vector2 objects.
function testAllTriBalance(wtc) {
    testOriginAndVertices(wtc);    
    testGetRayToLineSegmentIntersection(wtc);
    testTriBalance2to3(wtc);
    testInversion(wtc);
    testInversionNegativeY(wtc);

    let wtcp = new THREE.Vector2(wtc[1].x,
                                 wtc[1].y);
    
    let small_radius = wtcp.length()/3;    
    testInversionWithACircle(wtc,L1,small_radius);
    testInversionWithACircle(wtc,L2,small_radius);    
    let large_radius = wtcp.length()*3;    
    testInversionWithACircle(wtc,L1,large_radius);
    testInversionWithACircle(wtc,L2,large_radius);    
}

