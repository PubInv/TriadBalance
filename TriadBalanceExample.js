//   Copyright 2019, Robert L. Read
//  
//   This file is part of TriadBalance.
//  
//   TriadBalance is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//  
//   TriadBalance is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with TriadBalance.  If not, see <https://www.gnu.org/licenses/>.

"use strict";

var CUR_BALANCE;
var CUR_TRIANGLE_COORDS;
var W;
var H;
var Whalf;
var Hhalf;
var NORM_TO_USE;
var LABELS;
var TRIAD_WORLD_TRIANGLE;

function set_norm_to_use(norm) {
  NORM_TO_USE = norm;
}

function set_labels(labels) {
  LABELS = labels;
}

function get_world_triangle(svg_id) {
  var svg_elt = document.getElementById(svg_id);

  // Now the question is can we make this work with any width element?
  let W = svg_elt.clientWidth;
  let H = svg_elt.clientHeight;
  var min = Math.min(W,H);

  // This could be a parameter...
  let SIDE_LENGTH_PIXEL = min * (3/5);
  
  const TRIANGLE_WIDTH = 1;
  const TRIANGLE_HEIGHT = Math.sqrt(3)/2;
  
  const SIDE_LENGTH_HEIGHT = SIDE_LENGTH_PIXEL * TRIANGLE_HEIGHT;
  const BASE = -(1/3) * SIDE_LENGTH_HEIGHT;  

  let wtc_vector = [new THREE.Vector2(-SIDE_LENGTH_PIXEL/2,BASE),
                    new THREE.Vector2(SIDE_LENGTH_PIXEL/2,BASE),
                    new THREE.Vector2(0,BASE+SIDE_LENGTH_HEIGHT)];
  return wtc_vector;
}

function vpy(y) { return (-y); }
function vpx(x) { return (x); }

function clear_markers(svg) {
  var x = document.getElementsByClassName("triad-marker");
    for(var i = 0; i < x.length; i++) {
      x[i].parentNode.removeChild(x[i]);
      console.log("removed");
    }
}
function rerender_marker(svg,bal_vec) {
  if (bal_vec) {
    // We have to find the current marker and remove it..
    clear_markers(svg);
    let tri_point = invertTriadBalance2to3(bal_vec,TRIAD_WORLD_TRIANGLE,NORM_TO_USE);
    let point = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
    point.setAttributeNS(null, 'cx', vpx(tri_point.x));
    point.setAttributeNS(null, 'cy', vpy(tri_point.y));
    point.setAttributeNS(null, 'r', 4);
    point.setAttributeNS(null,"class","triad-marker");
    point.ISMARKER = true;
    svg.appendChild(point);
  }
}

function render_svg(svg,fs_ratio_to_height) {
  console.assert(fs_ratio_to_height);
  let fs = H * fs_ratio_to_height;

  function append_text(svg,id,class_name,x,y,dy,text) {
    var newText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    newText.setAttributeNS(null,"x",x);      
    newText.setAttributeNS(null,"y",y);
    newText.setAttributeNS(null,"dy",dy);
    newText.setAttributeNS(null,"font-size",fs);          
    newText.setAttributeNS(null,"class",class_name);
    newText.setAttributeNS(null,"id",id);      
    newText.appendChild(document.createTextNode(text));
    svg.appendChild(newText);
  }

  // Note: if we wished to depend on jQueryUI or d3, there are other solutions:
  // https://stackoverflow.com/questions/3674265/is-there-an-easy-way-to-clear-an-svg-elements-contents
  while (svg.lastChild) {
    svg.removeChild(svg.lastChild);
  }
  
  var polygon = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
  polygon.setAttributeNS(null,"id","triad-balance-triangle");  
  svg.appendChild(polygon);

  for (let i = 0; i < 3; i++) {
    var point = svg.createSVGPoint();
    point.x = vpx(TRIAD_WORLD_TRIANGLE[i].x);
    point.y = vpy(TRIAD_WORLD_TRIANGLE[i].y);
    polygon.points.appendItem(point);
  }

  // These are ugly, they should really be computed from the text.
  // In fact, since this does not change, the whole thing could go into
  // HTML and css more profitably.
  // fs is the fontsize in pixels; hopefully we caculate this.
  function render_labels(svg,vertices,d_labels,fs) {
    append_text(svg,"vertex-2","triad-vertices",
                vpx(TRIAD_WORLD_TRIANGLE[2].x),vpy(TRIAD_WORLD_TRIANGLE[2].y),-fs/2,
                d_labels[2]
               );
    append_text(svg,"vertex-0","triad-vertices",
                vpx(TRIAD_WORLD_TRIANGLE[0].x),vpy(TRIAD_WORLD_TRIANGLE[0].y),fs,
                d_labels[0]              
               );
    append_text(svg,"vertex-1","triad-vertices",
                vpx(TRIAD_WORLD_TRIANGLE[1].x),vpy(TRIAD_WORLD_TRIANGLE[1].y),fs,
                d_labels[1]
               );
  }
  if (!fs) {
    debugger;
  }
  render_labels(svg,TRIAD_WORLD_TRIANGLE,LABELS,fs);


  // this is the center of the triangle...
  let origin = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
  origin.setAttributeNS(null, 'cx', vpx(0));
  origin.setAttributeNS(null, 'cy', vpy(0));
  origin.setAttributeNS(null, 'r', 2);
  origin.setAttributeNS(null,"id","triangle_origin");        
  svg.appendChild(origin);
  
//  rerender_marker(svg,CUR_POINT);
}


// This is tricky because click events on an SVG
// depend on which object inside the SVG are hit.
// We don't really want to do that, we have
// created a global triangle space. A solution is to
// use screen coordinates, but only to compute a difference.
function clicked(evt,container_id,fs,svg,labels,click_callback) {
  var br = document.getElementById(container_id).getBoundingClientRect();
  var x = evt.clientX - br.left;
  var y = evt.clientY - br.top;
  // x and y are in the coordinates of the
  // SVG system; we need to convert them
  // to the coordinates of our triangle.
  var xc = x + -Whalf;  
  var yc = -(y + -Hhalf);
  var triangle_coords = 

  // Note, we could balance and invert here to make sure we are inside the triangle!
      CUR_TRIANGLE_COORDS = new THREE.Vector2(xc,yc);
  console.log("COORDS",CUR_TRIANGLE_COORDS);
  CUR_BALANCE = TriadBalance2to3(CUR_TRIANGLE_COORDS,TRIAD_WORLD_TRIANGLE,NORM_TO_USE);

  var CUR_POINT = invertTriadBalance2to3(CUR_BALANCE,TRIAD_WORLD_TRIANGLE,NORM_TO_USE);
  
  rerender_marker(svg,CUR_BALANCE);  
  click_callback(CUR_TRIANGLE_COORDS,CUR_POINT,CUR_BALANCE);
}         

function initialize_triad_diagram(svg_id,norm_to_use,labels,click_callback,fs_ratio_to_height = (1/20)) {
  var svg_elt = document.getElementById(svg_id);

  // Now the question is can we make this work with any width element?
  function setSizeConstants(svg_elt) {
    W = svg_elt.clientWidth;
    H = svg_elt.clientHeight;
    Whalf = W/2;
    Hhalf = H/2;
    TRIAD_WORLD_TRIANGLE = get_world_triangle(svg_id);
  // This is convenient, but makes it hard for the client to
  // use this svg for their own purposes, which is probably okay
  // for this use...
    svg_elt.setAttribute("viewBox", `-${Whalf} -${Hhalf} ${W} ${H}`);
    render_svg(svg_elt,fs_ratio_to_height);
    rerender_marker(svg_elt,CUR_BALANCE);
  }
  
  set_labels(labels);
  set_norm_to_use(norm_to_use);

  setSizeConstants(svg_elt);

  svg_elt.addEventListener(
    "click",
    (evt) =>
      clicked(evt,svg_id,fs_ratio_to_height,svg_elt,labels,click_callback)
  );
  // This is insufficient; we need to recalculate the world triangle.
  // Maybe I should make it internal?
  window.addEventListener(
    "resize",
    (evt) => {
      setSizeConstants(svg_elt);
    });

}





