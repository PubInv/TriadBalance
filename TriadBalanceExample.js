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


var THREE_DIMENSIONS = { "MBS": ["Mind","Body","Spirit"],
                         "HHH": ["Head","Heart","Hands"],
                         "PSM": ["Physical","Social","Mental"],
                         "SEF": ["Strength","Endurance","Flexibility"],
                         "FSH": ["Father","Son","Holy Spirit"]
                       };
var CURRENT_D = "MBS";

var WORLD_TRIANGLE_COORDS;
const W = 500;
const H = 500;

const TRIANGLE_WIDTH = 1;
const TRIANGLE_HEIGHT = Math.sqrt(3)/2;
const SIDE_LENGTH_PIXEL = 300;
const SIDE_LENGTH_HEIGHT = SIDE_LENGTH_PIXEL * TRIANGLE_HEIGHT;
const BASE = -(1/3) * SIDE_LENGTH_HEIGHT;

function get_world_triangle() {
  let wtc_vector = [new THREE.Vector2(-SIDE_LENGTH_PIXEL/2,BASE),
                    new THREE.Vector2(SIDE_LENGTH_PIXEL/2,BASE),
                    new THREE.Vector2(0,BASE+SIDE_LENGTH_HEIGHT)];
  return wtc_vector;
}


// These function convert my abstract coordinates
// to svg coordinates.

var CUR_POINT;
var CUR_TRIANGLE_COORDS;

function append_text(svg,id,class_name,x,y,text) {
  var newText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
  newText.setAttributeNS(null,"x",x);      
  newText.setAttributeNS(null,"y",y);
  newText.setAttributeNS(null,"class",class_name);
  newText.setAttributeNS(null,"id",id);      
  newText.appendChild(document.createTextNode(text));
  svg.appendChild(newText);
}


function render_svg(svg,wtc) {

  function vpy(y) { return (-y); }
  function vpx(x) { return (x); }


  // Note: if we wished to depend on jQueryUI or d3, there are other solutions:
  // https://stackoverflow.com/questions/3674265/is-there-an-easy-way-to-clear-an-svg-elements-contents
  while (svg.lastChild) {
    svg.removeChild(svg.lastChild);
  }
  
  var polygon = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
  polygon.setAttributeNS(null,"id","triad-balance-triangle");  
  svg.appendChild(polygon);

  // These are the vertices....
  var array = [ [ vpx(wtc[0].x),vpy(wtc[0].y) ], 
                [ vpx(wtc[1].x),vpy(wtc[1].y) ],
                [ vpx(wtc[2].x),vpy(wtc[2].y) ] ];
  
  for (let i = 0; i < 3; i++) {
    value = array[i];
    var point = svg.createSVGPoint();
    point.x = value[0];
    point.y = value[1];
    polygon.points.appendItem(point);
  }
//  polygon.style.fill='lemonchiffon';

  // These are ugly, they should really be computed from the text.
  // In fact, since this does not change, the whole thing could go into
  // HTML and css more profitably.
  
  append_text(svg,"vertex-0","triad-vertices",array[2][0]-20,array[2][1]-5,
              THREE_DIMENSIONS[CURRENT_D][2]);
  append_text(svg,"vertex-1","triad-vertices",array[0][0]-40,array[0][1]+20,
              THREE_DIMENSIONS[CURRENT_D][0]);
  append_text(svg,"vertex-2","triad-vertices",array[1][0],array[1][1]+20,
              THREE_DIMENSIONS[CURRENT_D][1]
             );


  // this is the center of the triangle...
  let orign = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
  orign.setAttributeNS(null, 'cx', vpx(0));
  orign.setAttributeNS(null, 'cy', vpy(0));
  orign.setAttributeNS(null, 'r', 2);
  orign.setAttributeNS(null, 'style', 'fill: black; stroke: black; stroke-width: 1px;' );
  svg.appendChild(orign);
  
  function add_point(tri,c) {
    if (tri) {
      let point = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
      point.setAttributeNS(null, 'cx', vpx(tri.x));
      point.setAttributeNS(null, 'cy', vpy(tri.y));
      point.setAttributeNS(null, 'r', 4);
      point.setAttributeNS(null, 'style', 'fill: '+c+'; stroke: '+c+'; stroke-width: 1px;' );
      svg.appendChild(point);
      
    }
  }
  add_point(CUR_POINT,"red");
}


function main() {
  var svg = document.getElementById("create_svg");        
  render_svg(svg,WORLD_TRIANGLE_COORDS);

  function set_dimension_labels(cur) {
    $("#d0l").text(THREE_DIMENSIONS[cur][0] + ":");
    $("#d1l").text(THREE_DIMENSIONS[cur][1] + ":");
    $("#d2l").text(THREE_DIMENSIONS[cur][2] + ":");
  }
  function set_and_render(cur) {
    CURRENT_D = cur;
    set_dimension_labels(CURRENT_D);
    update_dimension_names();
    var svg = document.getElementById("create_svg");            
    render_svg(svg,WORLD_TRIANGLE_COORDS);        
  }
  
  $("#MBS").click(() => { set_and_render("MBS");});
  $("#HHH").click(() => { set_and_render("HHH");});
  $("#PSM").click(() => { set_and_render("PSM");});
  $("#SEF").click(() => { set_and_render("SEF");});
  $("#FSH").click(() => { set_and_render("FSH");});

  // Catch the norm change and recompute...
  $(':radio[name="norm"]').change(function() {
    setBalance(CUR_POINT);
  });
}


function getRadioValue(name) {
  var rates = document.getElementsByName(name);
  var rate_value;
  for(var i = 0; i < rates.length; i++){
    if(rates[i].checked){
      rate_value = i;
    }
  }
  return rate_value;
}


function setBalance(atp) {
  var norm_to_use = (getRadioValue("norm") == 0 ? L1 :L2);
  let bal = TriadBalance2to3(atp,WORLD_TRIANGLE_COORDS,norm_to_use);

  $( "#d0" ).text( (bal.x * 100).toFixed(0) +  "%" );
  $( "#d1" ).text( (bal.y * 100).toFixed(0) +  "%" );
  $( "#d2" ).text( (bal.z * 100).toFixed(0) +  "%" );
  return [bal.x,bal.y,bal.z];
}

// This is tricky because click events on an SVG
// depend on which object inside the SVG are hit.
// We don't really want to do that, we have
// created a global triangle space. A solution
// that doesn't force us to become dependent on the SVG model
// of objects rendered is to use screen coordinates.

function clicked(evt){
  var br = document.getElementById("container_to_have_global_coords_on_svg").getBoundingClientRect();
  var x = evt.originalEvent.clientX - br.left;
  var y = evt.originalEvent.clientY - br.top;
  // x and y are in the coordinates of the
  // SVG system; we need to convert them
  // to the coordinates of our triangle.
  var xc = x + -W/2;  
  var yc = -(y + -H/2);
  XC = xc;
  YC = yc;
  var triangle_coords = new THREE.Vector2(xc,yc);

  // Note, we could balance and invert here to make sure we are inside the trianble!
  CUR_TRIANGLE_COORDS = triangle_coords;
  var bal = setBalance(triangle_coords);
  var vec = new THREE.Vector3(bal[0],bal[1],bal[2]);
  
  var norm_to_use = (getRadioValue("norm") == 0 ? L1 :L2);
  var triangle_coords_inside_triangle = invertTriadBalance2to3(vec,WORLD_TRIANGLE_COORDS,norm_to_use);

  CUR_POINT = triangle_coords_inside_triangle;
  var svg = document.getElementById("create_svg");          
  render_svg(svg,WORLD_TRIANGLE_COORDS);
}         

$("#container_to_have_global_coords_on_svg").click(clicked);

