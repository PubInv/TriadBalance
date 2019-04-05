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

var m = require("./TriadBalanceMath.js");
var vec = require("../js/vec.module.js");

// The TriadBalanceState is stashed in the svg element object
// given us to insure uniqueness. I use the annoying long name "triad_balance_state"
// specifically to make the chance of name collision very small.

class TriadBalanceState {
  constructor(svg,ccb,lbls,
              twt,w,h,cb,ctc,
              ntu = m.L1,
              s_t_m_b = 7/10,
              fs_r_t_h = 1/20,
              ydpc = 10,
              pu = true) {

    // REQUIRED
    // The SVG HTML ELMENT
    this.SVG_ELT = svg;
    // callback to receive numeric data on click
    // This routine is called ccb(tp,tipi,bal) upon clicks
    // on the SVG elemement.
    // tp = The triangle coordinates where the click occured,
    // tpi = The point on the same line guaranteed to be "in" the triangle
    // bal = The normalized 3-vector representing the attributes.
    this.CLICK_CALLBACK = ccb;
    // An array of 3 titles for the 3 attributes
    this.LABELS = lbls;

    // THESES ARE COMPUTED VALUES
    // The World Triangle coordinates
    this.TRIAD_WORLD_TRIANGLE = twt;
    // The width of the SVG element
    this.W = w;
    // The Height  of the SVG element
    this.H = h;

    // These are convenient storage
    // The most recent "balance vector", a 3-attribute vector
    this.CUR_BALANCE = cb;
    // The most recent click point in 2-dimensional space of the SVG triangle
    this.CUR_TRIANGLE_COORDS = ctc;

    // OPTIONAL
    // Either the L1 or L2 norm (see TriadBalanceMath.js)
    this.NORM_TO_USE = ntu;
    // For confiruation, the size of the font relative to the total height
    this.FONT_SIZE_RATIO_TO_HEIGHT = fs_r_t_h;
    // The percent of the height to lower the orign to create a balanced appearance
    this.Y_DISP_PER_CENT = ydpc;
    // Whether the triangle points up or not;
    this.POINTS_UPWARD = pu;
    // The length of a triangle side to SVG elment size (ratio)
    this.SIDE_TO_MIN_BOUND = s_t_m_b;
  }
  get Hhalf() {
    return this.H/2;
  }
  get Whalf() {
    return this.W/2;
  }
}


// svg is the HTML Element.
// norm is either L1 or L2
function set_norm_to_use(svg,norm) {
  svg.triad_balance_state.NORM_TO_USE = norm;
}

// svg is the HTML Element.
// labels is an array of three strings
function set_labels(svg,labels) {
  svg.triad_balance_state.LABELS = labels;
  render_svg(svg,svg.triad_balance_state.FONT_SIZE_RATIO_TO_HEIGHT);
  rerender_marker(svg,svg.triad_balance_state.CUR_BALANCE);
}

// svg is the HTML Element.
// labels is an array of three strings
// Create an upward pointing triangle

function get_triangle(upward,svg,s_to_m_b = 7/10) {
  let W = svg.clientWidth;
  let H = svg.clientHeight;
  // we compute against whichever dimension is smaller
  var min = Math.min(W,H);

  const TRIANGLE_WIDTH = 1;
  const TRIANGLE_HEIGHT = Math.sqrt(3)/2;

  // This could be a parameter...
  let SIDE_LENGTH_PIXEL = min * s_to_m_b;

  const SIDE_LENGTH_HEIGHT = SIDE_LENGTH_PIXEL * TRIANGLE_HEIGHT;
  const BASE = -(1/3) * SIDE_LENGTH_HEIGHT;

  const UF = upward ? 1 : -1;

  let wtc_vector = [[-SIDE_LENGTH_PIXEL/2,BASE*UF],
                    [SIDE_LENGTH_PIXEL/2,BASE*UF],
                    [0,(BASE+SIDE_LENGTH_HEIGHT)*UF]];
// This is the "shield" formation (Scutum Fidei).
  return wtc_vector;
}

// remove the markers from the svg element (there may be only one or none.)
function clear_markers(svg) {
  var x = document.getElementsByClassName("triad-marker");
  for(var i = x.length -1; i >= 0; i--) {
    x[i].parentNode.removeChild(x[i]);
  }
}

// Graphics systems make y positve point downward, but our virtual triangle space has y upward
function vpy(y) { return -y; }
function vpx(x) { return x; }

// svg is the HTML Element
// bal_vec is the balance vector (a 3-vector in the unit attribute space.)
function rerender_marker(svg,bal_vec) {
  if (bal_vec) {
    // We have to find the current marker and remove it..
    clear_markers(svg);
    let tri_point = m.invertTriadBalance2to3(bal_vec,
                                           svg.triad_balance_state.TRIAD_WORLD_TRIANGLE,
                                           svg.triad_balance_state.NORM_TO_USE);
    let point = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
    point.setAttributeNS(null, 'cx', vpx(tri_point[0]));
    point.setAttributeNS(null, 'cy', vpy(tri_point[1]));
    // This value can be overridden in CSS
    point.setAttributeNS(null, 'r', 3);
    point.setAttributeNS(null,"class","triad-marker");
    point.ISMARKER = true;
    svg.appendChild(point);
  }
}

// svg is the HTML Element
// fs_ratio_to_height is the ratio of the font size of the labels to the total height
function render_svg(svg,fs_ratio_to_height) {
  let fs = svg.triad_balance_state.H * fs_ratio_to_height;

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
    point.x = vpx(svg.triad_balance_state.TRIAD_WORLD_TRIANGLE[i][0]);
    point.y = vpy(svg.triad_balance_state.TRIAD_WORLD_TRIANGLE[i][1]);
    polygon.points.appendItem(point);
  }

  function render_labels(svg,vertices,d_labels,fs) {
    // This is just what looks good to me, perhaps this should be
    // configurable or stylable.
    let pa = 4; // pixel_adjustment
    let vertical_adjustments = [fs+pa,fs+pa,-(fs/2 + pa)];
    for(let i = 0; i < 3; i++) {
//      console.log("DY",svg.triad_balance_state.POINTS_UPWARD ? vertical_adjustments[i] : -vertical_adjustments[i]+fs/2);
      append_text(svg,"triad-vertex-label-"+i,
                  "triad-vertices-labels",
                  vpx(svg.triad_balance_state.TRIAD_WORLD_TRIANGLE[i][0]),
                  vpy(svg.triad_balance_state.TRIAD_WORLD_TRIANGLE[i][1]),

                  // WARNING: This is somewhat arbitrary.
                  svg.triad_balance_state.POINTS_UPWARD ? vertical_adjustments[i] : -vertical_adjustments[i]+fs/2,
                  d_labels[i]
                 );
    }
  }

  render_labels(svg,svg.triad_balance_state.TRIAD_WORLD_TRIANGLE,
                svg.triad_balance_state.LABELS,
                fs);

  // this is the center of the triangle...
  let origin = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
  origin.setAttributeNS(null, 'cx', vpx(0));
  origin.setAttributeNS(null, 'cy', vpy(0));
  origin.setAttributeNS(null, 'r', 2);
  origin.setAttributeNS(null,"id","triangle_origin");
  svg.appendChild(origin);
}


// This is tricky because click events on an SVG
// depend on which object inside the SVG are hit.
// We don't really want to do that, we have
// created a global triangle space. A solution is to
// use screen coordinates, but only to compute a difference.

// evt is the mouse event
// fs is the font_size
// svg is the html SVG element
// labels is the set of strings
// click_callback is called when a click occurs with the position data
function clicked(evt,fs,svg,labels,click_callback) {
  var br = svg.getBoundingClientRect();
  var x = evt.clientX - br.left;
  var y = evt.clientY - br.top;
  // x and y are in screen coordinates of the
  // SVG ; we need to convert them
  // to the coordinates of our triangle.
  var xc = x + -svg.triad_balance_state.Whalf;

  let oriented_ydpc = (svg.triad_balance_state.POINTS_UPWARD ?
                       svg.triad_balance_state.Y_DISP_PER_CENT :
                       -svg.triad_balance_state.Y_DISP_PER_CENT)
  var yc = -(y + -svg.triad_balance_state.Hhalf
             + -svg.triad_balance_state.H * oriented_ydpc/100.0 );
  // Note, we balance and invert here to make sure we are inside the triangle!
  svg.triad_balance_state.CUR_TRIANGLE_COORDS = [xc,yc];

  svg.triad_balance_state.CUR_BALANCE =
    m.TriadBalance2to3(svg.triad_balance_state.CUR_TRIANGLE_COORDS,
                     svg.triad_balance_state.TRIAD_WORLD_TRIANGLE,
                     svg.triad_balance_state.NORM_TO_USE);

  var tri_point = m.invertTriadBalance2to3(
    svg.triad_balance_state.CUR_BALANCE,
    svg.triad_balance_state.TRIAD_WORLD_TRIANGLE,
    svg.triad_balance_state.NORM_TO_USE);

  rerender_marker(svg,svg.triad_balance_state.CUR_BALANCE);

  click_callback(svg.triad_balance_state.CUR_TRIANGLE_COORDS,
                 tri_point,
                 svg.triad_balance_state.CUR_BALANCE);
}
// svg is the HTML SVG element
// norm_to_use is either L1 or L2
// labels is an array of 3 strings
// click_callback is the callback that send the balance vector back on a click
// fs_ratio_to_height is the ration of the font_size to the height of svg
// s_to_m_b is the ratio of a side of the triangle to the minimum bound of the svg
// ydpc  is the Y displacement percent (downward) to make it look balanced
function initialize_triad_diagram(tbs) {
  tbs.SVG_ELT.triad_balance_state = tbs;

  // We need this as a separate function to handle resize events.
  function setSizeConstants(svg) {
    let W = svg.triad_balance_state.W = svg.clientWidth;
    let H = svg.triad_balance_state.H = svg.clientHeight;
    let Whalf = svg.triad_balance_state.Whalf;
    let Hhalf = svg.triad_balance_state.Hhalf;
    svg.triad_balance_state.TRIAD_WORLD_TRIANGLE =
      get_triangle(svg.triad_balance_state.POINTS_UPWARD,
                   svg,
                   svg.triad_balance_state.SIDE_TO_MIN_BOUND);
    // This is convenient, but makes it hard for the client to
    // use this svg for their own purposes, which is probably okay
    // for this use...
    // Here is were we can adjust the Hhalf value to move the triangle
    // down a bit, I think. However, this will like make us
    // have to make our inverstion functioon formal.
    let oriented_ydpc = (svg.triad_balance_state.POINTS_UPWARD ?
                         svg.triad_balance_state.Y_DISP_PER_CENT :
                         -svg.triad_balance_state.Y_DISP_PER_CENT)

    svg.setAttribute("viewBox",
                     `-${Whalf} -${Hhalf+H*oriented_ydpc/100.0} ${W} ${H}`);
    render_svg(tbs.SVG_ELT,svg.triad_balance_state.FONT_SIZE_RATIO_TO_HEIGHT);
    rerender_marker(tbs.SVG_ELT,svg.triad_balance_state.CUR_BALANCE);
  }

  setSizeConstants(tbs.SVG_ELT);

  tbs.SVG_ELT.addEventListener(
    "click",
    (evt) =>
      clicked(evt,tbs.SVG_ELT.triad_balance_state.FONT_SIZE_RATIO_TO_HEIGHT,
              tbs.SVG_ELT,tbs.LABELS,tbs.CLICK_CALLBACK)
  );

  window.addEventListener(
    "resize",
    (evt) => {
      setSizeConstants(tbs.SVG_ELT);
    });

}

module.exports = {
  vec: vec,
  m: m,
  initialize_triad_diagram: initialize_triad_diagram,
  get_triangle: get_triangle,
  TriadBalanceState: TriadBalanceState,
  set_norm_to_use: set_norm_to_use,
  set_labels: set_labels,
  TriadBalance2to3: m.TriadBalance2to3,
  invertTriadBalance2to3: m.invertTriadBalance2to3
};

console.log(vec);
