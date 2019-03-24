<!--
    Copyright 2019, Robert L. Read

    This file is part of TriBalance.

    TriBalance is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TriBalance is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TriBalance.  If not, see <https://www.gnu.org/licenses/>.
-->

<section>
<link rel="stylesheet" href="https://code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">
<script src="./js/three.js"></script>
<!-- Bootstrap CSS -->
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">
<script src="https://code.jquery.com/jquery-1.12.1.js"></script>
<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
<script src="https://code.jquery.com/jquery-1.12.4.js"></script>
<script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
<script src="./js/three.js"></script>
<!-- Prism for typesetting code -->
<link href="css/prism.css" rel="stylesheet" />
</section>    
<body>
<script src="./js/prism.js"></script>
<br>
  <div class="container">
    <h1 class="section-title h1">TriBalance</h1>
    <h2>First Three Dimensions</h2>
    <div>
      <button type="button" class="btn btn-primary" id="MBS">Mind/Body/Spirit</button>
      <button type="button" class="btn btn-secondary" id="HHH">Head/Heart/Hands</button>
      <button type="button" class="btn btn-success" id="PSM">Physical/Social/Mental</button>
      <button type="button" class="btn btn-danger" id="SEF">Strength/Endurance/Flexibility</button>
      <button type="button" class="btn btn-warning" id="FSH">Father/Son/Holy Spirit</button>
    </div>
    
    <fieldset>
    <legend>Balancing Norm to Use:</legend>
    <div style="display: flex; justify-content: space-around; font-size: large;">

    <span>
      <label for="l1norm">L1 Norm (Absolute Value):</label>
    <input type="radio" name="norm" id="l1norm" checked="checked">
    </span>
    <span>
    <label for="l2norm">L2 Norm (Sqrt of Squares):</label>
    <input type="radio" name="norm" id="l2norm">
    </span>    
    </div>
    </fieldset>


    <legend>Attributes Calculated from Triangle Position:</legend>
    <div style="display: flex; justify-content: space-around; font-size: x-large; font-weight: bold">
    <div>
      <label for="d0" id="d0l">Mind: </label>
    <span type="text" id="d0">0%</span>
    </div>
    <div>    
      <label for="d1" id="d1l">Body: </label>
    <span type="text" id="d1">0%</span>
    </div>    
    <div>    
      <label for="d2" id="d2l">Spirit: </label>
    <span type="text" id="d2">0%</span>
    </div>    
    </div>

</section>

<section>
<div style="display:flex; justify-content: center">
    <div id="container_to_have_global_coords_on_svg">
        <svg id="create_svg"  height="500px" width="500px" style="background: silver;" viewbox="-250 -250 500 500" font-size="24px"> </svg>
    </div>
</div>
</section>
<section>
  <h1>
    TriBalance Diagrams
</h1>
    <p>
    This is an explanation of what may be a new graphical user interface element---the "tri-balance" diagram. The purpose is to use a
triangle to enter three values such as "Mind/Body/Sprit" or "Strength/Endurance/Flexibility" that are presumed to in some sense
be "in balance" or where it is interesting if they are in balance. 
    </p>
    <p>
    This is an <a href="https://github.com/PubInv/TriBalance">open source demonstration project</a> of <a href="https://pubinv.github.io/PubInv/">Public Invention</a>. This grew out of the project <a href="https://github.com/PubInv/https://github.com/PubInv/SocialTetrahedrons">open source demonstration project</a>.
    Altough I searched a good bit, I may have missed a previous implementation if the tri-balance diagram has been implmented previously please inform me at &lt;read.robert@gmail.com&gt;.
  </p>
    <h2>Usage</h2>
Click near the yellow triangle and are red dot will appear within the triangle and will determine three values. These three values
will either sum to 100% or the square of the square roots will sum to 100%. 
The closer to one corner of the triangle you are, the more strongly you reflect that
  attribute. The center of the trinagle reflects all attributes equally--that is, "perfect balance."
  The weights of these attributes will be recorded with your assertion, and displayed (temporarily) below the
triangle.
    <h2>Purpose</h2>
    <p>
    Clicking on the triangle gives you a two-dimensional specification. In some cases, that may be all you need.
    </p>
    <p>
    However, if you want to execute some algorithm to determine, for example, which of two points has the higher "mind" value,
then calculating a 3-dimensional "balance vector" from the 2-dimensions makes more sense. It is only possible to specify three values
from a two-dimensional coordinate by assuming the values are in some sense balance and not completely independent of each other.
    </p>
    <p> Additionally, to be useful, the function that produces the "balance vector" from a point on the triangle must have an inverse.
    That is, there must exist code which when given a "balance vector" and the shape of a triangle as input, computes the coordinates
gives you back the same point you started with. This allows you to, for example, store balance vectors in a database and reconstuct
points on the GUI  at will. Moving back and forth between the model and the view is essential.
  <h2>L1 Norm vs. L2 Norm</h2>
  <p>
    To be able to write algorithms about the similarity or complimentarity of,
    for example, a person to a project in terms of Mind/Body/Spirit balance,
    we need a function to translate a position on the triangle into a mathematical
    object of some kind. We also must be able to invert this function, to
    obtain a position on the triangle from a mathematical object.
  </p>
  <p>
    The most obvious way to do this is to treat a given assertion as a vector
    with dimensions labeled, for example, Mind/Body/Spirit. Then a "balanced"
    vector maybe thought of an any vector of unit length. There are different
    ways to express the length of a vector. The L1 Norm is simply the sum
    of the absolute values, so by choosing that all dimensions sum to 1.0.
  </p>
  <p>
    However, it is also reasonable to use the L2 Norm, which is the "Euclidean distance." This is rather like saying every vector is a point on the surface
    of sphere of radius 1.0.  The components of the L2 Norm will sum to something
    usually greater than 1.0.
    </p>
    <h2> Basic Geometric Algorithm</h2>
    <img src="./Geometry of Tribalance.svg">
    <figure>
  <img src=""./Geometry of Tribalance.svg" alt="A diagram of the geometry of the algorithm with variables" style="width:100%">
  <figcaption>Fig.1 - Geometric quantities represented by variables in the algorithm TriBalance2to3</figcaption>
</figure>
<p>
    Please refer to the diagram above when considering Algorithm TriBalance2to3 below.
In the algorithm below, I use the <a href="https://threejs.org/">THREE.js</a> library, which provides a few useful vector operations (.add,.multiplyScalar,.clampLength). However, these could be easily replaced with a lighter library or those functions could be coded by hand.
    We should all be grateful to the contributors to THREE.js, but it is coded in an
imperative, rather than a functional style, which creates extra lines of code and
a bit of awkwardness in the code below. (I am fan of purely functional programming whenever possible.)
    Additionally, I use
a slight modification of a routine "GetRayToLineSegmentIntersection" <a href="https://stackoverflow.com/questions/14307158/how-do-you-check-for-intersection-between-a-line-segment-and-a-line-ray-emanatin">provided</a>
    on Stack Overflow by <a href="https://stackoverflow.com/users/71689/ezolotko">ezolotko</a>is called.
    </p>
    <p>
    The basic geometric algorithm to compute the three values from a point on the triangle is straightforward. The basic insight is that any point specifies an amount of "balance",
proportionaly to how close it is to the center of the triangle,
represented in the vector "bal" and an amount of imbalance, depending on where a ray
intersects the edge, represented by the variable
"imb". The final result is a linear combination of these two unit-length vectors, which is guaranteed by
linearity to have a unit-length norm. Note that the function is parametrized by the
norming function, we can use either the L1 or L2 norm (or any other.)
</p>
    <p> Before we begin, let us note that special case of the input point being precisely
the origin or a vertex is in fact handled correctly, although this is perhaps not obvious;
this is tested in our unit test suite.
    </p>
    <p> The TriBalance2to3 algorithm takes as input a triangle point, three vertices, and
a selection of which norm to use (L1 and L2 norms are provided).</p>
    <pre><code class="language-js">
// tp is the a point in the 2-dimensional triangle space
// wtc are the three vertices of an eqilateral triangle whose centroid is the origin
// LXnorm_and_length is a pair of functions to to normalize a vector and compute the length
// return the corresponding 3-vector in the attribute space
    function TriBalance2to3(tp,wtc,LXnorm_and_length = L2) {
        </code></pre>
            <p> For clarity, we provide a name for the normalization function.</p>
    <pre><code class="language-js">
    let LXnormalize = LXnorm_and_length[0];            
        </code></pre>        
    <p>
            If the point is the triangle exact center of the triangle, we cannot
        construct a ray to an edge, and we return a perfectly balanced vector.
             </p>
            <pre><code class="language-js">
    if (near(tp.lengthSq(),0,1e-5)) {
        return LXnormalize(new THREE.Vector3(1,1,1));
    }
        </code></pre>        
            
            <p> Find the point "point_on_edge" where the line from the center of the triangle through the point intersects a triangle edge.
            Such a point must exists, and is on one edge or the intersection of two edges. Imagine a ray from the origin through point p to point E.
            The point E is along the ray OP and length is at most the distance from the origin to the edge along this ray (to keep the point inside the triangle.)
        </p>
            <pre><code class="language-js">
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
        </code></pre>        
            
    <p>
    Compute the length from P to E, and the total length from O to P.
            </p>
            <pre><code class="language-js">
    let total_distance_to_edge = origin.distanceTo(point_on_edge);
    // If the point is outside the triangle, we clamp (truncate if needed)
    // it's length so that it is precisely on the edge.
    p.clampLength(0,total_distance_to_edge);

    let distance_to_p_o_e = p.distanceTo(point_on_edge);            
        </code></pre>        
            <p> Determine the ratio (p_to_edge) of PE to OP. If this is close to zero, we are close to perfectly balanced; if it is close to 1, we are close to having zero for one attribute.
              Compute the a balanced vector BAL  with all values the same, scaled by the ratio to the p_to_edge. 
            </p>
    <pre><code class="language-js">        
    var ratio_p_to_edge =  distance_to_p_o_e/total_distance_to_edge;
    
    let bal = LXnormalize(new THREE.Vector3(1,1,1));
    bal.multiplyScalar(ratio_p_to_edge);
        </code></pre>
            <p>
            Now to compute the "imbalance", we linearly divide the edge
        into the distance to the two vertices. Possibly one of these lengths will be zero.
            </p>
        <pre><code class="language-js">
    // Now the remainder of the contribution
    // to the unit vector should come from the two
    // points on the edge, in linear proportion.
    // These coordinates are fe_idx and (fe_idx+1) % 3.
    var d1 = wtc[fe_idx].distanceTo(point_on_edge);
    var d2 = wtc[(fe_idx+1) % 3].distanceTo(point_on_edge);
        </code></pre>
            <p>We assign these two distances to to the corresponding attributes
            in a 3-vector. The remaining one is zero. We then normalize the vector,
        driving it&#39;s length to one, and then scale it by an amount chosen such
that by linearity the sum with the balance vector will be one: (1 - ratio_p_to_edge).
        <pre><code class="language-js">    
    let vs = [0,0,0];
    vs[fe_idx] = d2;
    vs[(fe_idx+1) % 3] = d1;
    
    let imb = LXnormalize(new THREE.Vector3(vs[0],vs[1],vs[2]));
    imb.multiplyScalar(1 - ratio_p_to_edge);
        </code></pre>            
    <p>
Now the vector bal and the vector imb added together are the return value.
            </p>
    <pre><code class="language-js">        
    // now construct a balanced vector proportional
    // to the length from the edge to the p towards the axis
    // so that this be a unit vector if p is the origin.
    return new THREE.Vector3().add(imb).add(bal);
}
        </code></pre>        
    <p>
    Return the sum of BAL and IMB.
     </p>
    <h2> The Inverse Function</h2>
    <p> The function TriBalance2to3 calculates a 3-vector from the clicked point on the triangle; the inverse function invertTriBlance2to3 gives you back the clicked point when given the 3-vector. The utility of producing a 3-vector to represent the 3 attributes in balance depends entirely on having such an inverse function; no alternative TriBalance2to3 function which does not have an inverse can be useful.
    </p>
    <p>
    This fundamental insight of this algorithm is that since the ray trough the triangle point
touches at most two edges, there is always a component of the attribute vector that
is purely contributed by the "balance" component, and it is always a minimum of the
values in the vector. We can therefore find this minimum, construct a balance vector
with those values, subtract it from the input, and have a vector representing pure
imbalance which has at least one zero in it. The two non-zero values represent a linear
interpolation along the edge (or point, if two edges) that the ray strikes. Since
we constucted the attribut value as a ratio of lengths between the vertices, we can
use this fact to reconstruct, via linear interpretation, a point on the edge.
    A second linear interpretation from the point based on the imbalance length gives
us the point in two space.
    </p>
    <pre><code class="language-js">
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
    </code></pre>
        <p>
        Now it is the case that by construction, imb has has at least one zero
    (whichever attributes were minimal), and bal has all attributes equal.
        </p>
    <pre><code class="language-js">    
    // Now that we have balance, we need to compute it's length,
    // which is dependent on the norm we chose!

    let imb_r = length(imb);
    let bal_r = length(bal);
    console.assert(Math.abs((bal_r+imb_r) - 1) <   1e-5);
    </code></pre>
        <p> The ratio computed below is a ratio of values in attribute space;
    it is not obious therefore that this can be used to perform a linear
    interpreation in the <em>triangle</em> coordinate space. However,
    we intensionally constructed the attribute vector in proportion to
    the distance in the triangle space, so this works no matter which norm
    we use in the attribute space.
        </p>
    <pre><code class="language-js">        
    // Now we have the ratios. We need to determine the direction.
    // this is a function of the imbalance vector. We could determine
    // which side we are on, and then compute our position along that
    // to determine a point on the triangle, and then multiply by the imb_r
    // to obtain the actual point.
    // At least one value of imb will be zero.
    var from, to, ratio;
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
    </code></pre>
        <p>The THREE.js library provides is a linear interpolation betwen
    two vectors, named "lerpVectors", out of the box. We use this once
    to interpolate along the edge, and then once to interpolate from th
    origin towared this point.
    <pre><code class="language-js">        
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
        </code></pre>        

    <h2>Alternative Approaches</h2>
    <p>
    I am not a user interface expert; nonetheless I offer these comments. The obvious alternative is to use three sliders. However, this at a minimum requires the user to make two clicks. It also does not visually represent the concept of interdependence.
    </p>
    <p>
    A possible alternative which may be more attractive and in a since more elegant is to use a circle, with a the same triangle inscribed. This would be similar, but make better use of the screen space corresponding to the bounding box. The math would be slightly different, and probably simpler, than what is presented here.
    </p>
    <p>
    Finally, it might be possible to define a GUI element that allows you to select from not 3, but four or more dimenstions, such as the elements of activity: Earth/Air/Water/Fire. However, this would necessarily impose an even further constraint than "balance" on what is possible. For example, if any number of elements can be
arranged radially in a natural way so that being close to one completely excludes being close to those that are 2 or more elments away, this could be a natural system.
    </p>
<h2> How to Best Reuse </h2>
This project is contributing javascript and math. To make an easily reusable GUI component requires configuration and design skill. I would love to have somebody take this 
code and make a <a href="https://d3js.org/"> d3</a> or <a href="https://reactjs.org/">React</a> component out of it, so that it could be enjoyed by others as a plug-in as easy using another d3 GUI element, for example. To perform vector operations, these algorithms use <a href="https://threejs.org/">THREE.js</a>. One could selectively replace the small number of operations used from that library to make it a pure javascript project with a light footprint (a few hundred lines of code at most for the math part.)
    
  <h2>License</h2>
  You are free to reuse this software under the terms
  of the GNU General Public License.
</section>
</div>

</div>

    <!-- Optional JavaScript -->
    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js" integrity="sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy" crossorigin="anonymous"></script>




<!-- This is brings in the math and the norm defintions we need -->
    <script src="./TriBalance.js"></script>

    <script>

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
     let wtc = [[-SIDE_LENGTH_PIXEL/2,BASE],
                             [SIDE_LENGTH_PIXEL/2,BASE],
                [0,BASE+SIDE_LENGTH_HEIGHT]];
    let wtc_vector = [new THREE.Vector2(wtc[0][0],wtc[0][1]),
               new THREE.Vector2(wtc[1][0],wtc[1][1]),
               new THREE.Vector2(wtc[2][0],wtc[2][1])];
    return wtc_vector;
}


// These function convert my abstract coordinates
// to canvas coordinates.
function vph(h,y) { return (-y); }
function vpw(w,x) { return (x); }

var CUR_POINT;

var CUR_TRIANGLE_COORDS;

function append_text(svg,x,y,text) {
    var newText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    newText.setAttributeNS(null,"x",x);      
    newText.setAttributeNS(null,"y",y);
    newText.appendChild(document.createTextNode(text));
    svg.appendChild(newText);
}


function render_svg() {
    var svg = $("#create_svg")[0];
    $("#create_svg").empty();
    var polygon = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
    svg.appendChild(polygon);

    let wtc = WORLD_TRIANGLE_COORDS;

    var array = [ [ vpw(W,wtc[0].x),vph(H,wtc[0].y) ], 
                  [ vpw(W,wtc[1].x),vph(H,wtc[1].y) ],
                  [ vpw(W,wtc[2].x),vph(H,wtc[2].y) ] ];
    
    for (value of array) {
        var point = svg.createSVGPoint();
        point.x = value[0];
        point.y = value[1];
        polygon.points.appendItem(point);
    }
    polygon.style.fill='lemonchiffon';

    

    // These are ugly, they should really be computed from the text.
    // In fact, since this does not change, the whole thing could go into
    // HTML and css more profitably.
    
    append_text(svg,array[2][0]-20,array[2][1]-5,
                THREE_DIMENSIONS[CURRENT_D][2]);
    append_text(svg,array[0][0]-40,array[0][1]+20,
                THREE_DIMENSIONS[CURRENT_D][0]);
    append_text(svg,array[1][0],array[1][1]+20,
                THREE_DIMENSIONS[CURRENT_D][1]
               );
    
    var circle = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
    circle.setAttributeNS(null, 'cx', vpw(W,0));
    circle.setAttributeNS(null, 'cy', vph(H,0));
    circle.setAttributeNS(null, 'r', 2);
    circle.setAttributeNS(null, 'style', 'fill: black; stroke: black; stroke-width: 1px;' );
    svg.appendChild(circle);
    
    function add_triangle(tri,c) {
        if (tri) {
            var polygon = document.createElementNS("http://www.w3.org/2000/svg", 'circle');
            polygon.setAttributeNS(null, 'cx', vpw(W,tri.x));
            polygon.setAttributeNS(null, 'cy', vph(H,tri.y));
            polygon.setAttributeNS(null, 'r', 4);
            polygon.setAttributeNS(null, 'style', 'fill: '+c+'; stroke: '+c+'; stroke-width: 1px;' );
            svg.appendChild(polygon);
            
        }
    }
    add_triangle(CUR_POINT,"red");
}


function main() {
    render_svg();

    function set_dimension_labels(cur) {
        $("#d0l").text(THREE_DIMENSIONS[cur][0] + ":");
        $("#d1l").text(THREE_DIMENSIONS[cur][1] + ":");
        $("#d2l").text(THREE_DIMENSIONS[cur][2] + ":");
    }
    function set_and_render(cur) {
        CURRENT_D = cur;
        set_dimension_labels(CURRENT_D);
        render_svg();        
    }
    
    $("#MBS").click(() => { set_and_render("MBS");});
    $("#HHH").click(() => { set_and_render("HHH");});
    $("#PSM").click(() => { set_and_render("PSM");});
    $("#SEF").click(() => { set_and_render("SEF");});
    $("#FSH").click(() => { set_and_render("FSH");});
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
    let bal = TriBalance2to3(atp,WORLD_TRIANGLE_COORDS,norm_to_use);

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
    var yc = -(y + -H/2) ;
    var xc = x + -W/2;
    var triangle_coords = new THREE.Vector2(xc,yc);

    // Note, we could balance and invert here to make sure we are inside the trianble!
    CUR_TRIANGLE_COORDS = triangle_coords;
    var bal = setBalance(triangle_coords);
    var vec = new THREE.Vector3(bal[0],bal[1],bal[2]);
    
    var norm_to_use = (getRadioValue("norm") == 0 ? L1 :L2);
    var triangle_coords_inside_triangle = invertTriBalance2to3(vec,WORLD_TRIANGLE_COORDS,norm_to_use);

    CUR_POINT = triangle_coords_inside_triangle;
    render_svg();
}         

$("#container_to_have_global_coords_on_svg").click(clicked);

$( document ).ready(function() {
    
    WORLD_TRIANGLE_COORDS = get_world_triangle();
    
    testAllTriBalance(WORLD_TRIANGLE_COORDS);
    
    main();
});

</script>
</body>
</html>
