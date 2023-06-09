// bubble radius:
R = 0.25;
Rn = 10*R;
Rk = 25*R;

// axisymmetri domain:
zmin = -30.0 * R;
zmax = 30.0 * R;
rmin = 0.;
rmax = 30. * R;

//refinement value:
bc = .01;
nc = 0.1;
kc = 0.5;
dc = 1;

// domain boundary points:
Point(1) = {zmin, rmin, 0, dc};
Point(2) = {0, rmin, 0, bc};
Point(3) = {zmax, rmin, 0, dc};
Point(4) = {zmax, rmax, 0, dc};
Point(5) = {0, rmax, 0, dc};
Point(6) = {zmin, rmax, 0, dc};

// bubble arc points:
Point(7) = {-R, 0, 0, bc};  
Point(8) = {0, R, 0, bc};   
Point(9) = {R, 0, 0, bc};   
Point(10) = {-R*0.70710678118, R*0.70710678118, 0, bc};
Point(11) = {R*0.70710678118, R*0.70710678118, 0, bc};


// bubble neighborhood arc points:
Point(12) = {-Rn, 0, 0, nc};  
Point(13) = {0, Rn, 0, nc};   
Point(14) = {Rn, 0, 0, nc};   
Point(15) = {-Rn*0.70710678118, Rn*0.70710678118, 0, nc};
Point(16) = {Rn*0.70710678118, Rn*0.70710678118, 0, nc};


// Kirchhoff arc points:
Point(17) = {-Rk, 0, 0, kc};  
Point(18) = {0, Rk, 0, kc};   
Point(19) = {Rk, 0, 0, kc};   
Point(20) = {-Rk*0.70710678118, Rk*0.70710678118, 0, kc};
Point(21) = {Rk*0.70710678118, Rk*0.70710678118, 0, kc};

 
//line connecting domain boundary points:
Line(1) = {1, 17};
Line(2) = {17, 12};
Line(3) = {12, 7};
Line(4) = {7, 2};
Line(5) = {2, 9};
Line(6) = {9, 14};
Line(7) = {14, 19};
Line(8) = {19, 3};
Line(9) = {3, 4};
Line(10) = {4, 5};
Line(11) = {5, 6};
Line(12) = {6, 1};

Curve Loop(1) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(1) = {1};
Recombine Surface {1};
