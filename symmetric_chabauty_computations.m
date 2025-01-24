//////// Attach file "symmetric_chabauty.m" ////////


P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);


////// Model of XD10 and the Atkin-Lehner invoslution //////

ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);
// w11 := iso<XD10 -> XD10 | [-P.1, -P.2, -P.3, P.4, P.5, -P.6],[-P.1, -P.2, -P.3, P.4, P.5, -P.6]>;


// G := AutomorphismGroup(XD10,[w11]);
// C, XD10ToC := CurveQuotient(G);
// CQQ := [C![1, -3, 1], C![1, -1, 1], C![0, -1, 1], C![0, 0, 1], C![1, -1, 0], C![1, 0 , 0]];
// JC := Jacobian(C);
// GJC, GJCToJC := MordellWeilGroup(JC);
// JCToGJC := Inverse(GJCToJC);


////// Finite index subgroup of JD10(Q) //////

R := CoordinateRing(P);
D1a := Divisor(XD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
Dinf := Divisor(XD10, Ideal([R.1*R.6 + 1/3*R.6^2, R.2*R.6, R.3*R.6 + 1/3*R.6^2, R.4^2 - 22/9*R.6^2, R.5^2 - 22/9*R.6^2]));
D1 := D1a - Dinf;

D2a := Divisor(XD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 + 11*R.6^2, R.5^2 + 11*R.6^2]));
D2 := D2a - Dinf;



// The set the set L of points on XD10^(2)(Q)

P1 := Divisor(XD10, Ideal([4*R.1*R.6 + 3*R.6^2, 4*R.2*R.6 - R.6^2, R.3*R.6, 4*R.4^2 - 77*R.6^2, R.5*R.6]));
P2 := Divisor(XD10, Ideal([4*R.1*R.6 - 3*R.6^2, 4*R.2*R.6 + 5*R.6^2, R.3*R.6, 4*R.4^2 - 77*R.6^2, R.5*R.6]));
P3 := Divisor(XD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 + 11*R.6^2, R.5^2 + 11*R.6^2]));
P4 := Divisor(XD10, Ideal([5*R.1*R.6 + 2*R.6^2, 5*R.2*R.6 - 2*R.6^2, 5*R.3*R.6 - R.6^2, 25*R.4^2 - 209*R.6^2, 25*R.5^2 - 209*R.6^2]));
P5 := Divisor(XD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
P6 := Divisor(XD10, Ideal([3*R.1*R.6 + R.6^2, R.2*R.6, 3*R.3*R.6 + R.6^2, 9*R.4^2 - 22*R.6^2, 9*R.5^2 - 22*R.6^2]));
S := [P1, P2, P3, P4, P5, P6];
Gens := SetToSequence({D1a, D2a});



//// Apply Symmetric Chabauty and Mordel-Weil Sieve ////
primes := PointsLieOnSingleResidueClass(XD10, S, [5, 7, 13, 17]);
print primes;
MordellWeilSieve(XD10, Gens, [5, 0], S, primes, 10);