// Attach Assaf's code. Attach spec AttachSpec("path_to_folder/ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");


///////////// Models Computations /////////////

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);


///// Equation for XD10 /////
ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);




///// Initial equation for XD10 using Assaf's method /////
gens := [[4,0,0,3], [0,1,1,0], [2,0,0,2]];
N := 11;
H_N := sub<GL(2, Integers(N)) | gens>;
H := PSL2Subgroup(H_N);
M := ModularSymbols(H, 2, Rationals(), 0);
S := CuspidalSubspace(M);
XD10nred, basisD10 := ModularCurve(H);
AssignNames(~XD10nred, ["x", "y", "z", "u", "v", "w"]);
isom, XD10ToXD10nred := IsIsomorphic(XD10, XD10nred);
assert isom;
// isom, XD10nredToXD10 := IsIsomorphic(XD10nred, XD10);



///// Equation for modular curve X0(121) by Galbraith's thesis /////
I121 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
-6*u^2 + 6*u*v - 3*v^2 - w^2 + 6*w*x - x^2 - 8*w*y + 10*x*y - 9*y^2 - 4*w*z + 10*x*z,
6*u^2 + 12*u*v + 12*v^2 - 17*w*x - 2*x^2 - 5*w*y + 4*x*y + 14*y^2 - 6*w*z - 7*x*z - 11*y*z,
-9*v^2 - 8*w^2 + w*x + 9*x^2 + 7*w*y - 10*x*y + y^2 - 7*w*z + 27*x*z - 11*y*z,
-6*u^2 - 12*u*v - 12*v^2 - 3*w^2 + 7*w*x - 6*x^2 + 11*w*y + 12*x*y + 10*y^2 + 4*w*z + 17*x*z - 11*y*z - 11*z^2];
X121 := Curve(P, I121);


//// Compute the Atkin-Lehner involution ////

// Atkin-Lehner involution of X0(121) is defined over Q.
// The automorphisms of X121 defined over Q is a group of order 2, hence the map [x, y, z, u, v, w] -> [-x, -y, -z, u, v, -w] is 
// the Atkin-Lehner involution.

AutX121 := AutomorphismGroup(X121);
assert Order(AutX121) eq 2;


///// Equation for X0+(121) and XpD10 /////

w11 := iso< X121 -> X121 | [-x, -y, -z, u, v, -w], [-x, -y, -z, u, v, -w]>;
Xp121, X121ToXp121 := CurveQuotient(AutomorphismGroup(X121, [w11]));

w11 := iso< XD10 -> XD10 | [-x, -y, -z, u, v, -w], [-x, -y, -z, u, v, -w]>;
XpD10, XD10ToXpD10 := CurveQuotient(AutomorphismGroup(XD10, [w11]));

assert XpD10 eq Xp121;

