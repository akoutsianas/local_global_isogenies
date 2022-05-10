// AttachSpec("/cecm/home/akoutsia/Programming/Magma/ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);


// Initial equation for XD10 using Assaf's method and code
gens := [[4,0,0,3], [0,1,1,0], [2,0,0,2]];
N := 11;
H_N := sub<GL(2,Integers(N)) | gens>;
H := PSL2Subgroup(H_N);
M := ModularSymbols(H, 2, Rationals(), 0);
S := CuspidalSubspace(M);
XD10_not_reduced<[x]>, fs := ModularCurve(H);
AssignNames(~XD10_not_reduced, ["x0", "x1", "x2", "x3", "x4", "x5"]);


// Reduced equation for XD10
I := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, I);
isom, XD10toXD10_not_reduced := IsIsomorphic(XD10, XD10_not_reduced);


// Equation for modular curve X0(121) by Galbraith's thesis
I := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
-6*u^2 + 6*u*v - 3*v^2 - w^2 + 6*w*x - x^2 - 8*w*y + 10*x*y - 9*y^2 - 4*w*z + 10*x*z,
6*u^2 + 12*u*v + 12*v^2 - 17*w*x - 2*x^2 - 5*w*y + 4*x*y + 14*y^2 - 6*w*z - 7*x*z - 11*y*z,
-9*v^2 - 8*w^2 + w*x + 9*x^2 + 7*w*y - 10*x*y + y^2 - 7*w*z + 27*x*z - 11*y*z,
-6*u^2 - 12*u*v - 12*v^2 - 3*w^2 + 7*w*x - 6*x^2 + 11*w*y + 12*x*y + 10*y^2 + 4*w*z + 17*x*z - 11*y*z - 11*z^2];
X121 := Curve(P, I);






// Equation for X0+(121)
N := 121;
M := ModularForms(Gamma0(N), 2);
S := CuspidalSubspace(M);
basis := Basis(S);
A121 := AtkinLehnerOperator(S, 121);
I := ScalarMatrix(Dimension(S), 1);
B := A121 - I;
Eigenvecs := Basis(Nullspace(B));
prec := 100;
vf1 := Eigenvecs[1]*A121;
vf2 := Eigenvecs[2]*A121;
f1 := &+[vf1[i]*basis[i] : i in [1..6]];
f1 := qExpansion(f1, prec);
f2 := &+[vf2[i]*basis[i] : i in [1..6]];
f2 := qExpansion(f2, prec);
X := f1/f2;
Y := (Parent(X).1 * Derivative(X)) / f2;
Y2 := Y^2;
coef_bound := 50;
coefY := Vector(Rationals(), [Coefficient(Y2, i) : i in [-6..coef_bound]]);
coefX := [];
for k in [0..6] do
	Xk := X^k;
	coefXk := [Coefficient(Xk, i): i in [-6..coef_bound]];
	Append(~coefX, coefXk);
end for;
W := Matrix(coefX);
poly_coef := Solution(W, coefY);
_<x> := PolynomialRing(Rationals());
X0p := HyperellipticCurve(&+[x^i * poly_coef[i+1] : i in [0..6]]);
J0p := Jacobian(X0p);





// Equation for modular curve X0(121) by Assaf's method and code


