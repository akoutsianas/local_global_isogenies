// AttachSpec("/cecm/home/akoutsia/Programming/Magma/ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);


// Initial equation for XD10 using Assaf's method and code
gens := [[4,0,0,3], [0,1,1,0], [2,0,0,2]];
N := 11;
H_N := sub<GL(2,Integers(N)) | gens>;
H := PSL2Subgroup(H_N);
M := ModularSymbols(H, 2, Rationals(), 0);
S := CuspidalSubspace(M);
XD10nred<[x]>, basisD10 := ModularCurve(H);
//AssignNames(~XD10_not_reduced, ["x", "y", "z", "u", "v", "w"]);


// Reduced equation for XD10
ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);
isom, XD10ToXD10nred := IsIsomorphic(XD10, XD10nred);
isom, XD10nredToXD10 := IsIsomorphic(XD10nred, XD10);
// assert isom;



// Equations for the modular curve X0(121)

//P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);
N := 121;
M := ModularForms(Gamma0(N), 2);
S := CuspidalSubspace(M);
basis121 := Basis(S);



// Equation for modular curve X0(121) by Galbraith's thesis
I121Gal := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
-6*u^2 + 6*u*v - 3*v^2 - w^2 + 6*w*x - x^2 - 8*w*y + 10*x*y - 9*y^2 - 4*w*z + 10*x*z,
6*u^2 + 12*u*v + 12*v^2 - 17*w*x - 2*x^2 - 5*w*y + 4*x*y + 14*y^2 - 6*w*z - 7*x*z - 11*y*z,
-9*v^2 - 8*w^2 + w*x + 9*x^2 + 7*w*y - 10*x*y + y^2 - 7*w*z + 27*x*z - 11*y*z,
-6*u^2 - 12*u*v - 12*v^2 - 3*w^2 + 7*w*x - 6*x^2 + 11*w*y + 12*x*y + 10*y^2 + 4*w*z + 17*x*z - 11*y*z - 11*z^2];
X121Gal := Curve(P, I121Gal);



// Equation for X0(121) using Basis function
prec := 100;
basis_sq := [];
vars_sq := [];
for i in [1..6] do
	for j in [1..6] do
		fij := basis121[i] * basis121[j];
		varij := P.i * P.j;
		if fij notin basis_sq then
			Append(~basis_sq, fij);
			Append(~vars_sq, varij);
		end if;
	end for;
end for;


W := [];
for f in basis_sq do
	qcoefs := [Coefficient(f, i) : i in [0..prec]];
	Append(~W, qcoefs);
end for;
W := Matrix(W);
I := [];
for v in Basis(Nullspace(W)) do
	poly := &+[v[i] * vars_sq[i] : i in [1..#vars_sq]];
	Append(~I, poly);
end for;
X121 := Curve(P, I);
// isom, X121toX121Gal := IsIsomorphic(X121, X121Gal);
// assert isom;



// Equation for X0+(121)
A121 := AtkinLehnerOperator(S, 121);
Id := ScalarMatrix(Dimension(S), 1);
B := A121 - Id;
Eigenvecs := Basis(Nullspace(B));
prec := 100;
vf1 := Eigenvecs[1]*A121;
vf2 := Eigenvecs[2]*A121;
f1 := &+[vf1[i]*basis121[i] : i in [1..6]];
f1 := qExpansion(f1, prec);
f2 := &+[vf2[i]*basis121[i] : i in [1..6]];
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
_<t> := PolynomialRing(Rationals());
f := &+[t^i * poly_coef[i+1] : i in [0..6]];
X0p := HyperellipticCurve(f);
J0p := Jacobian(X0p);


// The map X0(121) -> X0+(121)
Xdenom := &+[vf2[i] * P.i : i in [1..6]];
R := Parent(Xdenom);
RI := R/ideal<R | I>;
Ynum := Y * f2^3;
basis_cub := [];
vars_cub := [];
for i in [1..6] do
        for j in [1..6] do
		for k in [1..6] do
	                fijk := basis121[i] * basis121[j] * basis121[k];
        	        varijk := P.i * P.j * P.k;
                	if fijk notin basis_cub then
                        	Append(~basis_cub, fijk);
                        	Append(~vars_cub, varijk);
                	end if;
		end for;
        end for;
end for;

W := [];
for f in basis_cub do
        qcoefs := [Coefficient(f, i) : i in [0..coef_bound]];
        Append(~W, qcoefs);
end for;
W := Matrix(W);
coefYnum := Vector(Integers(), [Coefficient(Ynum, i) : i in [0..coef_bound]]);
Ynum_coef := Solution(W, coefYnum);
Ynum_var := &+[Ynum_coef[i] * vars_cub[i] : i in [1..Degree(Ynum_coef)]];
Xmap := &+[vf1[i] * P.i : i in [1..6]] / &+[vf2[i] * P.i : i in [1..6]];
Ymap := Ynum_var / (&+[vf2[i] * P.i : i in [1..6]])^3;
Ymap := R!(RI!Numerator(Ymap)) / Denominator(Ymap);
X121ToX0p := map<X121 -> X0p | [Xmap, Ymap, 1]>;


// Map XD10 -> X121

PK<x,y,z,u,v,w> := ProjectiveSpace(K,5);
L := BaseRing(Parent(basisD10[1]));
K<r11> := QuadraticField(-11);
bol, KtoL := IsSubfield(K, L);
LtoK := Inverse(KtoL);

Mbasis := Matrix([[Coefficient(fi, i): i in [1..99]]: fi in basis121]);
Mbasis := ChangeRing(Mbasis, L);
MD10nred := Matrix([[Coefficient(fi,i): i in [1..99]] : fi in basisD10]);
A := Solution(MD10nred, Mbasis);
A :=  ChangeRing(A, K, LtoK);
XD10nredK := ChangeRing(XD10nred, K);
X121K := ChangeRing(X121, K);
XD10K := ChangeRing(XD10, K);
X121GalK := ChangeRing(X121Gal, K);
XD10nredToX121 := map<XD10nredK -> X121K | [&+[ri[i]*PK.i: i in [1..6]] : ri in Rows(A)]>;
polys := [Evaluate(fi, [PK.1, PK.2, PK.3, PK.4, PK.5, PK.6]) : fi in DefiningPolynomials(XD10ToXD10nred)];
XD10ToXD10nred := map<XD10K -> XD10nredK | polys>;
// XD10ToX121 := XD10ToXD10nred * XD10nredToX121;
XD10ToX121 := map<XD10K -> X121K | DefiningPolynomials(XD10ToXD10nred * XD10nredToX121)>;


// Map X121 -> XD10
B := Solution(Mbasis, MD10nred);
B :=  ChangeRing(B, K, LtoK);
X121ToXD10nred := map<X121K -> XD10nredK | [&+[ri[i]*PK.i: i in [1..6]] : ri in Rows(B)]>;
polys := [Evaluate(fi, [PK.1, PK.2, PK.3, PK.4, PK.5, PK.6]) : fi in DefiningPolynomials(XD10nredToXD10)];
XD10nredToXD10 := map<XD10nredK -> XD10K | polys>;
//X121ToXD10 := X121ToXD10nred * XD10nredToXD10;
X121ToXD10 := map<X121K -> XD10K | DefiningPolynomials(X121ToXD10nred * XD10nredToXD10)>;


// Find the basis of cusp forms of X121Gal
R<x,y,z,u,v,w> := PolynomialRing(Rationals(), 6);
RK<x,y,z,u,v,w> := PolynomialRing(K, 6);
XD10ToX121Gal := map<XD10K -> X121GalK | [r11*PK.1, r11*PK.2, r11*PK.3, PK.4, PK.5, r11*PK.6]>;
X121GalToXD10 := map<X121GalK -> XD10K | [PK.1, PK.2, PK.3, r11*PK.4, r11*PK.5, PK.6]>;
X121GalToX121 := X121GalToXD10 * XD10ToX121;
X121GalToX121 := map<X121Gal -> X121 | [R!(f/r11) : f in DefiningPolynomials(X121GalToX121)]>;
X121ToX121Gal := X121ToXD10 * XD10ToX121Gal;
polys := [R!f : f in DefiningPolynomials(X121ToX121Gal)];
X121ToX121Gal := map< X121 -> X121Gal | polys>;
LS := PowerSeriesRing(Rationals());
basis121exp := [(LS!qExpansion(f, prec)) : f in basis121];
basis121Galexp := [(LS!Evaluate(f, basis121exp)) : f in polys];
Mbasis := Matrix(Rationals(),[[Coefficient(fi, i): i in [1..90]]: fi in basis121exp]);
MGalbasis := Matrix(Rationals(), [[Coefficient(fi, i): i in [1..90]]: fi in basis121Galexp]);
// Solution(Mbasis, MGalbasis);



// Equation for X121Gal^+
// Compute action of Atkin-Lehner involution on X121GAl
vAGalexp := [&+[r[i] * basis121exp[i]: i in [1..6]] : r in Rows(A121)];
A121Galexp := [Evaluate(f, vAGalexp) : f in polys];
A121Galexp := Matrix(Rationals(), [[Coefficient(fi, i) : i in [1..90]] : fi in A121Galexp]); 
A121Gal := Solution(MGalbasis, A121Galexp);
B := A121Gal - Id;
Eigenvecs :=  Basis(Nullspace(B));
phi := iso<X121Gal -> X121Gal | [-P.1, -P.2, -P.3, P.4, P.5, -P.6],[-P.1, -P.2, -P.3, P.4, P.5, -P.6]>;
G := AutomorphismGroup(X121Gal,[phi]);
C, X121GalToC := CurveQuotient(G);



polys := [RK!f : f in DefiningPolynomials(X121ToXD10)];
LSK :=  PowerSeriesRing(K);
basis121exp := [(LSK!qExpansion(f, prec)) : f in basis121];
basisD10exp := [(LSK!Evaluate(f, basis121exp)) : f in polys];
Mbasis := Matrix(K,[[Coefficient(fi, i): i in [1..90]]: fi in basis121exp]);
MD10basis := Matrix(K, [[Coefficient(fi, i): i in [1..90]]: fi in basisD10exp]);
Solution(Mbasis, MD10basis);




for i1,i2,i3,i4,i5,i6 in [1,r11] do
	I121new := [Evaluate(f, [i1*PK.1, i2*PK.2, i3*PK.3, i4*PK.4, i5*PK.5, i6*PK.6]) : f in I121];
	for f in I121new do
		try
			w := R!f;
		catch e
			//print "No element in R";
			break;
		end try;
		print "Success";
	end for; 
end for;
