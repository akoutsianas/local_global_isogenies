/*

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);

ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);
w11 := iso<XD10 -> XD10 | [-P.1, -P.2, -P.3, P.4, P.5, -P.6],[-P.1, -P.2, -P.3, P.4, P.5, -P.6]>;
G := AutomorphismGroup(XD10,[phi]);
C, XD10ToC := CurveQuotient(G);
CQQ := [C![1, -3, 1], C![1, -1, 1], C![0, -1, 1], C![0, 0, 1], C![1, -1, 0], C![1, 0 , 0]];
JC := Jacobian(C);
GJC, GJCToJC := MordellWeilGroup(JC);
JCToGJC := Inverse(GJCToJC);


// Finite index subgroup of JD10(Q)

R := CoordinateRing(P);
D1a := Divisor(XD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
Dinf := Divisor(XD10, Ideal([R.1*R.6 + 1/3*R.6^2, R.2*R.6, R.3*R.6 + 1/3*R.6^2, R.4^2 - 22/9*R.6^2, R.5^2 - 22/9*R.6^2]));
D1 := D1a - Dinf;

D2a := Divisor(XD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 + 11*R.6^2, R.5^2 + 11*R.6^2]));
D2 := D2a - Dinf;


// A the set S of points on XD10^(2)(Q)

P1 := Divisor(XD10, Ideal([4*R.1*R.6 + 3*R.6^2, 4*R.2*R.6 - R.6^2, R.3*R.6, 4*R.4^2 - 77*R.6^2, R.5*R.6]));
P2 := Divisor(XD10, Ideal([4*R.1*R.6 - 3*R.6^2, 4*R.2*R.6 + 5*R.6^2, R.3*R.6, 4*R.4^2 - 77*R.6^2, R.5*R.6]));
P3 := Divisor(XD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 + 11*R.6^2, R.5^2 + 11*R.6^2]));
P4 := Divisor(XD10, Ideal([5*R.1*R.6 + 2*R.6^2, 5*R.2*R.6 - 2*R.6^2, 5*R.3*R.6 - R.6^2, 25*R.4^2 - 209*R.6^2, 25*R.5^2 - 209*R.6^2]));
P5 := Divisor(XD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
P6 := Divisor(XD10, Ideal([3*R.1*R.6 + R.6^2, R.2*R.6, 3*R.3*R.6 + R.6^2, 9*R.4^2 - 22*R.6^2, 9*R.5^2 - 22*R.6^2]));
S := [P1, P2, P3, P4, P5, P6];

*/

// We apply Symmetric Chabauty

intrinsic AnnihilatorSpaceOfDifferentialsModp(Xp::Crv, phi::MapAutSch) -> SeqEnum
{
	Xp: a curve over a finite field
	phi: the Atkin-Lehner involution on X (not on the reduction curve Xp)

	Output: A sequence of diffentials that are the annihilator of JXp(Fp) 
}
	Rp := CoordinateRing(AmbientSpace(Xp));
	phi_modp_equations := [Evaluate(f, [Rp.i: i in [1..Rank(Rp)]]) : f in DefiningPolynomials(phi)];
	phi_modp := iso<Xp->Xp | phi_modp_equations, phi_modp_equations>;

	Vom, VomToOmega := SpaceOfDifferentialsFirstKind(Xp);
	vanishing_differentials_modp := hom<Vom -> Vom | [(Pullback(phi_modp, VomToOmega (Vom.i)))@@VomToOmega - Vom.i : i in [1..Genus(Xp)]]>;
	vanishing_differentials_modp := Image(vanishing_differentials_modp);
	vanishing_differentials_modp := [VomToOmega(omega) : omega in Basis(vanishing_differentials_modp)];
	
	return vanishing_differentials_modp;
end intrinsic;


intrinsic SinglePointOnResidueClass(divisor::DivCrvElt, phi::MapAutSch, p::RngIntElt) -> BoolElt
{
	divisor: a degree 2 divisor
	phi: the Atkin-Lehner involution
	p: a rational prime;

	Output: True if Pt is the only point in the residue class of the Symmetric Power modulo p, else False
}
	assert Degree(divisor) le 2;
	QQ := RationalsAsNumberField();
	X := Curve(Parent(divisor));
	Pt := Decomposition(divisor)[1][1];
	K := ResidueClassField(Pt);
	L := QuadraticField(SquareFree(Discriminant(K)));
	iso, KtoL := IsIsomorphic(K, L);
	Pt := [KtoL(ai) : ai in Eltseq(RepresentativePoint(Pt))];
	pts := [[emb(ai) : ai in Pt] : emb in Automorphisms(L)];
	assert iso;

	OL := RingOfIntegers(L);
	pr := Factorization(p*OL)[1][1];
	uni_pr := UniformizingElement(pr);
	fpr := InertiaDegree(pr);
	Fp, OKtoFp := ResidueClassField(pr);
	Xp := ChangeRing(X, Fp);
	V := AnnihilatorSpaceOfDifferentialsModp(Xp, phi);	
	A := [];
	for pt in pts do
		min_val := Minimum([Valuation(ai, pr) : ai in pt | not IsZero(ai)]);
		Qred := [uni_pr^(-min_val)*ai : ai in pt];
		Qt := Xp![OKtoFp(ai) : ai in Qred];
		tQ := UniformizingParameter(Qt);
		rowA := [(omega/Differential(tQ))(Qt) : omega in V];
		Append(~A, rowA);
	end for;
		
	A := Matrix(A);
	if Rank(A) eq 1 then
		return true;
	end if;
	
	return false;
end intrinsic;


intrinsic PointsLieOnSingleResidueClass(X::Crv, phi::MapAutSch, sym_pts::SeqEnum, primes::SeqEnum) -> SeqEnum
{
	X: a curve over QQ
	phi: the Atkin-Lehner involution
	sym_pts: a sequence of rational degree 2 divisors that represent points on X^(2)(Q)
	primes: a sequence of primes over we will apply the relative symmetric Chabauty method.

	Output: A sequence of pairs [p, bol], when p lies in primes and bol is True if every point in S lies in a single residue class modulo p. 
}
	chabauty_primes := [];
	for p in primes do
		bol := true;
		for sym_pt in sym_pts do
			if not SinglePointOnResidueClass(sym_pt, phi, p) then
				bol := false;
				break;
			end if;
		end for;
		Append(~chabauty_primes, p);
	end for;

	return chabauty_primes;	
end intrinsic;



// Confirm computations 
//pairs := PointsLieOnSingleResidueClass(XD10, w11, S, [5, 7, 13]);


