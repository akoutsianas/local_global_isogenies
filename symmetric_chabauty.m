/*

// Attach("/cecm/home/akoutsia/Programming/Magma/local_global_isogenies/symmetric_chabauty.m");

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);

ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);
w11 := iso<XD10 -> XD10 | [-P.1, -P.2, -P.3, P.4, P.5, -P.6],[-P.1, -P.2, -P.3, P.4, P.5, -P.6]>;
G := AutomorphismGroup(XD10,[w11]);
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
Gens := SetToSequence({D1a, D2a});

*/

// We apply Relative Symmetric Chabauty

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



// Mordell-Weil Sieve


intrinsic DegreeTwoDivisorsOfXp(Xp::Crv) -> SeqEnum
{
	Xp: the reduction of X mod p

	Output: The degree 2 divisors of Xp.
}
	p := Characteristic(BaseRing(Xp));
	pts := [Divisor(Place(pt)) : pt in Points(Xp,GF(p^2)) | not IsSingular(pt)];
	divs_deg1 := [pt : pt in pts | Degree(pt) eq 1];
	divs_deg2 := [pt : pt in pts | Degree(pt) eq 2];
	
	Xp_sym2 := [pt1 + pt2 : pt1, pt2 in divs_deg1];
	for pt in divs_deg2 do
		if pt notin Xp_sym2 then
			Append(~Xp_sym2, pt);
		end if;
	end for;
	return Xp_sym2;
end intrinsic;


intrinsic ReduceDivisorModp(D::DivCrvElt, Xp::Crv) -> DivCrvElt
{
	D: a divisor of X
	Xp: the reduction of X mod p
	
	Output: The divisor of D mod p.
}
	X := Curve(Parent(D));
	p := Characteristic(BaseRing(Xp));
	I := Ideal(D);
	poly_gens := Basis(I) cat DefiningPolynomials(X);	

	assert Rank(CoordinateRing(Xp)) eq Rank(CoordinateRing(X));
	Zx := PolynomialRing(Integers(), Rank(CoordinateRing(Xp)));
	
	poly_gens_reduced := [];
	for f in poly_gens do
		fred := f * LCM([Denominator(c) : c in Coefficients(f)]);
		vp := Minimum([Valuation(c, p) : c in Coefficients(fred)]);
		fred := fred * p^(-vp);
		Append(~poly_gens_reduced, fred);
	end for;

	Ired := ideal<Zx | poly_gens_reduced>;
	Fpx := CoordinateRing(AmbientSpace(Xp));
	poly_gens_modp := [Evaluate(f, [Fpx.i : i in [1..Rank(Fpx)]]) : f in poly_gens_reduced];
	Imodp := ideal<Fpx | poly_gens_modp>;
	Dp := Divisor(Xp, Imodp);
	
	return Dp;
end intrinsic;


intrinsic MordellWeilSieveSinglePrimeInfo(X::Crv, Gens::SeqEnum, G::GrpAb, sym_pts::SeqEnum, p::RngIntElt, N::RngIntElt) -> SeqEnum, GrpAb
{
	X: a curve over QQ
	Gens: a sequence of divisors [P1, P2, ..., Pn] such that [P1 - Pinf, P2 - Pinf, ..., Pn - Pinf] generate a subgroup of J(Q) of finite index.
	G: the finite index subgroup of J(Q) as an abstract abelian group.
	sym_pts: a sequence of rational degree 2 divisors that represent points on X^(2)(Q). The last divisor is the point at infinity.
	p: a prime number
	N: a positive integer such that N*J(Q) lies in the subgroup generated by Gens

	Output: Representatives of the cosets that correspond to preimages of Mp under phi_p and the kernel of phi_p map.
}
	Pinf := sym_pts[#sym_pts];	
	Fp := GF(p);
	Xp := ChangeRing(X, Fp);
	ClXp, ClXpToDivXp, DivXpToClXp := ClassGroup(Xp);
	Z := AbelianGroup([0]);
	deg_map := hom<ClXp -> Z | [Degree(ClXpToDivXp(a)) : a in OrderedGenerators(ClXp)]>;
	JFp := Kernel(deg_map);
	sym_pts_modp := [ReduceDivisorModp(D, Xp) : D in sym_pts];
	Gens_modp := [ReduceDivisorModp(D, Xp) : D in Gens];
	Pinf_modp := sym_pts_modp[#sym_pts_modp];
	phi_p := hom<G -> JFp | [JFp!(DivXpToClXp(D) - DivXpToClXp(Pinf_modp)) : D in Gens_modp]>;
	Gmodp := Image(phi_p);

	// Determination of the set Mp
	Xp_sym2 := DegreeTwoDivisorsOfXp(Xp);
	im_ip := [N * JFp!(DivXpToClXp(D) - DivXpToClXp(Pinf_modp)) : D in Xp_sym2];

	Mp := [];	
	for pt in Xp_sym2 do
		if pt notin sym_pts_modp then
			Append(~Mp, pt);
		end if;
	end for;
	ip_Mp := [N * JFp!(DivXpToClXp(D) - DivXpToClXp(Pinf_modp)) : D in Mp];
	Mp_cosets := [pt@@phi_p : pt in ip_Mp | pt in Image(phi_p)];
	ker_phip := Kernel(phi_p);
		
	return Mp_cosets, ker_phip;
end intrinsic;


intrinsic MordellWeilSieve(X::Crv, Gens::SeqEnum, Gens_order::SeqEnum, sym_pts::SeqEnum, primes::SeqEnum, N::RngIntElt) -> BoolElt
{
	X: a curve over QQ
	Gens: a sequence of divisors [P1, P2, ..., Pn] such that [P1 - Pinf, P2 - Pinf, ..., Pn - Pinf] generate a subgroup of J(Q) of finite index.
	Gens_order: the order of each point of Gens
	sym_pts: a sequence of rational degree 2 divisors that represent points on X^(2)(Q). The last divisor is the point at infinity.
	primes: a sequence of primes over the relative symmetric Chabauty works for the set of points in sym_pts
	N: a positive integer such that N*J(Q) lies in the subgroup generated by Gens

	Output: True if we can apply the Mordell-Weil sieve for the set of quadratic points sym_pts and the given primes
}

	G := AbelianGroup(Gens_order);
	W, H := MordellWeilSieveSinglePrimeInfo(X, Gens, G, sym_pts, primes[1], N);
	if #W eq 0 then
		return true;
	end if;
	_, HtoG := sub<G | H>;
	
	for i in [2..#primes] do
		p := primes[i];
		printf "p: %o\n", p;
		Wp, Hp := MordellWeilSieveSinglePrimeInfo(X, Gens, G, sym_pts, p, N);
		_, HpToG := sub<G | Hp>;
		Hnew := H meet Hp;
		quoHnew, GToquoHnew := quo<G | Hnew>;
		quoH := GToquoHnew(H);
		quoHp := GToquoHnew(Hp);
		cosetsW := {a + b : a in Set(quoH), b in [GToquoHnew(c) : c in W]};
		cosetsWp := {a + b : a in Set(quoHp), b in [GToquoHnew(c) : c in Wp]};
		W := [pt@@GToquoHnew : pt in cosetsW meet cosetsWp];
		H := Hnew;
		printf "W: %o\n", #W;
		printf "quoHnew: %o\n", quoHnew;
		if #W eq 0 then
			return true;
		end if;
	end for;

	return false;
end intrinsic;


// The code someone should run
// primes := PointsLieOnSingleResidueClass(XD10, w11, S, [5, 7, 13, 17]);
// MordellWeilSieve(XD10, Gens, [5, 0], S, primes, 2);





