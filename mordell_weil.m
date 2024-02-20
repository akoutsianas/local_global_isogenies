// In this file we determine a finite index subgroup of the Mordell-Weil group of the Jacobian of the modular curves.



// Useful functions


function ZetaFunctionNumerator(X, p, genus)
// intrinsic ZetaFunctionNumerator(X::Crv, p::RngIntElt, genus::RngIntElt) -> RngUPolElt
//
//	X: A model of the curve over the Q
//	p: a prime of good reduction
//	genus: the genus of X


	////////// We may not need it ///////////

	PS<t> := PowerSeriesRing(Rationals(), genus + 1);
	R<T> := PolynomialRing(Rationals());

	f := PS!0;
	for i in [1..genus] do
		bi := #Points(ChangeRing(X, FiniteField(p^i)));
		f +:= bi * t^i/i;
		//printf "%o\n", f;
	end for;
	F := (1 - t) * (1 - p*t) * Exp(f);
	printf "F: %o\n", F;	
	
	P1 := R!1 + p^genus * T^(2*genus);
	for i in [1..genus] do
		ai := Coefficient(F, i);
		P1 +:= ai*T^i;
		if i ne genus then
			P1 +:= ai * p^(genus-i) * T^(2*genus-i);		
		end if;
	end for;

	return P1;
end function;
//end intrinsic;


function TorsionOrderBound(X, primes)

	L := [];
	for p in primes do
		Xp := ChangeRing(X, GF(p));
		if IsNonSingular(Xp) then
			P1 := Numerator(ZetaFunction(ChangeRing(X, FiniteField(p))));
			Append(~L, Evaluate(P1, 1));
		else
			printf "X is not singular at %o.\n", p;
		end if;
	end for;

	return GCD(L);

end function;



function TorsionSubgroupOfJacobian(X, p)

	Xp := ChangeRing(X, GF(p));
	is_singular := false;
	for i in [1..NumberOfAffinePatches(Xp)] do
		if IsSingular(AffinePatch(Xp, i)) then
			is_singular := true;
			break;
		end if;
	end for;
	if not is_singular then
		H, phi, psi := ClassGroup(Xp);
		C1 := FreeAbelianGroup(1);
		deg := hom<H -> C1 | [Degree(phi(a)) * C1.1 : a in OrderedGenerators(H)]>;
		Jp := Kernel(deg); 

		return Jp;
	else
		printf "The curve is not singular at %o.\n", p;
		return false;
	end if;
end function;



intrinsic ReduceDivModp(Xp::Crv, D::DivSchElt) -> DivCrvElt
{
	
	Xp: a curve over Fp
	D: a 
}
	p := Characteristic(BaseRing(Xp));
	Rp := CoordinateRing(AmbientSpace(Xp));
	Dp := DivisorGroup(Xp) ! 0;
	Fp := GF(p);
	for Di in IdealFactorisation(D) do
		I := Di[1];
		ai := Integers()!Di[2];
		Ip := [];
		for f in Generators(I) do
			val := AbsoluteValue(Min([Valuation(ci, p) : ci in Coefficients(f)]));
			if val gt 0 then
				Append(~Ip, Rp!(p^val * f));
			else
				Append(~Ip, Rp!f);
			end if;
		end for;

		Dip := Divisor(Xp, Ideal(Ip));
		Dp +:= ai * Dip;
	end for;
	
	return Dp;
end intrinsic;


intrinsic LinearIndependent(X::Crv, divs::SeqEnum, primes::SeqEnum, torsion_order::RngIntElt) -> BoolElt
{
	X: a curve over Q
	divs: a sequence of 0 rational divisors
	primes: the primes of good reduction that we use
	torsion_order: the order of the torsion divides this number

	True in case we probe that divs a Z-linear idependent in J(Q)/Tor, otherwise false
}
	
	G := FreeAbelianGroup(#divs);
	T := G;
	torsion_order_divs := PrimeFactors(torsion_order);

	for p in primes do
		Xp := ChangeRing(X, GF(p));
		divsp := [ReduceDivModp(Xp, D) : D in divs];
		Cp, psi, phi := ClassGroup(Xp);
		GtoCp := hom<G -> Cp | [phi(Dip) : Dip in divsp]>;
		T := T meet Kernel(GtoCp);
		T := sub<G | [G!g : g in Generators(T)]>;
		M := GCD([GCD(Eltseq(G!g)) : g in Generators(T)]);
		printf "M:%o, T: %o\n", M, T;

		if not(Set(PrimeDivisors(M)) subset torsion_order_divs) then
			return true;
		end if;
	end for;

	return false;

end intrinsic;


intrinsic ReducePtModp(Xp::Crv, P::Pt) -> Pt
{
	Xp: the reduction mod p of a projective curve X
	P: a rational point of X(Q)

	The reduction of P mod p
}

	Fp := BaseRing(Xp);
	p := Characteristic(Fp);
	val := AbsoluteValue(Min([Valuation(c, p) : c in Coordinates(P)]));
	return Xp![Fp!(c*p^val) : c in Coordinates(P)];
end intrinsic;


intrinsic LinearIndependentFromPts(X::Crv, pts::SeqEnum, primes::SeqEnum, torsion_order::RngIntElt, r::RngIntElt) -> SeqEnum, BoolElt
{
	X: a curve over Q
	pts: points on X(Q)
	primes: a set of primes of good reduction that we use
	torsion_order: the order of the torsion divides this number
	r: The number of Pi - Pj that we want to prove that are Z-linear independent
	
}
	subs := SetToSequence(Subsets(Set(pts), r + 1));
	G := FreeAbelianGroup(r);
	T := G;
	torsion_order_divs := PrimeFactors(torsion_order);

	SubsRec := recformat<points: SeqEnum, ker: GrpAb>;
	CpRec := recformat<Xp: Crv, Cp: GrpAb, psi: Map, phi: Map>;

	subs_recs := [];
	for sub in subs do
		rec := rec<SubsRec | points := SetToSequence(sub), ker := T>;
		Append(~subs_recs, rec);
	end for;

	jac_recs := [];
	for p in primes do
		Xp := ChangeRing(X, GF(p));
		Cp, psi, phi := ClassGroup(Xp);
		rec := rec<CpRec | Xp := Xp, Cp := Cp, psi := psi, phi := phi>;
		Append(~jac_recs, rec);
	end for;
	
	for sub in subs_recs do
		printf "sub: \n\n", sub;
		points := sub`points;
		T := sub`ker;
		for jac in jac_recs do
			Xp := jac`Xp;
			printf "p: %o\n", Characteristic(BaseRing(Xp));
			Cp := jac`Cp;
			phi := jac`phi;
			ptsmodp := [ReducePtModp(Xp, P) : P in points];
			divsp := [Place(ptsmodp[i]) - Place(ptsmodp[r + 1]) : i in [1..r]];
			GtoCp := hom<G -> Cp | [phi(Dip) : Dip in divsp]>;
			T := T meet Kernel(GtoCp);
			T := sub<G | [G!g : g in Generators(T)]>;
			M := GCD([GCD(Eltseq(G!g)) : g in Generators(T)]);
			printf "M:%o\n", M;
			printf "T: %o\n\n", T;
			
			if not(Set(PrimeDivisors(M)) subset torsion_order_divs) then
				return points, true;
			end if;
		end for;
	end for;

	return [], false;
end intrinsic;



function ComputeFracs(B)
        fracs := [];
        for a in [-B..B] do
                for b in [1..B] do
                        if GCD(a,b) eq 1 then
                                Append(~fracs, a/b);
                        end if;
                end for;
        end for;
        return fracs;
end function;

function LiftPoints(a, phi, fracs)
        p := Characteristic(Codomain(phi));
        rk := NumberField(Domain(phi)).1;
        ha := a@@phi;
        ha0 := ha[1];
        ha1 := ha[2];
        Ha0 := [];
        Ha1 := [];
        for frac in fracs do
                Append(~Ha0, ha0 + p*frac);
                Append(~Ha1, ha1 + p*frac);
        end for;
        Ha := [];
        for a0 in Ha0 do
                for a1 in Ha1 do
                        Append(~Ha, a0 + rk*a1);
                end for;
        end for;
        return Ha;
end function;

function LiftPointsJ(J, p : B := 2, d := 2)
	K := BaseField(J);
	_<s> := PolynomialRing(K);
	OK := RingOfIntegers(K);
	fp := SetToSequence(Support(p*OK))[1];
	Fp, OKToFp := ResidueClassField(fp);
	Jp := BaseChange(J, Fp);
	ptsJp := Points(Jp);
	printf "ptsJp: %o\n", #ptsJp;
	
	fracs := ComputeFracs(B);
	printf "fracs: %o\n", #fracs;

	pts := [];
	// lift points of Jp to J
	for pt in ptsJp do
		print pt;
		if IsZero(pt) then
			continue;
		end if;
		ptx := Coefficients(pt[1]);
		if #ptx eq 1 then
			Append(~ptx, 0);
		end if;
		//return ptx[1], OKToFp, fracs;
		a0s := LiftPoints(ptx[1], OKToFp, fracs);
		a1s := LiftPoints(ptx[2], OKToFp, fracs);
		printf "a0s: %o, a1s: %o", #a0s, #a1s;
		for a0 in a0s do
			for a1 in a1s do
				f := s^2 + a1*s + a0;
				ptsf := RationalPoints(J, f, d);
				if not IsEmpty(ptsf) then
					for P in ptsf do
						Append(~pts, P);
					end for;
				end if;
			end for;
		end for;
	end for;
	return pts;
end function;

