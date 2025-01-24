/////// In this file we determine a finite index subgroup of the Mordell-Weil group of the Jacobian of XD10 ///////
////// Load the model of XD10 from file models_and_maps_computations.m //////



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


////// Torsion part //////

Gs := [TorsionSubgroupOfJacobian(XD10, p) : p in [5, 7, 13, 17, 19]];
bound_order := GCD([Order(G) : G in Gs]);
assert bound_order eq 25;

//  We check that the case C25 can not happen //
invs := TorsionInvariants(Gs[1]);
assert inv in invs | (inv mod 25) eq 0] eq 0;