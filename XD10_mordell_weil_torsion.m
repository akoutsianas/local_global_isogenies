/////// In this file we determine a finite index subgroup of the Mordell-Weil group of the Jacobian of XD10 ///////

/////// Equation for XD10 /////
ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);




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
