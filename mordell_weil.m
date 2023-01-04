// In this file we determine a finite index subgroup of JD10 orver the rationals.

P<x,y,z,u,v,w> := ProjectiveSpace(Rationals(),5);

ID10 := [u*w - 2*v*w + 2*u*x - 6*v*x + 2*u*y + 2*v*y + u*z,
u*w + v*w + 2*u*x - 2*v*x + 2*u*y - 10*v*y - 5*u*z + 11*v*z,
- 6*u^2 + 6*u*v - 3*v^2 + 11*w^2 - 66*w*x + 11*x^2 + 88*w*y  - 110*x*y + 99*y^2 + 44*w*z - 110*x*z,
6*u^2 + 12*u*v + 12*v^2 + 187*w*x + 22*x^2 + 55*w*y - 44*x*y - 154*y^2 + 66*w*z + 77*x*z  + 121*y*z,
- 9*v^2 + 88*w^2- 11*w*x -99*x^2 - 77*w*y + 110*x*y - 11*y^2 + 77*w*z - 297*x*z  + 121*y*z,
- 6*u^2 - 12*u*v - 12*v^2 + 33*w^2 - 77*w*x + 66*x^2 - 121*w*y - 132*x*y - 110*y^2 - 44*w*z - 187*x*z + 121*y*z  + 121*z^2];
XD10 := Curve(P, ID10);
SD10 := Scheme(P, ID10);
phi := iso<XD10 -> XD10 | [-P.1, -P.2, -P.3, P.4, P.5, -P.6],[-P.1, -P.2, -P.3, P.4, P.5, -P.6]>;
G := AutomorphismGroup(XD10,[phi]);
C, XD10ToC := CurveQuotient(G);
CQQ := [C![1, -3, 1], C![1, -1, 1], C![0, -1, 1], C![0, 0, 1], C![1, -1, 0], C![1, 0 , 0]];
JC := Jacobian(C);
GJC, GJCToJC := MordellWeilGroup(JC);
JCToGJC := Inverse(GJCToJC);



// Compute points on JC over K = Q(sqrt(-11))

K<r11> := QuadraticField(-11);
CK := BaseChange(C, K);
JK := Jacobian(CK);
_<s> := PolynomialRing(K);
Csim, CToCsim := SimplifiedModel(C);
Jsim := Jacobian(Csim);
G2 := CK![1,0,0] - CK![1,-1,0];
G3 := JK![s^2 + s + 1, (-r11 - 1)/2];

CsimK := BaseChange(Csim, K);
JsimK := Jacobian(CsimK);
G2sim := CsimK!(CToCsim(C![1,0,0])) - CsimK!(CToCsim(C![1,-1,0]));
G3sim := JsimK![s^2 + s + 1, -r11];



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
pts := LiftPointsJ(JK, 7);


// Prove that the rank of JD10 over K is 2



// Finite index subgroup of the JD10(Q)


R := CoordinateRing(P);
D1a := Divisor(SD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
Dinf := Divisor(SD10, Ideal([R.1*R.6 + 1/3*R.6^2, R.2*R.6, R.3*R.6 + 1/3*R.6^2, R.4^2 - 22/9*R.6^2, R.5^2 - 22/9*R.6^2]));
D1 := D1a - Dinf;

D2a := Divisor(SD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 - 11*R.6^2, R.5^2 - 11*R.6^2]));
D2 := D2a - Dinf;
