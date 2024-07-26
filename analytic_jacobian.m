####### Computations in Jacobians of modular curves using analytic methods #######








##### test #####

prec := 50;
N := 67;
S := CuspidalSubspace(ModularSymbols(Gamma0(N), 2));
basisq := qIntegralBasis(S, prec);
// S := CuspidalSubspace(ModularForms(Gamma0(N), 2));

// X0(67) model
P<x,y,z,w,t> := ProjectiveSpace(Rationals(),4);
R := PolynomialRing(Rationals(), 5);


mons := [];
monsq := [];
for i in [1..5] do
	for j in [i..5] do
		Append(~mons, R.i * R.j);
		Append(~monsq, basisq[i] * basisq[j]);
	end for;
end for;

A := [];
for fq in monsq do
	Append(~A, [Coefficient(fq, i): i in [1..prec]]);
end for;
ker := Basis(Kernel(Matrix(A)));

I := [&+[v[i]*mons[i] : i in [1..#mons]]: v in ker];
C := Curve(P, I);
pts := PointSearch(C, 100);

OmC := SpaceOfHolomorphicDifferentials(C);
OmCbasis := BasisOfHolomorphicDifferentials(C);
KC<a,b,c,d> := FunctionField(C);



function qExpandFunctionFieldElement(f, basisq)
/*
	INPUT:
		f: a function field element
		basisq: the q-expansion of the modular forms

	OUTPUT:
		The q-expansion of f.	
*/
	fproj := ProjectiveFunction(f);

	numfproj := Numerator(fproj);
	mons_numfproj := Monomials(numfproj);
	coefs_numfproj := Coefficients(numfproj);
	numfq := &+[coefs_numfproj[i] * Evaluate(mons_numfproj[i], basisq) : i in [1..#mons_numfproj]];

	denomfproj := Denominator(fproj);
	mons_denomfproj := Monomials(denomfproj);
	coefs_denomfproj := Coefficients(denomfproj);
	denomfq := &+[coefs_denomfproj[i] * Evaluate(mons_denomfproj[i], basisq) : i in [1..#mons_denomfproj]];	

	fq := numfq / denomfq;

	return fq;
end function;


function qExpandHolomorphicDifferential(om, basisq)
/*
	INPUT:
		om: a differential
		basisq: the q-expansion of the modular forms
		
	OUTPUT:
		The q-expansion of om.
*/
	KC := FunctionField(Curve(om));
	var1 := KC.1;
	f1 := om / Differential(var1);	

	var1q := qExpandFunctionFieldElement(KC.1, basisq);
	f1q := qExpandFunctionFieldElement(f1, basisq);

	omq := f1q * Derivative(var1q);
	return omq;
end function;




function HolomorphicDifferentialChangeBasisMatrix(basisq, OmCbasis : Precision := 30)
/*
	INPUT:
		basisq: the q-expansion of the modular forms
		OmCbasis: the basis of homolomorphic differentials of C given by the function BasisOfHolomorphicDifferentials.


	OUTPUT:
		A square matrix that gives the change of basis of holomorphic differentials between basisq and the one we get from 
		Magma's function BasisOfHolomorphicDifferentials. A*OmCbasis = basisq.
*/	
	OmCbasisq := [qExpandHolomorphicDifferential(om, basisq) : om in OmCbasis];
	MOmCbasisq := Matrix(Rationals(), [[Coefficient(f, i) : i in [0..Precision]] : f in OmCbasisq]);
	Mbasisq := Matrix(Rationals(), [[Coefficient(f, i) : i in [1..Precision+1]] : f in basisq]);
	
	A := Solution(MOmCbasisq, Mbasisq);

	return A;
end function;	


function IntegralOfDifferentialAroundP0(om, P0, P : Precision := 30)
/*
	INPUT:
		om: a differential of C
		P0: a given point of C over the integration will take place
		P: a given point of C	

	OUTPUT:
		The integral \sum_P^{P0}om
*/
	C := Curve(om);
	KC := FunctionField(C);
	KCP0, KCtoKCP0 := Completion(KC, Place(P0) : Precision := Precision);
	uni := UniformizingParameter(P0);
	uniP := Evaluate(uni, P);

	assert uniP ne Infinity();

	f := om / Differential(uni);
	funi := KCtoKCP0(f);
		
	return Evaluate(Integral(funi), uniP);
end function;


function AbelJacobiMap(OmCbasis, P0, P : Precision := 30)
/*
	INPUT:
		OmCbasis: a basis of homolomorphic differentials of C.
		P0: a given point of C over the integration will take place
		P: a given point of C

	OUTPUT:
		The Abel-Jacobi map of [P0 - P] with respect to OmCbasis.
*/
	integrals := [IntegralOfDifferentialAroundP0(om, P0, P : Precision := Precision) : om in OmCbasis];
	return integrals;
end function;



function AbelJacobiWithqExpansion(basisq, P0, P : Precision := 30)
/*
	INPUT:
		basisq: The q-expansion of the space of modular forms. The cuspforms have to be with rational coefficients.
		P0: a given point of C over the integration will take place
		P: a given point of C

	OUTPUT:
		The Abel-Jacobi map of [P0 - P] with respect to basisq
*/
	C := Curve(P0);
	OmCbasis := BasisOfHolomorphicDifferentials(C);
	A := HolomorphicDifferentialChangeBasisMatrix(basisq, OmCbasis : Precision := Precision);
	AJ := AbelJacobiMap(OmCbasis, P0, P : Precision := Precision);
	AJ := Vector(AJ) * A;

	return AJ;
end function;






// Compute C^g/L and the image of [P - oo] where P is a rational point
periods := Periods(S, 300);

function ComputeTau(P, basisq : prec := 100, lam := 10^(-3))
/*
	INPUT:
		P: is a point on C
		basisq: the q-expansion of the integral basis of cuspforms

	OUTPUT:
		A point on the upper half plane that corresponds to P

	NOTE:
		P should not be a cusp.	
*/
	n := #basisq;
	for i in [1..n] do
		if not IsZero(P[i]) then
			i0 := i;
			break;
		end if;
	end for;
	// P0 := [P[i]/P[i0] : i in [1..n] | i ne i0];
	Fs := [basisq[i] * P[i0] - basisq[i0] * P[i] : i in [1..n] | i ne i0];
	//Fs := [basisq[i] /basisq[i0] -  P[i]/P[i0] : i in [1..n] | i ne i0];
	print([Valuation(f) : f in Fs]);

	// We apply Newton-Raphson method
	CC := ComplexField(prec);
	z0 := 0.5 + CC.1/10;
	Fsdiv := [Derivative(f) : f in Fs];
	k := 50; //Floor(prec/2);
	pi := Pi(CC);
	Ez := function(z) return Exp(2 * pi * CC.1 * z); end function;
	for i in [1..k] do
		Jz0 := Matrix(n-1, 1, [Evaluate(f, Ez(z0)) * (2 * pi * CC.1 * Ez(z0))  : f in Fsdiv]);
		Jplus := Inverse(Transpose(Jz0) * Jz0) * Transpose(Jz0);
		Fsz0 := Matrix(n-1, 1, [Evaluate(f, Ez(z0)) : f in Fs]);
		//print(Jplus * Fsz0);
		z0 := z0 - (Jplus * Fsz0)[1][1];
		printf "Sum of norms: %o\n\n", &+[Norm(Evaluate(fs, Ez(z0))) : fs in Fs];
		printf "z0: %o\n\n", z0;
	end for;
	
	//printf "Values at Fs: %o\n", &+[Evaluate(fs, z0) : fs in Fs];
	return z0;
end function;



