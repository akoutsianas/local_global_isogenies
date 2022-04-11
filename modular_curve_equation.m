// AttachSpec("/cecm/home/akoutsia/Programming/Magma/ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");



gens := [[4,0,0,3], [0,1,1,0], [2,0,0,2]];
N := 11;
H_N := sub<GL(2,Integers(N)) | gens>;
H := PSL2Subgroup(H_N);
M := ModularSymbols(H, 2, Rationals(), 0);
S := CuspidalSubspace(M);
// D := Decomposition(S, HeckeBound(S));
// B := qExpansionBasis(S, prec);
B := qIntegralBasis(S, 50);

//K<r11> := QuadraticField(-11);
R := ProjectiveSpace(Rationals(), 5);
AssignNames(~R, ["x0", "x1", "x2", "x3", "x4", "x5"]);
X := [R.i : i in [1..6]];
I := [];

for d in [2..5] do
	Cart := CartesianPower([1..#B], d);
	Bd := [];
	Mons := [];
	for v in Cart do
		f := &*[B[i] : i in v];
		mon := &*[X[i] : i in v];
		if f notin Bd then
			Append(~Bd, f);
			Append(~Mons, mon);
		end if;
	end for;

	//printf "AbsolutePrecision %o\n", [AbsolutePrecision(f) : f in Bd];

	min_prec := Min([AbsolutePrecision(f) : f in Bd]);
	//printf "min_prec %o\n", min_prec;
	Bd := [ChangePrecision(f, min_prec) : f in Bd];

	Q := [];
	for f in Bd do
               	coeffs := [Coefficient(f,i) : i in [1..min_prec-1]];
               	Append(~Q,coeffs);
       	end for;
	//M := Matrix(Integers(), #Bd, min_prec-1, Q);
	M := Matrix(Rationals(), #Bd, min_prec-1, Q);
	N := Nullspace(M);
	printf "d: %o, Nullspace %o\n", d, Dimension(N);
		
	I := [];
	for vN in Basis(N) do
		Append(~I, &+[(vN[i])*Mons[i] : i in [1..#Mons]]);
	end for;
	C := Scheme(R, I);
	printf "Dim(C) %o\n", Dimension(C);
end for;



