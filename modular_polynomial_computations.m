/////// We prove that the modular polynomial always has a linear factor //////

F11 := ClassicalModularPolynomial(11);
_<x,y> := Parent(F11);

Js := [-3375, 8000, -884736, 16581375, -884736000];

for j0 in Js do
	assert Degree(Factorization(Evaluate(F11, [x, j0]))[1][1]) eq 1;
end for;

