/////// We compute the Mordell-Weil group and the rational points on XpD10 ///////


R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(x^6 - 6*x^5 + 11*x^4 - 8*x^3 + 11*x^2 - 6*x + 1);
J := Jacobian(C);


//// Mordell-Weil group ////

G, GtoJ, finite_index, proof := MordellWeilGroup(J);
JtoG := Inverse(GtoJ);
print G;
printf "Torsion generator: %o, Free part generator: %o\n", GtoJ(G.1), GtoJ(G.2); 
assert finite_index;
assert proof;

// A basis of the Mordell-Weil group over Q //
assert Order(JtoG(C![0, -1] - C![1, -1, 0])) eq 5;
assert JtoG(C![1, 1, 0] - C![1, -1, 0]) eq G.2;


// Compute the rational points of XpD10 using the Chabauty method //

ptsChabauty := Chabauty(GtoJ(G.2));
assert # ptsChabauty eq 6;



//// Computations over K=Q(sqrt(-11)) ////
// Attach Siksek's code AttachSpec(path_to_file/g2-jac.m) and the code Attach("path_to_file/add.m");

K<r11> := NumberField(x^2 + 11);
_<y> := PolynomialRing(K);
CK := ChangeRing(C, K);
JK := BaseChange(J, K);

assert TorsionBound(JK, 10) eq 5;
assert RankBound(JK) eq 2;

// We prove that G2 and G3 are linear independent //
G2 := CK![1, 1, 0] - CK![1, -1, 0];
G3 := JK![y^2 + y + 1, -r11];

// G3 has infinite order //
assert Order(G3) eq 0;

// We prove that G2 and G3 are Z-linear independent //
assert independentModTorsion([G2, G3]);



