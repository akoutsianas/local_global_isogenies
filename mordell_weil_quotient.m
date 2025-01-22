/////// We compute the Mordell-Weil group and the rational points on XpD10 ///////


R<x> := PolynomialRing(Rationals());
XpD10 := HyperellipticCurve(-2*x^5 + 2*x^4 - 3*x^3 + 2*x^2 - 2*x, x^3 + x^2 + x + 1);
JpD10 := Jacobian(XpD10);


//// Mordell-Weil group ////

G, GtoJpD10, finite_index, proof := MordellWeilGroup(JpD10);
JpD10toG := Inverse(GtoJpD10);
print G;
printf "Torsion generator: %o, Free part generator: %o\n", GtoJpD10(G.1), GtoJpD10(G.2); 
assert finite_index;
assert proof;

// A basis of the Mordell-Weil group //
assert Order(JpD10toG(XpD10![0, -1] - XpD10![1, -1, 0])) eq 5;
assert JpD10toG(XpD10![1, 0, 0] - XpD10![1, -1, 0]) eq G.2;


// Compute the rational points of XpD10 using the Chabauty method //

C, XpD10toC := SimplifiedModel(XpD10); 
J := Jacobian(C);
G, GtoJ, _, _ := MordellWeilGroup(J);
ptsChabauty := Chabauty(GtoJ(G.2));
assert # ptsChabauty eq 6;
