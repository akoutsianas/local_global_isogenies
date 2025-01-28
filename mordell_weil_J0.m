////// In this file we compute the rank of J0 over Q and K using the isogeny J0 to the prodcut of Efi //////


// Conductor 121 //
Ef1 := EllipticCurve("121a");
Ef2 := EllipticCurve("121b");
Ef3 := EllipticCurve("121c");
Ef4 := EllipticCurve("121d");

// Conductor 11 //
Ef5 := EllipticCurve("11a");


// Rank computations over Q //
r1, bol1 := Rank(Ef1);
assert bol1;
assert r1 eq 0;

r2, bol2 := Rank(Ef2);
assert bol2;
assert r2 eq 1;

r3, bol3 := Rank(Ef3);
assert bol3;
assert r3 eq 0;

r4, bol4 := Rank(Ef4);
assert bol4;
assert r4 eq 0;

r5, bol5 := Rank(Ef5);
assert bol5;
assert r5 eq 0;


// Rank computations over K //
K<r11> := QuadraticField(-11);

Ef1K := ChangeRing(Ef1, K);
Ef2K := ChangeRing(Ef2, K);
Ef3K := ChangeRing(Ef3, K);
Ef4K := ChangeRing(Ef4, K);
Ef5K := ChangeRing(Ef5, K);

r1, bol1 := Rank(Ef1K);
assert bol1;
assert r1 eq 0;

r2, bol2 := Rank(Ef2K);
assert bol2;
assert r2 eq 2;

r3, bol3 := Rank(Ef3K);
assert bol3;
assert r3 eq 0;

r4, bol4 := Rank(Ef4K);
assert bol4;
assert r4 eq 0;

r5, bol5 := Rank(Ef5K);
assert bol5;
assert r5 eq 0;





