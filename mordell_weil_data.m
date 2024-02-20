// Attach("/cecm/home/akoutsia/Programming/Magma/local_global_isogenies/mordell_weil.m");

// In this file we have the Mordell-Weil Data of the curves


// The modular curve XD10

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



// Finite index subgroup of the JD10(Q)


R := CoordinateRing(P);
D1a := Divisor(SD10, Ideal([R.1*R.6 + R.6^2, R.2*R.6 - 7*R.6^2, R.3*R.6 - 5*R.6^2, R.4^2 - 473*R.6^2, R.5^2 - 473*R.6^2]));
Dinf := Divisor(SD10, Ideal([R.1*R.6 + 1/3*R.6^2, R.2*R.6, R.3*R.6 + 1/3*R.6^2, R.4^2 - 22/9*R.6^2, R.5^2 - 22/9*R.6^2]));
D1 := D1a - Dinf;

D2a := Divisor(SD10, Ideal([R.1*R.6 - R.6^2, R.2*R.6 - R.6^2, R.3*R.6 - R.6^2, R.4^2 - 11*R.6^2, R.5^2 - 11*R.6^2]));
D2 := D2a - Dinf;



// The modular curve Xsp19

P8<x,y,z,w,t,u,v,r,s> := ProjectiveSpace(Rationals(), 8);

model := [8*x^2+2*x*y+y^2+2*x*z-y*z-2*x*w+2*x*t-y*t+z*t+2*x*u+y*u-w*u+t*u-u^2-3*x*v-2*y*v+2*z*v-2*t*v+2*u*v-v^2+x*r-4*y*r+z*r+2*t*r+v*r+2*r^2+z*s-v*s-s^2,x^2-4*x*y+3*x*z+y*z+z^2-2*x*w-2*y*w-z*w-3*x*t+y*t+z*t-3*x*u-2*y*u+w*u-t*u+u^2-x*v-y*v+2*t*v-2*u*v+v^2+x*r+4*y*r+z*r-5*t*r+3*u*r-4*r^2-6*x*s+y*s+t*s-u*s-v*s+2*r*s,3*x^2-x*y+3*y^2-3*x*z-y*z-z^2+4*x*w+3*y*w-3*x*t-2*z*t+4*x*u+2*y*u-u^2-y*v+2*z*v+t*v+2*u*v+v^2-x*r+4*y*r-t*r+2*u*r-v*r-r^2-2*x*s+3*y*s-5*z*s+2*t*s-4*u*s+3*v*s-r*s+s^2,4*x^2-x*y+y^2+2*x*z-y*z-z^2+y*w+2*z*w-4*x*t+z*t-x*u-y*u+x*v-5*y*v+2*t*v+5*x*r-5*y*r+2*t*r+2*u*r+3*v*r+r^2+3*y*s+3*z*s-3*t*s-u*s-4*v*s,x^2-x*y+y^2+4*x*z-2*y*z+3*z^2-y*w-2*z*w-7*x*t-4*y*t+2*z*t+2*x*u-z*u+w*u+t*u-2*x*v-2*y*v-u*v+2*x*r-5*y*r+2*z*r-2*t*r+u*r+3*v*r-r^2-4*x*s+z*s+t*s+u*s-4*v*s+r*s,x^2-x*z-y*z-z^2-4*x*w+4*x*t+2*y*t+t^2-x*u-y*u-z*u-w*u+t*u-u^2-3*x*v-3*y*v-z*v-t*v+2*u*v-v^2+4*x*r+y*r-3*z*r-4*u*r+v*r+x*s-y*s+3*z*s+t*s+3*u*s-2*v*s+r*s,x^2-x*y+y^2-2*x*z+y*z-2*z^2+y*w+3*z*w+2*x*t-3*y*t+3*z*t-t^2-x*u+z*u+2*x*v-2*t*v+2*u*v-2*v^2+7*x*r-3*y*r-z*r+3*t*r+2*u*r-2*v*r+2*r^2-x*s+y*s+2*z*s-2*t*s-2*u*s,x*y-2*y^2+4*x*z+2*y*z+z^2+z*w+x*t+5*y*t-2*z*t+t^2-3*y*u+3*z*u-w*u-2*t*u+u^2+2*y*v-z*v+t*v-4*u*v+v^2+y*r+z*r+t*r+2*u*r+v*r+5*x*s+2*z*s-3*t*s-u*s-s^2,3*x^2+2*x*z+2*y*z-z^2+x*w+y*w+z*w-4*x*t-2*y*t-w*t-t^2-2*x*u+y*u+2*w*u-t*u+u^2+4*x*v-2*y*v-z*v+3*t*v-u*v+v^2-x*r-y*r-2*z*r+v*r-x*s+2*y*s+z*s+w*s-2*t*s-v*s,x*y-2*x*z-z^2+2*x*w+4*y*w+3*z*w+y*t-3*z*t+3*y*u+z*u-w*u-t*u+2*x*v+y*v+4*z*v-t*v+u*v+4*x*r-y*r+5*t*r+3*u*r-v*r+3*r^2+3*y*s-2*z*s-w*s-3*t*s-5*u*s+3*v*s-3*r*s,x*y+y^2+x*z+2*y*z+z^2+x*w+y*w-3*z*w-x*t-2*y*t+2*z*t-w*t+t^2+5*x*u+2*y*u-2*z*u+t*u-u^2-5*x*v-2*t*v-u*v+v^2+2*x*r+3*y*r-2*t*r-u*r+v*r-r^2-2*x*s-2*y*s+z*s+t*s+u*s+v*s-s^2,2*x^2+x*y-2*y^2+x*z+2*y*z-x*w-y*w+2*z*w+4*x*t+2*y*t-2*z*t-w*t-4*x*u+y*u+w*u-t*u+u^2+x*v+y*v-z*v-t*v-2*u*v-x*r-w*r-2*t*r-3*u*r-r^2+2*x*s+z*s+w*s+2*u*s-v*s,x^2-y^2+3*x*z-y*z-3*x*w-2*y*w-z*w-2*x*t+4*y*t-4*z*t+w*t+t^2-x*u-w*u-2*x*v-2*y*v+z*v+w*v-t*v+u*v-z*r-w*r-t*r-u*r+v*r-r^2-2*x*s-y*s-w*s+3*t*s-v*s+r*s,6*x^2-x*y+2*y^2-3*x*z-3*y*z-z^2+2*z*w-4*x*t-2*z*t+3*x*u+z*u+t*u-u^2-3*x*v-5*y*v-w*v+3*t*v+u*v-2*x*r-y*r-z*r+w*r-u*r+2*v*r+y*s-z*s-w*s+2*t*s-2*v*s+s^2,2*x^2-2*y^2-x*z+y*z+x*w+3*y*w+3*z*w-2*x*t+y*t-3*z*t-w*t-4*x*u+y*u+z*u+w*u-2*t*u+u^2+2*x*v-y*v+z*v-2*w*v+4*t*v-u*v+y*r-z*r-w*r-t*r+2*u*r-v*r-2*r^2-2*x*s+4*y*s+z*s+w*s-2*t*s-2*u*s-v*s+s^2,2*y^2-x*z-2*y*z+z^2-x*w+3*y*w+z*w-x*t+3*y*t-3*z*t-w*t+t^2-3*x*u-t*u+2*x*v-3*y*v+w*v+t*v+u*v+v^2+3*x*r+2*y*r-z*r-2*t*r+u*r+v*r-2*r^2-x*s+4*y*s-2*u*s-v*s+r*s,3*x^2-2*x*y+3*y^2+2*x*z-2*z^2+3*x*w+w^2-2*x*t+y*t+z*t+3*x*u+w*u-t*u-3*x*v-3*y*v+2*z*v-t*v+u*v+5*x*r+2*z*r-w*r-t*r+4*u*r-r^2-3*x*s+2*y*s-z*s-4*u*s,4*x^2+x*y+2*y^2+4*x*z+z^2+x*w+3*y*w-3*z*w-3*x*t-2*y*t+z*t+w*t-t^2+2*x*u+3*y*u-z*u-4*x*v-3*y*v+3*z*v-w*v+t*v-u*v+v^2+x*r-2*y*r+2*z*r-2*w*r-2*t*r+2*u*r+2*v*r-r^2-5*x*s+3*y*s-z*s+2*w*s-2*u*s-v*s-r*s,6*x^2+x*y+4*x*z+2*y*z+z^2+x*w-w^2-y*t-z*t-t^2+2*x*u+4*y*u+z*u-w*u-3*x*v-y*v+2*z*v-w*v-u*v+v^2-4*x*r-y*r+3*z*r-3*t*r+u*r+v*r-r^2-4*x*s+2*t*s-u*s-s^2,5*x^2+y^2+x*z+z^2-2*x*w-2*w^2-x*t-z*t-x*u+y*u+z*u-w*u-3*x*v-3*y*v+2*z*v-w*v+u*v+v^2+x*r+y*r+w*r-t*r+2*v*r-4*x*s+z*s-w*s+t*s-s^2,2*x^2+3*x*y-x*z-2*y*z-2*y*w+x*t-2*z*t+2*w*t-t^2+x*u+z*u-2*x*v+2*y*v+3*z*v-w*v+2*u*v-v^2-3*x*r+y*r+2*w*r+t*r+u*r-2*v*r+r^2-y*s-4*z*s-w*s+2*t*s-2*u*s+3*v*s];

Xsp19 := Curve(P8, model);
Ssp19 := Scheme(P8, model);
pts := [Xsp19![-3/2, 1, 2, -3, -3/2, 1, -1/2, -1/2, 1], Xsp19![1/5, 3/5, 1/10, -7/10, -3/10, 9/10, 1/5, 2/5, 1], Xsp19![-1/10, 0, 1/10, -1/10, -2/5, -1/5, 3/10, 4/5, 1], Xsp19![-1/2, -4/5, -13/60, -53/60, -13/60, 11/20, -1/5, 1/15, 1], Xsp19![2, -3/2, 2, -1/2, -7/2, -2, -8, 15/2, 1], Xsp19![1/4, -1/4, 5/12, -1/6, 2/3, -1/2, 1/2, 1/3, 1]];


// Finite index subgroup of the JXsp19(Q)
R := CoordinateRing(P8);

// Divisors with respect the rational points
D1 := Divisor(Ssp19, Ideal([2*R.1 + 3*R.2, R.3 - 2*R.2, R.4 + 3*R.2, 2*R.5 + 3*R.2, R.6 - R.2, 2*R.7 + R.2, 2*R.8 + R.2, R.9 - R.2]));
D2 := Divisor(Ssp19, Ideal([R.2 - 3*R.1, 2*R.3 - R.1, 2*R.4 + 7*R.1, 2*R.5 + 3*R.1, 2*R.6 - 9*R.1, R.7 - R.1, R.8 - 2*R.1, R.9 - 5*R.1]));
D3 := Divisor(Ssp19, Ideal([R.1 + R.3, R.2, R.4 + R.3, R.5 + 4*R.3, R.6 + 2*R.3, R.7 - 3*R.3, R.8 - 8*R.3, R.9 - 10*R.3]));
D4 := Divisor(Ssp19, Ideal([30*R.2 - 48*R.1, 30*R.3 - 13*R.1, 30*R.4 - 53*R.1, 30*R.5 - 13*R.1, 30*R.6 + 33*R.1, 30*R.7 - 12*R.1, 30*R.8 + 4*R.1, R.9 + 2*R.1]));
D5 := Divisor(Ssp19, Ideal([2*R.1 - 4*R.9, 2*R.2 + 3*R.9, 2*R.3 - 4*R.9, 2*R.4 + R.9, 2*R.5 + 7*R.9, 2*R.6 + 4*R.9, 2*R.7 + 16*R.9, 2*R.8 - 15*R.9]));
D6 := Divisor(Ssp19, Ideal([12*R.1 - 3*R.9, 12*R.2 + 3*R.9, 12*R.3 - 5*R.9, 12*R.4 + 2*R.9, 12*R.5 - 8*R.9, 12*R.6 + 6*R.9, 12*R.7 - 6*R.9, 12*R.8 - 4*R.9]));

divs := [D1 - D2, D1 - D3];




// Prove that the rank of JD10 over K is 2

pts := LiftPointsJ(JK, 7);


