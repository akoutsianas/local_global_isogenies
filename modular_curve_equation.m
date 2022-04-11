// AttachSpec("/cecm/home/akoutsia/Programming/Magma/ModFrmGL2/ModFrmGL2/ModFrmGL2.spec");



gens := [[4,0,0,3], [0,1,1,0], [2,0,0,2]];
N := 11;
H_N := sub<GL(2,Integers(N)) | gens>;
H := PSL2Subgroup(H_N);
M := ModularSymbols(H, 2, Rationals(), 0);
S := CuspidalSubspace(M);
XD10<[x]>, fs := ModularCurve(H);
AssignNames(~XD10, ["x0", "x1", "x2", "x3", "x4", "x5"]);
