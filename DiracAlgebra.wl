(* ::Package:: *)

BeginPackage["DiracAlgebra`"];

DiracIdentity2::usage =
  "DiracIdentity2 is the 2 x 2 identity matrix.";
DiracIdentity4::usage =
  "DiracIdentity4 is the 4 x 4 identity matrix.";
PauliMatrixSet::usage =
  "PauliMatrixSet is {sigma1,sigma2,sigma3}.";
MinkowskiMetric::usage =
  "MinkowskiMetric is diag(1,-1,-1,-1).";

SigmaMatrix::usage =
  "SigmaMatrix[mu] gives sigma^mu=(I,sigma1,sigma2,sigma3), mu=0,...,3.";
SigmaBarMatrix::usage =
  "SigmaBarMatrix[mu] gives sigmaBar^mu=(I,-sigma1,-sigma2,-sigma3).";
GammaMatrix::usage =
  "GammaMatrix[mu] gives the contravariant Weyl-basis gamma^mu, mu=0,...,3.";
GammaLowerMatrix::usage =
  "GammaLowerMatrix[mu] gives gamma_mu=g_(mu nu) gamma^nu.";
Gamma5Matrix::usage =
  "Gamma5Matrix is i gamma^0.gamma^1.gamma^2.gamma^3 in the Weyl basis.";

MinkowskiDot::usage =
  "MinkowskiDot[a,b] evaluates a^mu g_(mu nu) b^nu.";
DiracSlash::usage =
  "DiracSlash[p] evaluates pslash=p^mu gamma_mu for a four-vector p.";
DiracAdjoint::usage =
  "DiracAdjoint[m] evaluates gamma^0.ConjugateTranspose[m].gamma^0.";
MatrixProduct::usage =
  "MatrixProduct[{m1,m2,...}] evaluates the ordered matrix product m1.m2....";

GammaAnticommutator::usage =
  "GammaAnticommutator[mu,nu] evaluates {gamma^mu,gamma^nu}.";
Gamma5Anticommutator::usage =
  "Gamma5Anticommutator[mu] evaluates {gamma5,gamma^mu}.";
CheckDiracAlgebra::usage =
  "CheckDiracAlgebra[] returns an Association of explicit Dirac-algebra checks.";
AssertDiracAlgebra::usage =
  "AssertDiracAlgebra[] returns the checks if all pass and $Failed otherwise.";

Begin["`Private`"];

DiracIdentity2 = IdentityMatrix[2];
DiracIdentity4 = IdentityMatrix[4];
zero2 = ConstantArray[0, {2, 2}];

PauliMatrixSet = {
   {{0, 1}, {1, 0}},
   {{0, -I}, {I, 0}},
   {{1, 0}, {0, -1}}
   };

sigmaMatrices = Prepend[PauliMatrixSet, DiracIdentity2];
sigmaBarMatrices = Prepend[-PauliMatrixSet, DiracIdentity2];

SigmaMatrix[mu_Integer /; 0 <= mu <= 3] :=
  sigmaMatrices[[mu + 1]];

SigmaBarMatrix[mu_Integer /; 0 <= mu <= 3] :=
  sigmaBarMatrices[[mu + 1]];

MinkowskiMetric = DiagonalMatrix[{1, -1, -1, -1}];

gammaMatrices = Table[
   ArrayFlatten[{
     {zero2, SigmaMatrix[mu]},
     {SigmaBarMatrix[mu], zero2}
     }],
   {mu, 0, 3}
   ];

GammaMatrix[mu_Integer /; 0 <= mu <= 3] :=
  gammaMatrices[[mu + 1]];

GammaLowerMatrix[mu_Integer /; 0 <= mu <= 3] :=
  MinkowskiMetric[[mu + 1, mu + 1]] GammaMatrix[mu];

gamma5FromDefinition =
  I GammaMatrix[0] . GammaMatrix[1] . GammaMatrix[2] . GammaMatrix[3];

Gamma5Matrix = gamma5FromDefinition;

MinkowskiDot[a_List, b_List] /; Length[a] == 4 && Length[b] == 4 :=
  a . MinkowskiMetric . b;

DiracSlash[p_List] /; Length[p] == 4 := Sum[
   p[[mu + 1]] GammaLowerMatrix[mu],
   {mu, 0, 3}
   ];

DiracAdjoint[matrix_?MatrixQ] :=
  GammaMatrix[0] . ConjugateTranspose[matrix] . GammaMatrix[0];

MatrixProduct[mats_List] /; Length[mats] >= 1 := Dot @@ mats;

GammaAnticommutator[
   mu_Integer /; 0 <= mu <= 3,
   nu_Integer /; 0 <= nu <= 3
   ] := GammaMatrix[mu] . GammaMatrix[nu] +
  GammaMatrix[nu] . GammaMatrix[mu];

Gamma5Anticommutator[mu_Integer /; 0 <= mu <= 3] :=
  Gamma5Matrix . GammaMatrix[mu] + GammaMatrix[mu] . Gamma5Matrix;

zeroMatrixQ[matrix_?MatrixQ] := And @@ Flatten[
    PossibleZeroQ /@ matrix
    ];

CheckDiracAlgebra[] := Module[
  {gamma5Definition, checks},

  gamma5Definition =
   I GammaMatrix[0] . GammaMatrix[1] . GammaMatrix[2] . GammaMatrix[3];

  checks = <|
    "MetricSignature" ->
     TrueQ[MinkowskiMetric == DiagonalMatrix[{1, -1, -1, -1}]],

    "PauliAnticommutator" -> And @@ Flatten@Table[
       zeroMatrixQ[
        PauliMatrixSet[[i]] . PauliMatrixSet[[j]] +
         PauliMatrixSet[[j]] . PauliMatrixSet[[i]] -
         2 KroneckerDelta[i, j] DiracIdentity2
        ],
       {i, 1, 3}, {j, 1, 3}
       ],

    "CliffordAnticommutator" -> And @@ Flatten@Table[
       zeroMatrixQ[
        GammaAnticommutator[mu, nu] -
         2 MinkowskiMetric[[mu + 1, nu + 1]] DiracIdentity4
        ],
       {mu, 0, 3}, {nu, 0, 3}
       ],

    "Gamma5Definition" -> TrueQ[Gamma5Matrix == gamma5Definition],

    "Gamma5ExplicitWeylForm" -> TrueQ[
      Gamma5Matrix == DiagonalMatrix[{-1, -1, 1, 1}]
      ],

    "Gamma5Squared" -> zeroMatrixQ[
      Gamma5Matrix . Gamma5Matrix - DiracIdentity4
      ],

    "Gamma5Anticommutator" -> And @@ Table[
       zeroMatrixQ[Gamma5Anticommutator[mu]],
       {mu, 0, 3}
       ],

    "Gamma5Hermitian" -> TrueQ[
      Gamma5Matrix == ConjugateTranspose[Gamma5Matrix]
      ],

    "Gamma0Hermitian" -> TrueQ[
      GammaMatrix[0] == ConjugateTranspose[GammaMatrix[0]]
      ],

    "SpatialGammaAntiHermitian" -> And @@ Table[
       TrueQ[GammaMatrix[mu] == -ConjugateTranspose[GammaMatrix[mu]]],
       {mu, 1, 3}
       ],

    "Gamma5Trace" -> PossibleZeroQ[Tr[Gamma5Matrix]],

    "GammaTraces" -> And @@ Table[
       PossibleZeroQ[Tr[GammaMatrix[mu]]],
       {mu, 0, 3}
       ]
    |>;

  Append[checks, "AllPassed" -> And @@ Values[checks]]
  ];

AssertDiracAlgebra::failed =
  "One or more Dirac-algebra checks failed. Inspect the returned Association.";

AssertDiracAlgebra[] := Module[{checks = CheckDiracAlgebra[]},
  If[TrueQ[checks["AllPassed"]],
   checks,
   Message[AssertDiracAlgebra::failed];
   $Failed
   ]
  ];

End[];
EndPackage[];

