(* ::Title:: *)
(*Master integral for all q qbar spin-density elements*)

(*
  Conventions follow EntangleCollinsFF.pdf:

    g^(mu nu) = diag(1,-1,-1,-1), epsilon^(0123) = +1,
    gamma^mu  = {{0,sigma^mu},{sigmaBar^mu,0}},
    gamma5    = diag(-I2,I2),
    (x,y,z)   = (-r,n,k).

  Notation used throughout:

    mq                         quark mass,
    pElectron, pPositron       incoming lepton momenta,
    kQuark, kAntiquark         outgoing parton momenta,
    quarkPolarizationVector    B_a,
    antiquarkPolarizationVector BBar_b,
    correlationMatrix          C_ab.

  A "Bar" suffix denotes the antiquark or the opposite jet; it does not
  denote complex conjugation.  Both spin indices use the common physical
  axes (x,y,z), including the antiquark index.

  The file does four things.

  1. It evaluates the transverse and longitudinal Dirac traces.
  2. It obtains the complete 3 x 3 correlation matrix C_ij.
  3. It constructs the general two-particle density matrix, including the
     quark vector B_i and antiquark vector BBar_j.
  4. It contracts B_i, BBar_j and every C_ij with the two dihadron spin
     analyzers.  The resulting masterIntegral can be used for arbitrary
     angular projections.

  The quantities h1AngleNormQ and h1AngleNormQBar keep the H1Angle normalization explicit;
  use dihadronNormalizationRules for the explicit correlator normalization,
  or absorbedH1NormalizationRules if it is included in H1Angle itself.
*)

ClearAll["Global`*"];

(* ---------------------------------------------------------------------- *)
(* 1. Import Dirac algebra and define e+ e- -> q qbar kinematics         *)
(* ---------------------------------------------------------------------- *)

diracPackagePath = FileNameJoin[{
    DirectoryName[ExpandFileName[$InputFileName]],
    "DiracAlgebra.wl"
    }];

Get[diracPackagePath];

diracAlgebraChecks = AssertDiracAlgebra[];
If[diracAlgebraChecks === $Failed, Abort[]];

(* Short local aliases keep the physics expressions compact; all matrix
   definitions and algebraic operations live in DiracAlgebra.wl. *)
id2 = DiracIdentity2;
id4 = DiracIdentity4;
pauli = PauliMatrixSet;
metric = MinkowskiMetric;
gammaUp = Table[GammaMatrix[mu], {mu, 0, 3}];
gammaDown = Table[GammaLowerMatrix[mu], {mu, 0, 3}];
gamma5 = Gamma5Matrix;

minkowskiDot[a_List, b_List] := MinkowskiDot[a, b];
diracSlash[p_List] := DiracSlash[p];
matrixProduct[mats_List] := MatrixProduct[mats];

kinematicAssumptions =
  Element[{s, mq, theta}, Reals] && s > 4 mq^2 && mq >= 0 &&
   0 <= theta <= Pi;

(* Keep energyCM and betaCM algebraic during the traces.  Substituting the
   square roots only after taking the traces is much faster than asking
   FullSimplify to manipulate radicals in every matrix element. *)

(* Incoming e-, incoming e+, outgoing q, outgoing qbar. *)
pElectron = energyCM {1, -Sin[theta], 0, Cos[theta]};
pPositron = energyCM {1, Sin[theta], 0, -Cos[theta]};
kQuark = energyCM {1, 0, 0, betaCM};
kAntiquark = energyCM {1, 0, 0, -betaCM};

kQuarkSlash = diracSlash[kQuark];
kAntiquarkSlash = diracSlash[kAntiquark];

(* L^(mu nu), with the irrelevant common normalization left out. *)
leptonicTensor = Table[
   pElectron[[mu]] pPositron[[nu]] +
    pElectron[[nu]] pPositron[[mu]] -
    metric[[mu, nu]] minkowskiDot[pElectron, pPositron],
   {mu, 1, 4}, {nu, 1, 4}
   ];

(*
  Sum L^(mu nu) Tr[gamma_mu left (kslash+m) gamma_nu right
                    (kbarslash-m)].
  "left" and "right" are the spin insertions in the amplitude and its
  conjugate.  This is the e+e- version of the generic Gamma traces in
  the transverse and longitudinal spin-trace definitions.
*)
contractedTrace[left_List, right_List] := Sum[
   leptonicTensor[[mu, nu]] Tr@matrixProduct@Join[
      {gammaDown[[mu]]},
      left,
      {kQuarkSlash + mq id4, gammaDown[[nu]]},
      right,
      {kAntiquarkSlash - mq id4}
      ],
   {mu, 1, 4}, {nu, 1, 4}
   ];

unpolarizedTrace = Factor[contractedTrace[{}, {}]];

(* i=1,2 are x,y; the third insertion is longitudinal. *)
transverseInsertion[i_Integer /; 1 <= i <= 2] :=
  {gammaUp[[i + 1]], gamma5};
longitudinalInsertion = {gamma5, gammaUp[[1]], gammaUp[[4]]};

correlationNumerator[i_Integer, j_Integer] /;
   1 <= i <= 3 && 1 <= j <= 3 := Which[
  i <= 2 && j <= 2,
  contractedTrace[transverseInsertion[i], transverseInsertion[j]],

  i == 3 && j <= 2,
  contractedTrace[longitudinalInsertion, transverseInsertion[j]],

  i <= 2 && j == 3,
  -contractedTrace[transverseInsertion[i], longitudinalInsertion],

  i == 3 && j == 3,
  -contractedTrace[longitudinalInsertion, longitudinalInsertion]
  ];

(* Fast on-shell reduction: betaCM^2=1-mq^2/energyCM^2 and s=4 energyCM^2. *)
kinematicReduce[expression_] := Module[
  {rational, numerator, denominator, reduced},
  rational = Together[TrigExpand[expression]];
  numerator = Expand[Numerator[rational]];
  denominator = Expand[Denominator[rational]];
  reduced = Factor@Cancel[
    (numerator/denominator) /. {
       betaCM^2 -> 1 - mq^2/energyCM^2,
       Sin[theta]^2 -> 1 - Cos[theta]^2
       } /. energyCM -> Sqrt[s]/2
    ];
  Cancel@Together[
    reduced /. Sin[theta]^2 -> 1 - Cos[theta]^2
    ]
  ];

correlationFromTraces = Table[
   kinematicReduce[correlationNumerator[i, j]/unpolarizedTrace],
   {i, 1, 3}, {j, 1, 3}
   ];

(*
  For a transverse antiquark insertion, the density-matrix matching sign
  cancels the sign in the transverse antiquark projector.  The longitudinal
  antiquark projector has no such sign, so the matching sign remains.
*)
quarkPolarizationNumerator[i_Integer] /; 1 <= i <= 3 := Which[
  i <= 2,
  contractedTrace[transverseInsertion[i], {}],
  i == 3,
  contractedTrace[longitudinalInsertion, {}]
  ];

antiquarkPolarizationNumerator[j_Integer] /; 1 <= j <= 3 := Which[
  j <= 2,
  contractedTrace[{}, transverseInsertion[j]],
  j == 3,
  -contractedTrace[{}, longitudinalInsertion]
  ];

quarkPolarizationFromTraces = Table[
   kinematicReduce[quarkPolarizationNumerator[i]/unpolarizedTrace],
   {i, 1, 3}
   ];

antiquarkPolarizationFromTraces = Table[
   kinematicReduce[antiquarkPolarizationNumerator[j]/unpolarizedTrace],
   {j, 1, 3}
   ];

(*
  Generic production interface.  initialTensor contracts the two vertex
  indices; leftVertices and rightVertices are the Dirac matrices on the
  conjugate-amplitude and amplitude sides.  This permits nonzero B vectors
  from parity-violating, polarized-initial-state, or absorptive dynamics.

  Example vertex list for real vector/axial couplings:
    vectorAxialVertices[gV,gA]

  Supply the corresponding (possibly nonsymmetric) initial-state tensor
  to spinParametersFromProduction.  Use Identity as the reducer for raw
  expressions or kinematicReduce for the present CM kinematics.
*)
vectorAxialVertices[gV_, gA_] := Table[
   gammaDown[[mu]] . (gV id4 - gA gamma5),
   {mu, 1, 4}
   ];

productionTrace[
   initialTensor_?MatrixQ,
   leftVertices_List,
   rightVertices_List,
   leftInsertion_List,
   rightInsertion_List
   ] := Sum[
   initialTensor[[mu, nu]] Tr@matrixProduct@Join[
      {leftVertices[[mu]]},
      leftInsertion,
      {kQuarkSlash + mq id4, rightVertices[[nu]]},
      rightInsertion,
      {kAntiquarkSlash - mq id4}
      ],
   {mu, 1, 4}, {nu, 1, 4}
   ];

spinParametersFromProduction[
   initialTensor_?MatrixQ,
   leftVertices_List,
   rightVertices_List,
   reduceFunction_: Identity
   ] := Module[{trace, denominator, bResult, bBarResult, cResult},

  trace[left_List, right_List] := productionTrace[
    initialTensor, leftVertices, rightVertices, left, right
    ];

  denominator = trace[{}, {}];

  bResult = Table[
    reduceFunction[
     Which[
       i <= 2, trace[transverseInsertion[i], {}],
       i == 3, trace[longitudinalInsertion, {}]
       ]/denominator
     ],
    {i, 1, 3}
    ];

  bBarResult = Table[
    reduceFunction[
     Which[
       j <= 2, trace[{}, transverseInsertion[j]],
       j == 3, -trace[{}, longitudinalInsertion]
       ]/denominator
     ],
    {j, 1, 3}
    ];

  cResult = Table[
    reduceFunction[
     Which[
       i <= 2 && j <= 2,
       trace[transverseInsertion[i], transverseInsertion[j]],

       i == 3 && j <= 2,
       trace[longitudinalInsertion, transverseInsertion[j]],

       i <= 2 && j == 3,
       -trace[transverseInsertion[i], longitudinalInsertion],

       i == 3 && j == 3,
       -trace[longitudinalInsertion, longitudinalInsertion]
       ]/denominator
     ],
    {i, 1, 3}, {j, 1, 3}
    ];

  <|
   "NormalizationTrace" -> reduceFunction[denominator],
   "B" -> bResult,
   "BBar" -> bBarResult,
   "C" -> cResult
   |>
  ];

(* ---------------------------------------------------------------------- *)
(* 2. Compact correlation matrix in the (x,y,z)=(-r,n,k) basis           *)
(* ---------------------------------------------------------------------- *)

angularFactor =
  1 + Cos[theta]^2 + (4 mq^2/s) Sin[theta]^2;

correlationMatrix = 1/angularFactor {
    {(1 + 4 mq^2/s) Sin[theta]^2,
     0,
     -(2 mq/Sqrt[s]) Sin[2 theta]},

    {0,
     -(1 - 4 mq^2/s) Sin[theta]^2,
     0},

    {-(2 mq/Sqrt[s]) Sin[2 theta],
     0,
     1 + Cos[theta]^2 - (4 mq^2/s) Sin[theta]^2}
    };

traceDerivationCheck = And @@ Flatten@Table[
     PossibleZeroQ[
      kinematicReduce[correlationFromTraces[[i, j]] - correlationMatrix[[i, j]]]
      ],
     {i, 1, 3}, {j, 1, 3}
     ];

(* For unpolarized e+e- -> gamma* -> q qbar at Born level, both local
   polarization vectors vanish.  They need not vanish for a general
   production vertex or initial state. *)
quarkPolarizationVector = quarkPolarizationFromTraces;
antiquarkPolarizationVector = antiquarkPolarizationFromTraces;

bornPolarizationCheck = And @@ Join[
    PossibleZeroQ /@ quarkPolarizationVector,
    PossibleZeroQ /@ antiquarkPolarizationVector
    ];

axes = {"x", "y", "z"};
correlationElements = Association@Flatten@Table[
    (axes[[i]] <> axes[[j]]) -> correlationMatrix[[i, j]],
    {i, 1, 3}, {j, 1, 3}
    ];

spinCorrelation[i_Integer, j_Integer] /;
   1 <= i <= 3 && 1 <= j <= 3 := correlationMatrix[[i, j]];

quarkPolarization[i_Integer] /; 1 <= i <= 3 := quarkPolarizationVector[[i]];
antiquarkPolarization[j_Integer] /; 1 <= j <= 3 := antiquarkPolarizationVector[[j]];

(* The general 4 x 4 two-particle spin-density matrix. *)
spinDensityMatrix[
   bVec_?VectorQ,
   bBarVec_?VectorQ,
   cMat_?MatrixQ
   ] := 1/4 (
    KroneckerProduct[id2, id2] +
     Sum[
      bVec[[i]] KroneckerProduct[pauli[[i]], id2],
      {i, 1, 3}
      ] +
     Sum[
      bBarVec[[j]] KroneckerProduct[id2, pauli[[j]]],
      {j, 1, 3}
      ] +
     Sum[
      cMat[[i, j]] KroneckerProduct[pauli[[i]], pauli[[j]]],
      {i, 1, 3}, {j, 1, 3}
      ]
    );

bSymbolic = Array[b, 3];
bBarSymbolic = Array[bBar, 3];
correlationSymbolic = Array[c, {3, 3}];

rhoGeneral = spinDensityMatrix[bSymbolic, bBarSymbolic, correlationSymbolic];

(* Explicit matrix in the ordered basis
   {|+,+>, |+,->, |-,+>, |-,->}. *)
bExplicit = {bX, bY, bZ};
bBarExplicit = {bBarX, bBarY, bBarZ};
correlationExplicit = {
   {cXX, cXY, cXZ},
   {cYX, cYY, cYZ},
   {cZX, cZY, cZZ}
   };

rhoGeneralExplicit = Expand[
   spinDensityMatrix[bExplicit, bBarExplicit, correlationExplicit]
   ];

photonDensityMatrix = spinDensityMatrix[quarkPolarizationVector, antiquarkPolarizationVector, correlationMatrix];

photonDensityElement[row_Integer, column_Integer] /;
   1 <= row <= 4 && 1 <= column <= 4 := photonDensityMatrix[[row, column]];

photonDensityChecks = <|
   "UnitTrace" -> PossibleZeroQ[Together[Tr[photonDensityMatrix] - 1]],
   "Hermitian" -> TrueQ[
     Simplify[
       photonDensityMatrix == ConjugateTranspose[photonDensityMatrix],
       Assumptions -> kinematicAssumptions
       ]
     ]
   |>;

(* The partial traces contain the B vectors:
     Tr_qbar rho = (I+B.sigma)/2,
     Tr_q    rho = (I+BBar.sigma)/2. *)
reducedDensityQ[bVec_?VectorQ] :=
  (id2 + Sum[bVec[[i]] pauli[[i]], {i, 1, 3}])/2;

reducedDensityQBar[bBarVec_?VectorQ] :=
  (id2 + Sum[bBarVec[[j]] pauli[[j]], {j, 1, 3}])/2;

spinDensityCoefficientChecks = <|
   "UnitTrace" -> PossibleZeroQ[Tr[rhoGeneral] - 1],
   "QuarkB" -> And @@ Table[
      PossibleZeroQ[
       Tr[rhoGeneral . KroneckerProduct[pauli[[i]], id2]] - b[i]
       ],
      {i, 1, 3}
      ],
   "AntiquarkBBar" -> And @@ Table[
      PossibleZeroQ[
       Tr[rhoGeneral . KroneckerProduct[id2, pauli[[j]]]] - bBar[j]
       ],
      {j, 1, 3}
      ],
   "Correlations" -> And @@ Flatten@Table[
      PossibleZeroQ[
       Tr[
         rhoGeneral . KroneckerProduct[pauli[[i]], pauli[[j]]]
         ] - c[i, j]
       ],
      {i, 1, 3}, {j, 1, 3}
      ]
   |>;

(* ---------------------------------------------------------------------- *)
(* 3. Dihadron analyzers and the all-C_ij master kernel                    *)
(* ---------------------------------------------------------------------- *)

(* epsilon^(xy)=+1; hence (epsilon.R)_x=R_y and
   (epsilon.R)_y=-R_x. *)
epsilon2 = {{0, 1}, {-1, 0}};

rVec = rT {Cos[phiR], Sin[phiR]};
pHVec = pHT {Cos[phiH], Sin[phiH]};
rVecBar = rTBar {Cos[phiRBar], Sin[phiRBar]};
pHVecBar = pHTBar {Cos[phiHBar], Sin[phiHBar]};

(*
  h1AngleNormQ and h1AngleNormQBar are deliberately symbolic.  For explicit normalization use
  h1AngleNormQ=M1+M2 and h1AngleNormQBar=M1Bar+M2Bar.  If those mass factors have already
  been absorbed into H1Angle, set both normalizations to one.
*)
dihadronNormalizationRules = {
   h1AngleNormQ -> m1 + m2,
   h1AngleNormQBar -> m1Bar + m2Bar
   };

absorbedH1NormalizationRules = {
   h1AngleNormQ -> 1,
   h1AngleNormQBar -> 1
   };

(* The six lower-case FF symbols below stand for FFs evaluated at the
   displayed kinematics.  Apply this optional rule list when explicit
   functional dependence is needed in an integration or a fit. *)
fragmentationFunctionRules = {
   d1Q -> d1QFunction[
     z, zeta, pHT, rT, Cos[phiR - phiH]
     ],
   h1AngleQ -> h1AngleQFunction[
     z, zeta, pHT, rT, Cos[phiR - phiH]
     ],
   g1PerpQ -> g1PerpQFunction[
     z, zeta, pHT, rT, Cos[phiR - phiH]
     ],
   d1QBar -> d1QBarFunction[
     zBar, zetaBar, pHTBar, rTBar, Cos[phiRBar - phiHBar]
     ],
   h1AngleQBar -> h1AngleQBarFunction[
     zBar, zetaBar, pHTBar, rTBar, Cos[phiRBar - phiHBar]
     ],
   g1PerpQBar -> g1PerpQBarFunction[
     zBar, zetaBar, pHTBar, rTBar, Cos[phiRBar - phiHBar]
     ]
   };

transverseAnalyzerQ = (epsilon2 . rVec) h1AngleQ/h1AngleNormQ;

(*
  The antiquark fragmentation correlator is built from the mirror
  light-cone structures of the quark's (it moves in the opposite
  direction), and its spin-dependent (G1Perp-, H1Angle-type) terms carry
  an extra relative minus sign under quark->antiquark that D1 does not --
  a C-conjugation effect (cf. EntangleCollinsFF.pdf, the remark after its
  Eq. (27), and compare its Eqs. (24)/(27) or (35)/(36) directly). A plain
  bar substitution of transverseAnalyzerQ/longitudinalAnalyzerQ is
  therefore not enough; both antiquark analyzers pick up an overall minus
  sign here.
*)
transverseAnalyzerQBar = -(epsilon2 . rVecBar) h1AngleQBar/h1AngleNormQBar;

(* epsilon^(ij) R_i k_j is invariant under the transverse boost to
   k_T=-P_hT/z. *)

longitudinalAnalyzerQ =
  rT pHT Sin[phiR - phiH] g1PerpQ/(z m1 m2);

longitudinalAnalyzerQBar =
  -rTBar pHTBar Sin[phiRBar - phiHBar] g1PerpQBar/
   (zBar m1Bar m2Bar);

analyzerQ = Join[transverseAnalyzerQ, {longitudinalAnalyzerQ}];
analyzerQBar =
  Join[transverseAnalyzerQBar, {longitudinalAnalyzerQBar}];

(*
  The chosen analyzer convention uses +C_ij for TT, TZ and ZT terms and
  -C_zz for the G1Perp product.  The signs are kept in a separate matrix.
*)
correlationSignMatrix = {
   {1, 1, 1},
   {1, 1, 1},
   {1, 1, -1}
   };

(* Single-spin analyzer signs. *)
quarkPolarizationSignVector = {1, 1, 1};
antiquarkPolarizationSignVector = {1, 1, 1};

applyVectorSigns[signs_?VectorQ, vector_?VectorQ] := Table[
   signs[[i]] vector[[i]],
   {i, 1, 3}
   ];

applyCorrelationSigns[signs_?MatrixQ, matrix_?MatrixQ] := Table[
   signs[[i, j]] matrix[[i, j]],
   {i, 1, 3}, {j, 1, 3}
   ];

fragmentationAnalyzerMatrix[d1_, analyzer_?VectorQ] :=
  d1 id2 + Sum[analyzer[[i]] pauli[[i]], {i, 1, 3}];

(*
  Direct Pauli-space realization of the general density matrix:

    Tr[rho (F_q x F_qbar)]
      = D1 D1Bar + D1Bar B_i A_i + D1 BBar_j ABar_j
        + C_ij A_i ABar_j.
*)
densityMatrixKernel[
   bVec_?VectorQ,
   bBarVec_?VectorQ,
   cMat_?MatrixQ
   ] := Expand@Tr[
    spinDensityMatrix[bVec, bBarVec, cMat] .
     KroneckerProduct[
      fragmentationAnalyzerMatrix[d1Q, analyzerQ],
      fragmentationAnalyzerMatrix[d1QBar, analyzerQBar]
      ]
    ];

densityKernelExpected = Expand[
   d1Q d1QBar +
    d1QBar Sum[b[i] analyzerQ[[i]], {i, 1, 3}] +
    d1Q Sum[bBar[j] analyzerQBar[[j]], {j, 1, 3}] +
    Sum[
     c[i, j] analyzerQ[[i]] analyzerQBar[[j]],
     {i, 1, 3}, {j, 1, 3}
     ]
   ];

densityKernelCheck = TrueQ[
   Expand[
     densityMatrixKernel[bSymbolic, bBarSymbolic, correlationSymbolic] -
      densityKernelExpected
     ] == 0
   ];

(* General master kernel, including B and BBar. *)
masterSpinKernel[
   bVec_?VectorQ,
   bBarVec_?VectorQ,
   cMat_?MatrixQ
   ] := densityMatrixKernel[
   applyVectorSigns[quarkPolarizationSignVector, bVec],
   applyVectorSigns[antiquarkPolarizationSignVector, bBarVec],
   applyCorrelationSigns[correlationSignMatrix, cMat]
   ];

(* Backward-compatible correlation-only form. *)
masterSpinKernel[cMat_?MatrixQ] :=
  masterSpinKernel[ConstantArray[0, 3], ConstantArray[0, 3], cMat];

(* Symbolic vectors and matrix expose the coefficient of every parameter. *)
spinDensityMasterGeneral =
  masterSpinKernel[bSymbolic, bBarSymbolic, correlationSymbolic];
correlationOnlyMaster = masterSpinKernel[correlationSymbolic];

quarkPolarizationWeights = Table[
   Coefficient[spinDensityMasterGeneral, b[i]],
   {i, 1, 3}
   ];

antiquarkPolarizationWeights = Table[
   Coefficient[spinDensityMasterGeneral, bBar[j]],
   {j, 1, 3}
   ];

spinDensityWeights = Table[
   Coefficient[spinDensityMasterGeneral, c[i, j]],
   {i, 1, 3}, {j, 1, 3}
   ];

quarkPolarizationWeight[i_Integer] /; 1 <= i <= 3 :=
  quarkPolarizationWeights[[i]];

antiquarkPolarizationWeight[j_Integer] /; 1 <= j <= 3 :=
  antiquarkPolarizationWeights[[j]];

spinDensityWeight[i_Integer, j_Integer] /;
   1 <= i <= 3 && 1 <= j <= 3 := spinDensityWeights[[i, j]];

(* Kernel specialized to the photon-exchange result derived above. *)
photonSpinKernel = masterSpinKernel[quarkPolarizationVector, antiquarkPolarizationVector, correlationMatrix];

hardPrefactor =
  nc Pi alphaEM^2/(2 qSquared) angularFactor chargeQ^2;

spinDensityDifferential[
   bVec_?VectorQ,
   bBarVec_?VectorQ,
   cMat_?MatrixQ
   ] := hardPrefactor masterSpinKernel[bVec, bBarVec, cMat];

spinDensityDifferential[cMat_?MatrixQ] :=
  spinDensityDifferential[ConstantArray[0, 3], ConstantArray[0, 3], cMat];

photonSpinDifferential =
  spinDensityDifferential[quarkPolarizationVector, antiquarkPolarizationVector, correlationMatrix];
photonSpinDifferentialNormalized =
  photonSpinDifferential /. dihadronNormalizationRules;

(* ---------------------------------------------------------------------- *)
(* 4. Master integral and element-by-element projections                  *)
(* ---------------------------------------------------------------------- *)

(* d^2 P_hT d^2 R_T d^2 P'_hT d^2 R'_T in polar coordinates. *)
polarJacobian = pHT rT pHTBar rTBar;

integralAssumptions = kinematicAssumptions &&
   Element[
    {pHT, rT, pHTBar, rTBar, z, zBar, m1, m2, m1Bar, m2Bar},
    Reals
    ] &&
   pHT >= 0 && rT >= 0 && pHTBar >= 0 && rTBar >= 0 &&
   z > 0 && zBar > 0 && m1 > 0 && m2 > 0 && m1Bar > 0 &&
   m2Bar > 0;

integrateOver[expression_, ranges_List] := Apply[
   Integrate,
   Join[
    {expression},
    ranges,
    {Assumptions -> integralAssumptions}
    ]
   ];

(*
  General observable:

    masterIntegral[B, BBar, C, weight, {{var1,min1,max1},...}]

  Examples:

    masterIntegral[quarkPolarizationVector,antiquarkPolarizationVector,correlationMatrix,1,
      {{phiR,0,2 Pi},{phiH,0,2 Pi},
       {phiRBar,0,2 Pi},{phiHBar,0,2 Pi}}]

    elementIntegral[1,3,1,
      {{phiR,0,2 Pi},{phiH,0,2 Pi},
       {phiRBar,0,2 Pi},{phiHBar,0,2 Pi}}]
*)
masterIntegral[
   bVec_?VectorQ,
   bBarVec_?VectorQ,
   cMat_?MatrixQ,
   weight_,
   ranges_List
   ] := integrateOver[
   polarJacobian weight spinDensityDifferential[bVec, bBarVec, cMat],
   ranges
   ];

(* Backward-compatible correlation-only integral. *)
masterIntegral[cMat_?MatrixQ, weight_, ranges_List] :=
  masterIntegral[
   ConstantArray[0, 3], ConstantArray[0, 3], cMat, weight, ranges
   ];

(* Integral of a selected B_i or BBar_j coefficient. *)
quarkPolarizationIntegral[
   i_Integer,
   weight_,
   ranges_List
   ] /; 1 <= i <= 3 := integrateOver[
   polarJacobian weight hardPrefactor quarkPolarizationWeight[i],
   ranges
   ];

antiquarkPolarizationIntegral[
   j_Integer,
   weight_,
   ranges_List
   ] /; 1 <= j <= 3 := integrateOver[
   polarJacobian weight hardPrefactor antiquarkPolarizationWeight[j],
   ranges
   ];

(* Integral of the coefficient of one chosen C_ij.  Calling this for
   i,j=1,2,3 provides the nine spin-density correlation channels. *)
elementIntegral[
   i_Integer,
   j_Integer,
   weight_,
   ranges_List
   ] /; 1 <= i <= 3 && 1 <= j <= 3 :=
  integrateOver[
   polarJacobian weight hardPrefactor spinDensityWeight[i, j],
   ranges
   ];

fullAzimuthRanges = {
   {phiR, 0, 2 Pi},
   {phiH, 0, 2 Pi},
   {phiRBar, 0, 2 Pi},
   {phiHBar, 0, 2 Pi}
   };

(* Compact summary printed when this file is evaluated as a script. *)
Print["Dirac algebra checks: ", diracAlgebraChecks];
Print["Trace derivation gives the compact correlation matrix: ", traceDerivationCheck];
Print["Born photon B=BBar=0 check: ", bornPolarizationCheck];
Print["C_ij in the (x,y,z)=(-r,n,k) basis:"];
Print[MatrixForm[correlationMatrix]];
Print["B_i from the production traces: ", quarkPolarizationVector];
Print["BBar_j from the production traces: ", antiquarkPolarizationVector];
Print["Photon density-matrix checks: ", photonDensityChecks];
Print["General density-matrix coefficient checks: ", spinDensityCoefficientChecks];
Print["Explicit general spin-density matrix rho:"];
Print[MatrixForm[rhoGeneralExplicit]];
Print["Density-matrix contraction check: ", densityKernelCheck];
Print["Quark B_i weights:"];
Print[MatrixForm[quarkPolarizationWeights]];
Print["Antiquark BBar_j weights:"];
Print[MatrixForm[antiquarkPolarizationWeights]];
Print["Weights multiplying all nine C_ij:"];
Print[MatrixForm[spinDensityWeights]];
