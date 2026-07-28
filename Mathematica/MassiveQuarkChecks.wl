(* ::Package:: *)

BeginPackage["MassiveQuarkChecks`"];

BornAssumptions::usage =
  "BornAssumptions[s,m,theta] gives the assumptions used in the symbolic reductions.";

MetricTensor::usage = "MetricTensor[] returns diag(1,-1,-1,-1).";
GammaUpper::usage = "GammaUpper[] returns {gamma^0,gamma^1,gamma^2,gamma^3}.";
GammaLower::usage = "GammaLower[] returns {gamma_0,gamma_1,gamma_2,gamma_3}.";
Gamma5Matrix::usage = "Gamma5Matrix[] returns gamma5 in the explicit Dirac representation.";
MinkowskiDot4::usage = "MinkowskiDot4[a,b] computes a.g.b.";
Slash4::usage = "Slash4[v] computes slash(v)=v^mu gamma_mu.";

BornKinematics::usage =
  "BornKinematics[s,m,theta] returns the exact center-of-mass kinematics as an Association.";
BornKernel::usage =
  "BornKernel[x,y,kin] computes L^{mu nu} Tr[x gamma_mu y gamma_nu].";
ExactBornCalculation::usage =
  "ExactBornCalculation[s,m,theta] explicitly evaluates the unpolarized Born trace.";
SpinCorrelationCalculation::usage =
  "SpinCorrelationCalculation[s,m,theta] explicitly evaluates the (n,r,k) spin-correlation matrix.";
TwistNumeratorCalculation::usage =
  "TwistNumeratorCalculation[s,m,theta] decomposes A_U C into twist-2, twist-3, and twist-4 terms.";

ProjectorCoefficient::usage =
  "ProjectorCoefficient[qName,aqName,theta] evaluates a normalized Born projector coefficient.";
ProjectorCoefficientTable::usage =
  "ProjectorCoefficientTable[theta] evaluates every coefficient quoted in the note.";

FreeBlockCalculation::usage =
  "FreeBlockCalculation[s,m,theta,kr,kn,kbr,kbn,eps] checks the 22+23+32+24+42+33 free-block expansion.";
ExactMassBlockCalculation::usage =
  "ExactMassBlockCalculation[s,m,theta] checks the exact free blocks through order m^2.";

WWNormalizationCalculation::usage =
  "WWNormalizationCalculation[s,m,kr,kn,z,mh] checks the gWW and on-shell closure prefactors.";
HarmonicCalculation::usage =
  "HarmonicCalculation[s,m,theta,phiS,phiSb,phiP,phiPb,pT,pTb,z,zb,mh,mhb] checks the TT and double-Collins harmonics.";

TTAppendixCheck::usage =
  "TTAppendixCheck[] reconstructs the displayed double-transverse appendix terms from projector coefficients.";
CollinsAppendixCheck::usage =
  "CollinsAppendixCheck[] reconstructs the displayed double-Collins appendix terms from projector coefficients.";

CoreVerificationTests::usage =
  "CoreVerificationTests[] returns VerificationTest objects for the central identities.";
RunAllChecks::usage =
  "RunAllChecks[] evaluates TestReport[CoreVerificationTests[]].";

Begin["`Private`"];

ClearAll["`*"];

$metric = DiagonalMatrix[{1, -1, -1, -1}];
$id2 = IdentityMatrix[2];
$id4 = IdentityMatrix[4];
$z2 = ConstantArray[0, {2, 2}];

$pauli = {
  {{0, 1}, {1, 0}},
  {{0, -I}, {I, 0}},
  {{1, 0}, {0, -1}}
};

block4[a_, b_, c_, d_] := ArrayFlatten[{{a, b}, {c, d}}];

$gammaU = {
  block4[$id2, $z2, $z2, -$id2],
  block4[$z2, $pauli[[1]], -$pauli[[1]], $z2],
  block4[$z2, $pauli[[2]], -$pauli[[2]], $z2],
  block4[$z2, $pauli[[3]], -$pauli[[3]], $z2]
};

$gammaD = MapThread[#1 #2 &, {Diagonal[$metric], $gammaU}];
$gamma5 = FullSimplify[I $gammaU[[1]].$gammaU[[2]].$gammaU[[3]].$gammaU[[4]]];

MetricTensor[] := $metric;
GammaUpper[] := $gammaU;
GammaLower[] := $gammaD;
Gamma5Matrix[] := $gamma5;

MinkowskiDot4[a_, b_] := FullSimplify[a.$metric.b];
Slash4[v_] := FullSimplify[Sum[v[[mu]] $gammaD[[mu]], {mu, 1, 4}]];

sigmaMatrix[a_, b_] := I/2 (a.b - b.a);
iSigmaGamma5[a_, b_] := FullSimplify[I sigmaMatrix[a, b].$gamma5];

BornAssumptions[s_, m_, th_] :=
  s > 0 && m > 0 && 4 m^2 < s && Element[th, Reals];

reduce[expr_, ass_: True] :=
  FullSimplify[TrigReduce[expr], Assumptions -> ass];

zeroQ[expr_, ass_: True] := TrueQ[reduce[expr, ass] == 0];

BornKinematics[s_, m_, th_] := Module[
  {rho, beta, e, nPlus, nMinus, p, pb, k, kb, q, L, ass},
  ass = BornAssumptions[s, m, th];
  rho = 4 m^2/s;
  beta = Sqrt[1 - rho];
  e = Sqrt[s]/2;
  nPlus = {1, 0, 0, 1}/Sqrt[2];
  nMinus = {1, 0, 0, -1}/Sqrt[2];
  p = e {1, Sin[th], 0, Cos[th]};
  pb = e {1, -Sin[th], 0, -Cos[th]};
  k = e {1, 0, 0, beta};
  kb = e {1, 0, 0, -beta};
  q = {Sqrt[s], 0, 0, 0};
  L = Table[
    reduce[
      p[[mu]] pb[[nu]] + p[[nu]] pb[[mu]]
        - $metric[[mu, nu]] MinkowskiDot4[p, pb],
      ass
    ],
    {mu, 1, 4}, {nu, 1, 4}
  ];
  <|
    "s" -> s, "m" -> m, "theta" -> th,
    "rho" -> rho, "beta" -> beta, "energy" -> e,
    "nPlus" -> nPlus, "nMinus" -> nMinus,
    "p" -> p, "pbar" -> pb, "k" -> k, "kbar" -> kb, "q" -> q,
    "L" -> L, "Assumptions" -> ass,
    "A0" -> 1 + Cos[th]^2, "B0" -> Sin[th]^2,
    "AU" -> 1 + Cos[th]^2 + rho Sin[th]^2
  |>
];

BornKernel[x_, y_, kin_Association] := reduce[
  Sum[
    kin["L"][[mu, nu]]
      Tr[x.$gammaD[[mu]].y.$gammaD[[nu]]],
    {mu, 1, 4}, {nu, 1, 4}
  ],
  kin["Assumptions"]
];

ExactBornCalculation[s_, m_, th_] := Module[
  {kin, raw, normalized, expected},
  kin = BornKinematics[s, m, th];
  raw = BornKernel[
    Slash4[kin["k"]] + m $id4,
    Slash4[kin["kbar"]] - m $id4,
    kin
  ];
  normalized = reduce[raw/s^2, kin["Assumptions"]];
  expected = kin["AU"];
  <|
    "RawTrace" -> raw,
    "NormalizedTrace" -> normalized,
    "Expected" -> expected,
    "Difference" -> reduce[normalized - expected, kin["Assumptions"]]
  |>
];

SpinCorrelationCalculation[s_, m_, th_] := Module[
  {kin, e, beta, qBasis, aqBasis, den, bilinear, singleQ, singleAQ,
   computed, expected, ass},
  kin = BornKinematics[s, m, th];
  ass = kin["Assumptions"];
  e = kin["energy"];
  beta = kin["beta"];
  qBasis = {
    {0, 0, 1, 0},
    {0, 1, 0, 0},
    {e beta/m, 0, 0, e/m}
  };
  aqBasis = {
    {0, 0, 1, 0},
    {0, 1, 0, 0},
    {-e beta/m, 0, 0, e/m}
  };
  den = BornKernel[
    Slash4[kin["k"]] + m $id4,
    Slash4[kin["kbar"]] - m $id4,
    kin
  ];
  bilinear[sv_, svb_] := BornKernel[
    (Slash4[kin["k"]] + m $id4).$gamma5.Slash4[sv],
    (Slash4[kin["kbar"]] - m $id4).$gamma5.Slash4[svb],
    kin
  ];
  singleQ[sv_] := BornKernel[
    (Slash4[kin["k"]] + m $id4).$gamma5.Slash4[sv],
    Slash4[kin["kbar"]] - m $id4,
    kin
  ];
  singleAQ[svb_] := BornKernel[
    Slash4[kin["k"]] + m $id4,
    (Slash4[kin["kbar"]] - m $id4).$gamma5.Slash4[svb],
    kin
  ];
  computed = Table[
    reduce[bilinear[qBasis[[a]], aqBasis[[b]]]/den, ass],
    {a, 1, 3}, {b, 1, 3}
  ];
  expected = 1/kin["AU"] {
    {-(1 - kin["rho"]) kin["B0"], 0, 0},
    {0, (1 + kin["rho"]) kin["B0"],
      2 Sqrt[kin["rho"]] Sin[th] Cos[th]},
    {0, 2 Sqrt[kin["rho"]] Sin[th] Cos[th],
      kin["A0"] - kin["rho"] kin["B0"]}
  };
  <|
    "BasisOrder" -> {"n", "r", "k"},
    "Computed" -> computed,
    "Expected" -> expected,
    "Difference" -> reduce[computed - expected, ass],
    "SingleQuarkSpinKernels" -> (reduce[singleQ[#], ass] & /@ qBasis),
    "SingleAntiquarkSpinKernels" -> (reduce[singleAQ[#], ass] & /@ aqBasis)
  |>
];

TwistNumeratorCalculation[s_, m_, th_] := Module[
  {kin, spin, t2, t3, t4, exact, ass},
  kin = BornKinematics[s, m, th];
  ass = kin["Assumptions"];
  spin = SpinCorrelationCalculation[s, m, th];
  exact = reduce[kin["AU"] spin["Computed"], ass];
  t2 = {
    {-kin["B0"], 0, 0},
    {0, kin["B0"], 0},
    {0, 0, kin["A0"]}
  };
  t3 = Sqrt[kin["rho"]] Sin[2 th] {
    {0, 0, 0},
    {0, 0, 1},
    {0, 1, 0}
  };
  t4 = kin["rho"] kin["B0"] {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, -1}
  };
  <|
    "ExactNumerator" -> exact,
    "Twist2" -> t2,
    "Twist3" -> t3,
    "Twist4" -> t4,
    "Difference" -> reduce[exact - t2 - t3 - t4, ass]
  |>
];

$nPlus = {1, 0, 0, 1}/Sqrt[2];
$nMinus = {1, 0, 0, -1}/Sqrt[2];
$gMinus = Slash4[$nPlus];
$gPlus = Slash4[$nMinus];
$gR = $gammaU[[2]];
$gN = $gammaU[[3]];
$vR = -$gR;
$vN = -$gN;

$qNames = {
  "V+", "A+", "Tr+", "Tn+",
  "S", "P", "Vr", "Vn", "Ar", "An", "Trn", "T-+",
  "V-", "A-", "Tr-", "Tn-"
};

$aqNames = {
  "V-", "A-", "Tr-", "Tn-",
  "S", "P", "Vr", "Vn", "Ar", "An", "Trn", "T+-",
  "V+", "A+", "Tr+", "Tn+"
};

$qGamma = <|
  "V+" -> $gMinus,
  "A+" -> $gMinus.$gamma5,
  "Tr+" -> iSigmaGamma5[$gR, $gPlus],
  "Tn+" -> iSigmaGamma5[$gN, $gPlus],
  "S" -> $id4,
  "P" -> I $gamma5,
  "Vr" -> $vR,
  "Vn" -> $vN,
  "Ar" -> $vR.$gamma5,
  "An" -> $vN.$gamma5,
  "Trn" -> iSigmaGamma5[$gR, $gN],
  "T-+" -> iSigmaGamma5[$gMinus, $gPlus],
  "V-" -> $gPlus,
  "A-" -> $gPlus.$gamma5,
  "Tr-" -> iSigmaGamma5[$gR, $gMinus],
  "Tn-" -> iSigmaGamma5[$gN, $gMinus]
|>;

$aqGamma = <|
  "V-" -> $gPlus,
  "A-" -> $gPlus.$gamma5,
  "Tr-" -> iSigmaGamma5[$gR, $gMinus],
  "Tn-" -> iSigmaGamma5[$gN, $gMinus],
  "S" -> $id4,
  "P" -> I $gamma5,
  "Vr" -> $vR,
  "Vn" -> $vN,
  "Ar" -> $vR.$gamma5,
  "An" -> $vN.$gamma5,
  "Trn" -> iSigmaGamma5[$gR, $gN],
  "T+-" -> iSigmaGamma5[$gPlus, $gMinus],
  "V+" -> $gMinus,
  "A+" -> $gMinus.$gamma5,
  "Tr+" -> iSigmaGamma5[$gR, $gPlus],
  "Tn+" -> iSigmaGamma5[$gN, $gPlus]
|>;

makeDualBasis[names_, gammaAssociation_] := Module[
  {mats, gram, inverseGram},
  mats = Lookup[gammaAssociation, names];
  gram = Table[
    FullSimplify[Tr[mats[[i]].mats[[j]]]/2],
    {i, Length[names]}, {j, Length[names]}
  ];
  inverseGram = Inverse[gram];
  AssociationThread[
    names,
    Table[
      FullSimplify[
        Sum[inverseGram[[i, j]] mats[[j]], {j, Length[names]}]
      ],
      {i, Length[names]}
    ]
  ]
];

$qDual = makeDualBasis[$qNames, $qGamma];
$aqDual = makeDualBasis[$aqNames, $aqGamma];

ProjectorCoefficient[qName_String, aqName_String, th_] := Module[
  {p, pb, L},
  p = {1, Sin[th], 0, Cos[th]}/2;
  pb = {1, -Sin[th], 0, -Cos[th]}/2;
  L = Table[
    FullSimplify[
      p[[mu]] pb[[nu]] + p[[nu]] pb[[mu]]
        - $metric[[mu, nu]] MinkowskiDot4[p, pb]
    ],
    {mu, 1, 4}, {nu, 1, 4}
  ];
  reduce[
    2 Sum[
      L[[mu, nu]]
        Tr[$qDual[qName].$gammaD[[mu]].$aqDual[aqName].$gammaD[[nu]]],
      {mu, 1, 4}, {nu, 1, 4}
    ],
    Element[th, Reals]
  ]
];

coefficientTargets[th_] := {
  {"22: V+ Vbar-", "V+", "V-", 1 + Cos[th]^2},
  {"22: A+ Abar-", "A+", "A-", 1 + Cos[th]^2},
  {"22: Tr+ Trbar-", "Tr+", "Tr-", Sin[th]^2},
  {"22: Tn+ Tnbar-", "Tn+", "Tn-", -Sin[th]^2},

  {"23: V+ Vbar-r", "V+", "Vr", -Sqrt[2] Sin[th] Cos[th]},
  {"23: A+ Abar-r", "A+", "Ar", -Sqrt[2] Sin[th] Cos[th]},
  {"23: Tr+ Tbar+-", "Tr+", "T+-", Sqrt[2] Sin[th] Cos[th]},
  {"23: Tn+ Tbar-rn", "Tn+", "Trn", Sqrt[2] Sin[th] Cos[th]},

  {"32: Vr Vbar-", "Vr", "V-", Sqrt[2] Sin[th] Cos[th]},
  {"32: Ar Abar-", "Ar", "A-", Sqrt[2] Sin[th] Cos[th]},
  {"32: Trn Tnbar-", "Trn", "Tn-", -Sqrt[2] Sin[th] Cos[th]},
  {"32: T-+ Trbar-", "T-+", "Tr-", -Sqrt[2] Sin[th] Cos[th]},

  {"24: V+ Vbar+", "V+", "V+", Sin[th]^2},
  {"24: A+ Abar+", "A+", "A+", Sin[th]^2},
  {"24: Tr+ Trbar+", "Tr+", "Tr+", Sin[th]^2},
  {"24: Tn+ Tnbar+", "Tn+", "Tn+", Sin[th]^2},

  {"42: V- Vbar-", "V-", "V-", Sin[th]^2},
  {"42: A- Abar-", "A-", "A-", Sin[th]^2},
  {"42: Tr- Trbar-", "Tr-", "Tr-", Sin[th]^2},
  {"42: Tn- Tnbar-", "Tn-", "Tn-", Sin[th]^2},

  {"33: S Sbar", "S", "S", -2},
  {"33: P Pbar", "P", "P", -2},
  {"33: Vr Vbar-r", "Vr", "Vr", -2 Sin[th]^2},
  {"33: Ar Abar-r", "Ar", "Ar", -2 Sin[th]^2},
  {"33: Trn Tbar-rn", "Trn", "Trn", 2 Cos[th]^2},
  {"33: T-+ Tbar+-", "T-+", "T+-", -2 Cos[th]^2}
};

ProjectorCoefficientTable[th_] := Module[
  {rows},
  rows = coefficientTargets[th];
  Prepend[
    Map[
      Function[row,
        With[{value = ProjectorCoefficient[row[[2]], row[[3]], th]},
          {
            row[[1]], value, row[[4]],
            reduce[value - row[[4]], Element[th, Reals]]
          }
        ]
      ],
      rows
    ],
    {"Coefficient", "Computed", "Expected", "Difference"}
  ]
];

FreeBlockCalculation[
  s_, m_, th_, kr_, kn_, kbr_, kbn_, eps_
] := Module[
  {kin, kLarge, kbLarge, kPerp, kbPerp, kSmall, kbSmall,
   kLC, kbLC, qBlocks, aqBlocks, pairs, full, blocks, scaledFull,
   scaledBlocks, fullSeries, blockSeries, order0, order1, order2,
   target1, target2, ass},
  kin = BornKinematics[s, m, th];
  ass = kin["Assumptions"] &&
    Element[{kr, kn, kbr, kbn, eps}, Reals];
  kLarge = Sqrt[s]/Sqrt[2];
  kbLarge = Sqrt[s]/Sqrt[2];
  kPerp = {0, kr, kn, 0};
  kbPerp = {0, kbr, kbn, 0};
  kSmall = (m^2 + kr^2 + kn^2)/(2 kLarge);
  kbSmall = (m^2 + kbr^2 + kbn^2)/(2 kbLarge);
  kLC = kLarge $nPlus + kSmall $nMinus + kPerp;
  kbLC = kbSmall $nPlus + kbLarge $nMinus + kbPerp;
  qBlocks = <|
    2 -> kLarge $gMinus,
    3 -> Slash4[kPerp] + m $id4,
    4 -> kSmall $gPlus
  |>;
  aqBlocks = <|
    2 -> kbLarge $gPlus,
    3 -> Slash4[kbPerp] - m $id4,
    4 -> kbSmall $gMinus
  |>;
  pairs = {{2, 2}, {2, 3}, {3, 2}, {2, 4}, {4, 2}, {3, 3}};
  full = BornKernel[Slash4[kLC] + m $id4, Slash4[kbLC] - m $id4, kin]/s^2;
  blocks = Total[
    BornKernel[qBlocks[#[[1]]], aqBlocks[#[[2]]], kin]/s^2 & /@ pairs
  ];
  scaledFull = full /. {
    m -> eps m, kr -> eps kr, kn -> eps kn,
    kbr -> eps kbr, kbn -> eps kbn
  };
  scaledBlocks = blocks /. {
    m -> eps m, kr -> eps kr, kn -> eps kn,
    kbr -> eps kbr, kbn -> eps kbn
  };
  fullSeries = reduce[
    Expand[Normal[Series[scaledFull, {eps, 0, 2}]]], ass
  ];
  blockSeries = reduce[
    Expand[Normal[Series[scaledBlocks, {eps, 0, 2}]]], ass
  ];
  order0 = reduce[Coefficient[fullSeries, eps, 0], ass];
  order1 = reduce[Coefficient[fullSeries, eps, 1], ass];
  order2 = reduce[Coefficient[fullSeries, eps, 2], ass];
  target1 = Sin[2 th] (kr - kbr)/Sqrt[s];
  target2 = Sin[th]^2 (
    kr^2 + kn^2 + kbr^2 + kbn^2 - 4 kr kbr
  )/s;
  <|
    "FullSeries" -> fullSeries,
    "BlockSeries" -> blockSeries,
    "BlockDifference" -> reduce[fullSeries - blockSeries, ass],
    "Order0" -> order0,
    "Order1" -> order1,
    "Order1Expected" -> target1,
    "Order1Difference" -> reduce[order1 - target1, ass],
    "Order2AtZeroMass" -> reduce[order2 /. m -> 0, ass],
    "Order2KTExpected" -> target2,
    "Order2KTDifference" -> reduce[(order2 /. m -> 0) - target2, ass],
    "LinearMassAtOrder2" -> reduce[Coefficient[order2, m, 1], ass]
  |>
];

ExactMassBlockCalculation[s_, m_, th_] := Module[
  {kin, kp, km, kbp, kbm, qBlocks, aqBlocks, pairs, blockKernel,
   exactKernel, ass, exactSeries, blockSeries},
  kin = BornKinematics[s, m, th];
  ass = kin["Assumptions"];
  kp = MinkowskiDot4[kin["k"], kin["nMinus"]];
  km = MinkowskiDot4[kin["k"], kin["nPlus"]];
  kbm = MinkowskiDot4[kin["kbar"], kin["nPlus"]];
  kbp = MinkowskiDot4[kin["kbar"], kin["nMinus"]];
  qBlocks = <|
    2 -> kp $gMinus,
    3 -> m $id4,
    4 -> km $gPlus
  |>;
  aqBlocks = <|
    2 -> kbm $gPlus,
    3 -> -m $id4,
    4 -> kbp $gMinus
  |>;
  pairs = {{2, 2}, {2, 3}, {3, 2}, {2, 4}, {4, 2}, {3, 3}};
  blockKernel = Total[
    BornKernel[qBlocks[#[[1]]], aqBlocks[#[[2]]], kin]/s^2 & /@ pairs
  ];
  exactKernel = ExactBornCalculation[s, m, th]["NormalizedTrace"];
  exactSeries = reduce[Normal[Series[exactKernel, {m, 0, 2}]], ass];
  blockSeries = reduce[Normal[Series[blockKernel, {m, 0, 2}]], ass];
  <|
    "ExactTraceThroughM2" -> exactSeries,
    "RetainedBlocksThroughM2" -> blockSeries,
    "Difference" -> reduce[exactSeries - blockSeries, ass]
  |>
];

WWNormalizationCalculation[s_, m_, kr_, kn_, z_, mh_] := Module[
  {pPlus, scalar, vectorR, small},
  pPlus = z Sqrt[s]/Sqrt[2];
  scalar = reduce[(mh/pPlus) (z m/mh),
    s > 0 && z > 0 && mh > 0];
  vectorR = reduce[(mh/pPlus) (kr/mh) z,
    s > 0 && z > 0 && mh > 0];
  small = reduce[
    (mh^2/pPlus^2)
      z^2 (m^2 + kr^2 + kn^2)/(2 mh^2),
    s > 0 && z > 0 && mh > 0
  ];
  <|
    "PhiSOverD1" -> scalar,
    "PhiSExpected" -> Sqrt[2] m/Sqrt[s],
    "PhiVrOverD1" -> vectorR,
    "PhiVrExpected" -> Sqrt[2] kr/Sqrt[s],
    "PhiVMinusOverD1" -> small,
    "PhiVMinusExpected" -> (m^2 + kr^2 + kn^2)/s,
    "Differences" -> {
      reduce[scalar - Sqrt[2] m/Sqrt[s]],
      reduce[vectorR - Sqrt[2] kr/Sqrt[s]],
      reduce[small - (m^2 + kr^2 + kn^2)/s]
    }
  |>
];

HarmonicCalculation[
  s_, m_, th_, phiS_, phiSb_, phiP_, phiPb_,
  pT_, pTb_, z_, zb_, mh_, mhb_
] := Module[
  {rho, b0, sr, sn, srb, snb, ttFrom, ttExpected,
   cq, caq, ccFrom, ccExpected, ass},
  ass = s > 0 && Element[
    {m, th, phiS, phiSb, phiP, phiPb, pT, pTb, z, zb, mh, mhb},
    Reals
  ] && z > 0 && zb > 0 && mh > 0 && mhb > 0;
  rho = 4 m^2/s;
  b0 = Sin[th]^2;
  sr = Cos[phiS];
  sn = Sin[phiS];
  srb = Cos[phiSb];
  snb = Sin[phiSb];
  ttFrom = b0 ((1 + rho) sr srb - (1 - rho) sn snb);
  ttExpected = b0 (
    Cos[phiS + phiSb] + rho Cos[phiS - phiSb]
  );
  cq = {
    -pT Sin[phiP]/(z mh),
     pT Cos[phiP]/(z mh)
  };
  caq = {
     pTb Sin[phiPb]/(zb mhb),
    -pTb Cos[phiPb]/(zb mhb)
  };
  ccFrom = b0 (
    (1 + rho) cq[[1]] caq[[1]]
      - (1 - rho) cq[[2]] caq[[2]]
  );
  ccExpected = b0 pT pTb/(z zb mh mhb) (
    Cos[phiP + phiPb] - rho Cos[phiP - phiPb]
  );
  <|
    "TTFromTensor" -> reduce[ttFrom, ass],
    "TTExpected" -> ttExpected,
    "TTDifference" -> reduce[ttFrom - ttExpected, ass],
    "CollinsFromTensor" -> reduce[ccFrom, ass],
    "CollinsExpected" -> ccExpected,
    "CollinsDifference" -> reduce[ccFrom - ccExpected, ass]
  |>
];

TTAppendixCheck[] := Module[
  {s, m, th, kr, kn, kbr, kbn, z, zb, mh, mhb, pQ, pA,
   sr, sn, srb, snb, rankR, rankN, rankBR, rankBN,
   h1, h1b, hp, hpb, g, gb, kt2, kbt2, kdots, kbdots,
   kcross, kbcross, hR, hN, hBR, hBN, gR, gBR, cmix,
   qTrn, qTLC, aTrn, aTLC, from3, target3, zeroRules,
   residual, qAr, aAr, fromA, targetA, fromT, targetT,
   aSmallR, aSmallN, qSmallR, qSmallN, fromSmall, targetSmall,
   ass},
  ass = s > 0 && z > 0 && zb > 0 && mh > 0 && mhb > 0 &&
    Element[
      {m, th, kr, kn, kbr, kbn, sr, sn, srb, snb,
       rankR, rankN, rankBR, rankBN, h1, h1b, hp, hpb, g, gb},
      Reals
    ];
  pQ = z Sqrt[s]/Sqrt[2];
  pA = zb Sqrt[s]/Sqrt[2];
  kt2 = kr^2 + kn^2;
  kbt2 = kbr^2 + kbn^2;
  kdots = kr sr + kn sn;
  kbdots = kbr srb + kbn snb;
  kcross = sr kn - kr sn;
  kbcross = srb kbn - kbr snb;
  hR = sr h1 - rankR hp/mh^2;
  hN = sn h1 - rankN hp/mh^2;
  hBR = srb h1b - rankBR hpb/mhb^2;
  hBN = snb h1b - rankBN hpb/mhb^2;
  gR = sr kt2 g/(2 mh^2) - rankR g/mh^2;
  gBR = srb kbt2 gb/(2 mhb^2) - rankBR gb/mhb^2;
  cmix = Sqrt[2] Sin[th] Cos[th];
  qTrn = (mh/pQ) z (kcross/mh) (h1 - kt2 hp/(2 mh^2));
  qTLC = (mh/pQ) z (kdots/mh)
    (h1 + kt2 hp/(2 mh^2) - m g/mh);
  aTrn = (mhb/pA) zb (kbcross/mhb)
    (h1b - kbt2 hpb/(2 mhb^2));
  aTLC = (mhb/pA) zb (kbdots/mhb)
    (h1b + kbt2 hpb/(2 mhb^2) - m gb/mhb);
  from3 = cmix (hR aTLC + hN aTrn) - cmix (qTrn hBN + qTLC hBR);
  target3 = cmix (mhb/pA) zb (
    -(kbdots/mhb) hR (-h1b - kbt2 hpb/(2 mhb^2) + m gb/mhb)
    +(kbcross/mhb) hN (h1b - kbt2 hpb/(2 mhb^2))
  ) - cmix (mh/pQ) z (
    (kcross/mh) (h1 - kt2 hp/(2 mh^2)) hBN
    -(kdots/mh) (-h1 - kt2 hp/(2 mh^2) + m g/mh) hBR
  );
  zeroRules = {
    kr -> 0, kn -> 0, kbr -> 0, kbn -> 0,
    rankR -> 0, rankN -> 0, rankBR -> 0, rankBN -> 0,
    hp -> 0, hpb -> 0, g -> 0, gb -> 0
  };
  residual[x_] := Expand[x - (x /. zeroRules)];
  qAr = (mh/pQ) z (gR + (m/mh) hR);
  aAr = (mhb/pA) zb (gBR + (m/mhb) hBR);
  fromA = residual[-2 Sin[th]^2 qAr aAr];
  targetA = -2 Sin[th]^2 (mh mhb/(pQ pA)) z zb (
    gR gBR + (m/mhb) gR hBR + (m/mh) hR gBR
      + (m^2/(mh mhb)) residual[hR hBR]
  );
  fromT = 2 Cos[th]^2 (qTrn aTrn - qTLC aTLC);
  targetT = 2 Cos[th]^2 (mh mhb/(pQ pA)) z zb/(mh mhb) (
    kcross kbcross
      (h1 - kt2 hp/(2 mh^2))
      (h1b - kbt2 hpb/(2 mhb^2))
    - kdots kbdots
      (h1 + kt2 hp/(2 mh^2) - m g/mh)
      (h1b + kbt2 hpb/(2 mhb^2) - m gb/mhb)
  );
  aSmallR = zb^2 (m^2 + kbt2)/(2 pA^2) srb h1b;
  aSmallN = zb^2 (m^2 + kbt2)/(2 pA^2) snb h1b;
  qSmallR = z^2 (m^2 + kt2)/(2 pQ^2) sr h1;
  qSmallN = z^2 (m^2 + kt2)/(2 pQ^2) sn h1;
  fromSmall =
    Sin[th]^2 residual[hR aSmallR + hN aSmallN]
    + Sin[th]^2 residual[qSmallR hBR + qSmallN hBN];
  targetSmall =
    Sin[th]^2 zb^2/(2 pA^2)
      residual[(m^2 + kbt2) (srb hR + snb hN) h1b]
    + Sin[th]^2 z^2/(2 pQ^2)
      residual[(m^2 + kt2) h1 (sr hBR + sn hBN)];
  <|
    "Twist3Difference" -> reduce[from3 - target3, ass],
    "Axial33Difference" -> reduce[fromA - targetA, ass],
    "Tensor33Difference" -> reduce[fromT - targetT, ass],
    "SmallComponentDifference" -> reduce[fromSmall - targetSmall, ass]
  |>
];

CollinsAppendixCheck[] := Module[
  {s, th, kr, kn, kbr, kbn, z, zb, mh, mhb, pQ, pA,
   h, hb, kt2, kbt2, qR, qN, aR, aN, qTrn, aTrn, cmix,
   leadFrom, leadTarget, t3From, t3Target,
   aSmallR, aSmallN, qSmallR, qSmallN, t4From, t4Target, ass},
  ass = s > 0 && z > 0 && zb > 0 && mh > 0 && mhb > 0 &&
    Element[{th, kr, kn, kbr, kbn, h, hb}, Reals];
  pQ = z Sqrt[s]/Sqrt[2];
  pA = zb Sqrt[s]/Sqrt[2];
  kt2 = kr^2 + kn^2;
  kbt2 = kbr^2 + kbn^2;
  qR = kn h/mh;
  qN = -kr h/mh;
  aR = -kbn hb/mhb;
  aN = kbr hb/mhb;
  cmix = Sqrt[2] Sin[th] Cos[th];
  leadFrom = Sin[th]^2 (qR aR - qN aN);
  leadTarget = Sin[th]^2 (kr kbr - kn kbn) h hb/(mh mhb);
  qTrn = z kt2 h/(pQ mh);
  aTrn = zb kbt2 hb/(pA mhb);
  t3From = cmix qN aTrn - cmix qTrn aN;
  t3Target = -cmix (
    (mhb/pA) zb kr kbt2/(mh mhb^2)
    + (mh/pQ) z kt2 kbr/(mh^2 mhb)
  ) h hb;
  aSmallR = -zb^2 kbt2 kbn hb/(2 pA^2 mhb);
  aSmallN = zb^2 kbt2 kbr hb/(2 pA^2 mhb);
  qSmallR = z^2 kt2 kn h/(2 pQ^2 mh);
  qSmallN = -z^2 kt2 kr h/(2 pQ^2 mh);
  t4From = 2 Cos[th]^2 qTrn aTrn
    + Sin[th]^2 (qR aSmallR + qN aSmallN)
    + Sin[th]^2 (qSmallR aR + qSmallN aN);
  t4Target = (
    2 Cos[th]^2 (mh mhb/(pQ pA)) z zb kt2 kbt2/(mh^2 mhb^2)
    - Sin[th]^2 zb^2 kbt2/(2 pA^2)
      (kr kbr + kn kbn)/(mh mhb)
    - Sin[th]^2 z^2 kt2/(2 pQ^2)
      (kr kbr + kn kbn)/(mh mhb)
  ) h hb;
  <|
    "LeadingDifference" -> reduce[leadFrom - leadTarget, ass],
    "Twist3Difference" -> reduce[t3From - t3Target, ass],
    "Twist4Difference" -> reduce[t4From - t4Target, ass]
  |>
];

CoreVerificationTests[] := Module[
  {s, m, th, kr, kn, kbr, kbn, eps, z, mh,
   born, spin, twist, coeff, blocks, exactBlocks, ww, harm, tt, cc,
   clifford, gamma5Anti, ass},
  ass = BornAssumptions[s, m, th];
  born = ExactBornCalculation[s, m, th];
  spin = SpinCorrelationCalculation[s, m, th];
  twist = TwistNumeratorCalculation[s, m, th];
  coeff = Rest[ProjectorCoefficientTable[th]];
  blocks = FreeBlockCalculation[s, m, th, kr, kn, kbr, kbn, eps];
  exactBlocks = ExactMassBlockCalculation[s, m, th];
  ww = WWNormalizationCalculation[s, m, kr, kn, z, mh];
  harm = HarmonicCalculation[
    s, m, th, phiS, phiSb, phiP, phiPb, pT, pTb, z, zb, mh, mhb
  ];
  tt = TTAppendixCheck[];
  cc = CollinsAppendixCheck[];
  clifford = Table[
    reduce[
      $gammaU[[mu]].$gammaU[[nu]]
        + $gammaU[[nu]].$gammaU[[mu]]
        - 2 $metric[[mu, nu]] $id4
    ],
    {mu, 1, 4}, {nu, 1, 4}
  ];
  gamma5Anti = Table[
    reduce[$gamma5.$gammaU[[mu]] + $gammaU[[mu]].$gamma5],
    {mu, 1, 4}
  ];
  {
    VerificationTest[
      And @@ (TrueQ[# == ConstantArray[0, {4, 4}]] & /@ Flatten[clifford, 1]),
      True, TestID -> "Clifford algebra"
    ],
    VerificationTest[
      TrueQ[reduce[$gamma5.$gamma5 - $id4] == ConstantArray[0, {4, 4}]],
      True, TestID -> "gamma5 squared"
    ],
    VerificationTest[
      And @@ (TrueQ[# == ConstantArray[0, {4, 4}]] & /@ gamma5Anti),
      True, TestID -> "gamma5 anticommutation"
    ],
    VerificationTest[
      zeroQ[born["Difference"], ass],
      True, TestID -> "exact unpolarized Born trace"
    ],
    VerificationTest[
      TrueQ[spin["Difference"] == ConstantArray[0, {3, 3}]],
      True, TestID -> "spin-correlation matrix"
    ],
    VerificationTest[
      TrueQ[twist["Difference"] == ConstantArray[0, {3, 3}]],
      True, TestID -> "twist numerator decomposition"
    ],
    VerificationTest[
      And @@ (zeroQ[#[[4]], Element[th, Reals]] & /@ coeff),
      True, TestID -> "projector coefficients"
    ],
    VerificationTest[
      zeroQ[blocks["BlockDifference"]] &&
      zeroQ[blocks["Order1Difference"]] &&
      zeroQ[blocks["Order2KTDifference"]] &&
      zeroQ[blocks["LinearMassAtOrder2"]],
      True, TestID -> "free transverse block expansion"
    ],
    VerificationTest[
      zeroQ[exactBlocks["Difference"], ass],
      True, TestID -> "exact free mass blocks through m^2"
    ],
    VerificationTest[
      And @@ (zeroQ /@ ww["Differences"]),
      True, TestID -> "gWW normalizations"
    ],
    VerificationTest[
      zeroQ[harm["TTDifference"]] && zeroQ[harm["CollinsDifference"]],
      True, TestID -> "TT and Collins harmonics"
    ],
    VerificationTest[
      And @@ (zeroQ /@ Values[tt]),
      True, TestID -> "double-transverse appendix"
    ],
    VerificationTest[
      And @@ (zeroQ /@ Values[cc]),
      True, TestID -> "double-Collins appendix"
    ]
  }
];

RunAllChecks[] := TestReport[CoreVerificationTests[]];

End[];
EndPackage[];
