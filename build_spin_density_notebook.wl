(* Build the interactive notebook interface for the spin-density code. *)

ClearAll[inputCell];

inputCell[code_String, options___] := Module[{held},
  held = ToExpression["(" <> code <> ")", InputForm, HoldComplete];
  Cell[
   BoxData[
    held /. HoldComplete[expression_] :>
      MakeBoxes[expression, StandardForm]
    ],
   "Input",
   options
   ]
  ];

sourceDirectory = DirectoryName[ExpandFileName[$InputFileName]];
targetNotebook = FileNameJoin[{
    sourceDirectory,
    "spin_density_master.nb"
    }];

notebook = Notebook[{
    Cell[
     "Spin-density master integral",
     "Title"
     ],

    Cell[
     "This notebook loads the separate DiracAlgebra` package, verifies the Dirac algebra, evaluates the production traces, constructs the general spin-density matrix with B, BBar and C, and builds the dihadron master integral.",
     "Text"
     ],

    Cell["Conventions and notation", "Section"],

    Cell[
     "The metric is diag(1,-1,-1,-1), gamma matrices use the Weyl representation, and the common physical spin basis is (x,y,z)=(-r,n,k). The Bar suffix labels the antiquark or opposite jet and is not complex conjugation. The code variables quarkPolarizationVector, antiquarkPolarizationVector and correlationMatrix represent B_a, BBar_b and C_ab, respectively.",
     "Text"
     ],

    Cell["Load and verify the Dirac algebra", "Section"],

    inputCell[
     "notebookDirectory = NotebookDirectory[];\nGet[FileNameJoin[{notebookDirectory, \"DiracAlgebra.wl\"}]];\ndiracAlgebraChecks = AssertDiracAlgebra[];\nIf[diracAlgebraChecks === $Failed, Abort[]];\ndiracAlgebraChecks",
     InitializationCell -> True
     ],

    Cell[
     "The returned association explicitly checks the Clifford anticommutator, the Pauli anticommutator, the definition and Weyl representation of gamma5, gamma5 squared, its anticommutator with every gamma matrix, Hermiticity and traces.",
     "Text"
     ],

    Cell["Evaluate the complete derivation", "Section"],

    inputCell[
     "Get[FileNameJoin[{NotebookDirectory[], \"spin_density_master.wl\"}]];",
     InitializationCell -> True
     ],

    Cell["Explicit general spin-density matrix", "Section"],

    Cell[
     "The basis is {|+,+>, |+,->, |-,+>, |-,->}. The matrix contains all components of B, BBar and C explicitly.",
     "Text"
     ],

    inputCell["rhoGeneralExplicit // MatrixForm"],

    Cell["Density-matrix consistency checks", "Subsection"],

    inputCell["spinDensityCoefficientChecks"],

    inputCell[
     "{reducedDensityQ[bExplicit] // MatrixForm, reducedDensityQBar[bBarExplicit] // MatrixForm}"
     ],

    Cell["Production result", "Section"],

    inputCell["quarkPolarizationVector"],
    inputCell["antiquarkPolarizationVector"],
    inputCell["correlationMatrix // MatrixForm"],
    inputCell["photonDensityMatrix // MatrixForm"],

    Cell["General production vertices", "Section"],

    Cell[
     "For a different production process, supply its initial-state tensor and the vertex matrices on the amplitude and conjugate-amplitude sides.",
     "Text"
     ],

    inputCell[
     "generalParameters = spinParametersFromProduction[initialTensor, leftVertices, rightVertices, Identity];\n{generalParameters[\"B\"], generalParameters[\"BBar\"], generalParameters[\"C\"]}"
     ],

    Cell["General spin-density kernel", "Section"],

    inputCell["spinDensityMasterGeneral"],
    inputCell["quarkPolarizationWeights // MatrixForm"],
    inputCell["antiquarkPolarizationWeights // MatrixForm"],
    inputCell["spinDensityWeights // MatrixForm"],

    Cell["Master integration", "Section"],

    inputCell[
     "masterIntegral[quarkPolarizationVector, antiquarkPolarizationVector, correlationMatrix, 1, fullAzimuthRanges]"
     ],

    Cell[
     "Use quarkPolarizationIntegral, antiquarkPolarizationIntegral and elementIntegral to project individual B_i, BBar_j and C_ij channels.",
     "Text"
     ]
    },
   WindowTitle -> "Spin-density master integral",
   StyleDefinitions -> "Default.nb",
   Saveable -> True
   ];

Put[notebook, targetNotebook];
Print[targetNotebook];
