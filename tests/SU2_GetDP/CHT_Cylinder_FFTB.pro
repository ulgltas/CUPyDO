//

// Group
Group{
  Inner = Region[101];
  Outer = Region[102];
  Omega = Region[103];
  AllOmega = Region[{Inner, Outer, Omega}];
}

// Function
Function{
  DefineConstant[
    Initialize = 1, // initialize (no provided forces and displacements)?
    OutputFiles = 0, // output data files?
    T1 = 0.0, // initial time
    T2 = 0.1 // final time
  ];

  T_in[] = 350.0;

  lambda = 46.0; //W/(mK)

  If(Initialize)
    f[] = Vector[0, 0, 0];
    uNm2[] = Vector[0, 0, 0];
    uNm1[] = Vector[0, 0, 0];
    Q[] = Vector[0,0,0];
    Qn[] = 0.0;
  Else
    If(Exists(nodalForce) &&
       Exists(nodalDisplacementNm2) &&
       Exists(nodalDisplacementNm1) &&
       Exists(nodalTemperature) &&
       Exists(nodalHeatFlux))
      Printf("Received nodal forces and previous displacements from API");
    Else
      nodalForce() = ListFromFile["nodalForce.txt"];
      nodalDisplacementNm2() = ListFromFile["nodalDisplacementNm2.txt"];
      nodalDisplacementNm1() = ListFromFile["nodalDisplacementNm1.txt"];
      nodalTemperature() = ListFromFile["nodalTemperature.txt"];
      nodalHeatFlux() = ListFromFile["nodalHeatFlux.txt"];
      nodalNormalHeatFlux() = ListFromFile["nodalNormalHeatFlux.txt"];
    EndIf
    f[] = VectorFromIndex[]{nodalForce()};
    uNm2[] = VectorFromIndex[]{nodalDisplacementNm2()};
    uNm1[] = VectorFromIndex[]{nodalDisplacementNm1()};
    Q[] = VectorFromIndex[]{nodalHeatFlux()};
    Qn[] = ValueFromIndex[]{nodalNormalHeatFlux()};
  EndIf
 
}

// Jacobian
Jacobian{
  { Name JVolume;
    Case  {
      { Region All; Jacobian Vol;}
    }
  }
  { Name JSurface;
    Case  {
      { Region All; Jacobian Sur;}
    }
  }
  { Name JLine;
    Case  {
      { Region All; Jacobian Lin;}
    }
  }
}

//Integration
Integration {
  { Name Int1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Line ;          NumberOfPoints  4 ; }
          { GeoElement Triangle ;      NumberOfPoints  6 ; }
          { GeoElement Quadrangle ;    NumberOfPoints  4 ; }
          { GeoElement Tetrahedron ;   NumberOfPoints  4 ; }
          { GeoElement Hexahedron ;    NumberOfPoints  6 ; }
          { GeoElement Prism ;         NumberOfPoints  9 ; }
        }
      }
    }
  }
}

//Constraint
Constraint{
  { Name DirichletBC_IN;
    Case {
      {Region Inner; Type Assign; Value T_in[];}
    }
  }
  { Name NeumannBC;
    Case {
      {Region Outer; Value -Qn[]/lambda;}
    }
  }
}

//FunctionSpace
FunctionSpace{
  { Name H_T; Type Form0;
    BasisFunction{
      { Name sTn; NameOfCoef Tn; Function BF_Node;
        Support AllOmega; Entity NodesOf[All];
      }
    }
    Constraint{
      { NameOfCoef Tn; EntityType NodesOf; NameOfConstraint DirichletBC_IN; }
    }
  }
  { Name H_GradT; Type Form0;
    BasisFunction{
      { Name sGradTn; NameOfCoef GradTn; Function BF_Node;
        Support Outer; Entity NodesOf[All];
      }
    }
    Constraint{
      { NameOfCoef GradTn; EntityType NodesOf; NameOfConstraint NeumannBC; }
    }
  }
}

//Formulation (weak form)
Formulation{
  { Name HeatEquation; Type FemEquation;
    Quantity{
      {Name T; Type Local; NameOfSpace H_T;}
      {Name GradT; Type Local; NameOfSpace H_GradT;}
    }
    Equation{
      Galerkin{[-Dof{Grad T}, {Grad T}]; In Omega; Jacobian JVolume; Integration Int1;}
      Galerkin{[Dof{GradT}, {T}]; In Outer; Jacobian JSurface; Integration Int1;}
    }
  }
}

//Resolution
Resolution{
  { Name HeatEquation;
    System{
      {Name Syst; NameOfFormulation HeatEquation;}
    }
    Operation {
      If(!Initialize)
        Generate[Syst]; Solve[Syst]; //SaveSolution[Syst];
      EndIf
      PostOperation[nodalFields];
      PostOperation[plot];
    }
  }
}

//Post processing
PostProcessing{
  { Name HeatEquation; NameOfFormulation HeatEquation;
    Quantity{
      {Name Pos; Value {Local{[ Vector[X[], Y[], 0]]; In AllOmega; Jacobian JVolume;}} }
      {Name Disp; Value {Local{[ Vector[0, 0, 0]]; In Omega; Jacobian JVolume;}} }
      {Name Vel; Value {Local{[ Vector[0, 0, 0]]; In AllOmega; Jacobian JVolume;}} }
      {Name Force; Value {Local{[ Vector[0, 0, 0]]; In Outer; Jacobian JVolume;}} } 
      {Name Temp; Value {Local{[{T}]; In AllOmega; Jacobian JVolume;}} }
      {Name HF; Value {Local{[-lambda*{Grad T}]; In AllOmega; Jacobian JVolume;}} } 
    }
  }  
}

//Post operations
PostOperation{
  { Name nodalFields; NameOfPostProcessing HeatEquation; Format NodeTable;
    Operation{
      Print [Pos, OnElementsOf Outer, LastTimeStepOnly, Name "nodalPosition",
        File StrChoice[OutputFiles, "nodalPosition.txt", ""]];
      Print [Disp, OnElementsOf Omega, LastTimeStepOnly, Name "nodalDisplacement",
        File StrChoice[OutputFiles, "nodalDisplacement.txt", ""]];
      Print [Disp, OnElementsOf Omega, LastTimeStepOnly, Name "nodalDisplacementNm1",
        File StrChoice[OutputFiles, "nodalDisplacementNm1.txt", ""]];
      Print [Disp, OnElementsOf Omega, LastTimeStepOnly, Name "nodalDisplacementNm2",
        File StrChoice[OutputFiles, "nodalDisplacementNm2.txt", ""]];
      Print [Vel, OnElementsOf Outer, LastTimeStepOnly, Name "nodalVelocity",
        File StrChoice[OutputFiles, "nodalVelocity.txt", ""]];
      Print [Force, OnElementsOf Outer, LastTimeStepOnly, Name "nodalForce",
        File StrChoice[OutputFiles, "nodalForce.txt", ""]];
      Print [Temp, OnGrid Outer, LastTimeStepOnly, Name "nodalTemperature",
        File StrChoice[OutputFiles, "nodalTemperature.txt", ""]];
        Print [HF, OnGrid Outer, LastTimeStepOnly, Name "nodalHeatFlux",
          File StrChoice[OutputFiles, "nodalHeatFlux.txt", ""]];
    }
  }
  { Name plot; NameOfPostProcessing HeatEquation;
    Operation {
      Print [Temp, OnElementsOf Omega, File "T.pos"];
      Print[HF, OnElementsOf Omega, File "Q.pos"];
    }
  }
}
