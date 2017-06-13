.. highlight:: none

Preparing for an MD calculation
===============================

The initial structure for starting a molecular dynamics simulation should
usually be structurally relaxed.

The geometry relaxation ::

  Geometry = GenFormat {
    <<< "geom.gen"
  }
  
  Driver = ConjugateGradient {
    MaxForceComponent = 1e-5
  }
  
  
  Hamiltonian = DFTB {
    SCC = Yes
    SCCTolerance = 1E-7
    Filling = Fermi {
      Temperature [Kelvin] = 400
    }
    SlaterKosterFiles = Type2FileNames {
      Prefix = "../../slako/"
      Separator = "-"
      Suffix = ".skf" 
    }
    MaxAngularMomentum = {
      C = "p"
      O = "p"
      H = "s"
    }
  }

Vibrational modes
=================
  
Vibrational modes ::

   Geometry = GenFormat {
     <<< "geom.gen"
   }
   
   Driver = SecondDerivatives {
     Delta = 1e-4
   }
   
   
   Hamiltonian = DFTB {
     SCC = Yes
     SCCTolerance = 1E-8
     Filling = Fermi {
       Temperature [Kelvin] = 400
     }
     SlaterKosterFiles = Type2FileNames {
       Prefix = "../../slako/"
       Separator = "-"
       Suffix = ".skf" 
     }
     MaxAngularMomentum = {
       C = "p"
       O = "p"
       H = "s"
     }
   }

Calculating the modes
~~~~~~~~~~~~~~~~~~~~~

The input for the modes code ::
   
   # Neededs the equilibrium geometry, at which the Hessian had been calculated
   Geometry = GenFormat { 
     <<< geom.gen
   }
   
   DisplayModes = {
    PlotModes = -20:-1          # Take the top 10 modes
    Animate = Yes               # make xyz files showing the atoms moving
   }
   
   # You need to specify the SK-files, as the mass of the elements is needed
   SlaterKosterFiles = Type2FileNames {
     Prefix = "../../slako/"
     Separator = "-"
     Suffix = ".skf"
   }
   
   # Include the Hessian, which was calculated by DFTB+
   Hessian = {
     <<< "hessian.out"
   }

