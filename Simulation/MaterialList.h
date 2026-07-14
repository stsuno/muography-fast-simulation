#ifndef MATERIALLIST_H
#define MATERIALLIST_H

#include <string>
#include <TColor.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/Shape.h"

static double GetDensity(const std::string& material_name) {
  // https://pdg.lbl.gov/2026/AtomicNuclearProperties/HTML/
  if (material_name == "Infinity")          return 1000.0;          // [g/cm3] PDG
                                                                    // PDG 2025
  if (material_name == "StandardRock")      return 2.650;           // [g/cm3] PDG
  if (material_name == "ShieldingConcrete") return 2.300;           // [g/cm3] PDG
  if (material_name == "Polyvinyltoluene")  return 1.030;           // [g/cm3] PDG  (Scintillator)
  if (material_name == "Polyethylene")      return 0.890;           // [g/cm3] PDG  PE
  if (material_name == "Iron")              return 7.874;           // [g/cm3] PDG  Fe
  if (material_name == "Water")             return 1.000;           // [g/cm3] PDG  H2O
  if (material_name == "Air")               return 1.205 / 1000.0;  // [g/cm3] PDG  (dry, 1atm)
  if (material_name == "HydrogenLiquid")    return 0.0708;          // [g/cm3] PDG 
                                                                    // Engineering / Custom Values
  if (material_name == "NormalSoil")        return 2.000;           // [g/cm3] Civil Engineering (JIS A 1210)
  if (material_name == "Polyvinyl")         return 1.400;           // [g/cm3] Typical Rigid PVC
  if (material_name == "CastIron")          return 7.200;           // [g/cm3] Typical Ductile Iron
  if (material_name == "Steel")             return 7.850;           // [g/cm3] Industrial standard (JIS G 3192)
  if (material_name == "MuddyWater")        return 1.200;           // [g/cm3] Engineering value
  Error("GetDensity", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

static double GetAtomicNumber(const std::string& material_name) {
  if (material_name == "Infinity")          return 26.0;            // PDG        

  if (material_name == "StandardRock")      return 11.0;            // PDG
  if (material_name == "ShieldingConcrete") return 8.55;            // PDG
  if (material_name == "Polyvinyltoluene")  return 3.37;            // PDG  (Scintillator)                    
  if (material_name == "Iron")              return 26.0;           //  PDG  Fe
  if (material_name == "Air")               return 7.26;            // PDG  (dry, 1atm)
  if (material_name == "HydrogenLiquid")    return 1.0;             // PDG
  if (material_name == "Water")             return 0.5550;           // PDG  H2O

  if (material_name == "Steel")             return 26.0;            // same Infinity
  Error("GetAtomicNumber", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

static double GetAtomicMass(const std::string& material_name) {
  if (material_name == "Infinity")          return 55.845;          // [g/mol] PDG           

  if (material_name == "StandardRock")      return 22.0;            // [g/mol]PDG
  if (material_name == "ShieldingConcrete") return 17.01;           // [g/mol]PDG
  if (material_name == "Polyvinyltoluene")  return 6.23;            // [g/mol] PDG  (Scintillator)                   
  if (material_name == "Iron")              return 55.845;           // [g/cm3] PDG  Fe
  if (material_name == "Air")               return 14.55;           // [g/mol] PDG  (dry, 1atm)
  if (material_name == "HydrogenLiquid")    return 1.008;           // [g/mol] PDG  
  if (material_name == "Water")             return 1.000;           // [g/mol] PDG  H2O

  if (material_name == "Steel")             return 55.845;          // [g/mol] same Infinity
  Error("GetAtomicMass", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

static double GetMeanExcitationEnergy(const std::string& material_name) {
  if (material_name == "Infinity")          return 2.86e-4;         // [Mev] PDG          

  if (material_name == "StandardRock")      return 1.364e-4;        // [Mev] PDG
  if (material_name == "ShieldingConcrete") return 1.352e-4;        // [Mev] PDG
  if (material_name == "Polyvinyltoluene")  return 6.47e-5;         // [Mev] PDG  (Scintillator)                   
  if (material_name == "Iron")              return 2.86e-4;         // [MeV] PDG  Fe
  if (material_name == "Air")               return 8.57e-5;         // [Mev] PDG  (dry, 1atm)
  if (material_name == "HydrogenLiquid")    return 2.18e-5;         // [Mev] PDG   
  if (material_name == "Water")             return 7.97e-5;         // [Mev] PDG   

  if (material_name == "Steel")             return 2.86e-4;         // [Mev] same Infinity
  Error("GetMeanExcitationEnergy", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

// Sternheimer coefficients for the density-effect correction delta(beta*gamma)
// used in the Bethe-Bloch formula (PDG Eq. 34.7).
//
// Source: PDG Atomic and Nuclear Properties of Materials (2026 edition)
//         https://pdg.lbl.gov/2026/AtomicNuclearProperties/
// The coefficients below are quoted verbatim from the header line
// "Sternheimer coef:  a     k=m_s   x_0    x_1    I[eV]   Cbar  delta0"
// of the muon energy-loss table (MUE/*.txt) for each material:
//   StandardRock:      MUE/muE_standard_rock.txt
//                      a=0.0830 k=3.4120 x0=0.0492  x1=3.0549 Cbar=3.7738  delta0=0.00 (I=136.3 eV)
//   ShieldingConcrete: MUE/muE_shielding_concrete.txt
//                      a=0.0751 k=3.5467 x0= 0.1301 x1=3.0466 Cbar= 3.9464 delta0=0.00 (I=135.2 eV)
//   Polyvinyltoluene : MUE/muE_polyvinyltoluene.txt
//                      a=0.1610 k=3.2393 x0= 0.1464 x1=2.4855 Cbar= 3.1997 delta0=0.00 (I= 64.7 eV)
//   Air (dry, 1 atm) : MUE/muE_air_dry_1_atm.txt
//                      a=0.1091 k=3.3994 x0= 1.7418 x1=4.2759 Cbar=10.5961 delta0=0.00 (I= 85.7 eV)
//   Iron (Fe)        : MUE/muE_iron_Fe.txt
//                      a=0.1468 k=2.9632 x0=-0.0012 x1=3.1531 Cbar= 4.2911 delta0=0.12 (I=286.0 eV)
//   HydrogenLiquid   : MUE/muE_hydrogen_liquid.txt
//                      a=0.1348 k=5.6249 x0= 0.4400 x1=1.8856 Cbar= 3.0977 delta0=0.00 (I= 21.8 eV)
//   Water            : MUE/muE_water_liquid.txt
//                      a=0.0912 k=3.4773 x0= 0.2400 x1=2.8004 Cbar= 3.5017 delta0=0.00 (I= 79.7 eV)
//
// Note: "Steel" and "Infinity" reuse the PDG iron (Fe) coefficients, consistent
// with GetAtomicNumber()/GetAtomicMass()/GetMeanExcitationEnergy() above.
// (PDG iron density is 7.874 g/cm3 vs. 7.850 g/cm3 used here for steel; the
// resulting error in the density-effect correction is negligible.)
static SternheimerParameters GetSternheimerParameters(const std::string& material_name) {
  //                                                          a       k       x0      x1      Cbar    delta0
  if (material_name == "Infinity")          return {0.1468, 2.9632, -0.0012, 3.1531,  4.2911, 0.12};  // PDG (Fe)

  if (material_name == "StandardRock")      return {0.0830, 3.4120,  0.0492, 3.0549,  3.7738, 0.00};  // PDG
  if (material_name == "ShieldingConcrete") return {0.0751, 3.5467,  0.1301, 3.0466,  3.9464, 0.00};  // PDG
  if (material_name == "Polyvinyltoluene")  return {0.1610, 3.2393,  0.1464, 2.4855,  3.1997, 0.00};  // PDG  (Scintillator)
  if (material_name == "Iron")              return {0.1468, 2.9632, -0.0012, 3.1531,  4.2911, 0.12};
  if (material_name == "Air")               return {0.1091, 3.3994,  1.7418, 4.2759, 10.5961, 0.00};  // PDG  (dry, 1atm)
  if (material_name == "HydrogenLiquid")    return {0.1348, 5.6249,  0.4400, 1.8856,  3.0977, 0.00};  // PDG
  if (material_name == "Water")             return {0.0912, 3.4773, 0.2400, 2.8004, 3.5017, 0.00}; // PDG

  if (material_name == "Steel")             return {0.1468, 2.9632, -0.0012, 3.1531,  4.2911, 0.12};  // same Infinity (Fe)
  Error("GetSternheimerParameters", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

static Color_t GetColor(const std::string& color_name) {
  return (Color_t)gInterpreter->Calc(color_name.c_str());
}

#endif

