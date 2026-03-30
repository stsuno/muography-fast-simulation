#ifndef MATERIALLIST_H
#define MATERIALLIST_H

#include <string>
#include <TColor.h>
#include <TError.h>
#include <TSystem.h>

static double GetDensity(const std::string& material_name) {
  if (material_name == "Infinity")          return 1000.0;          // [g/cm3] PDG
                                                                    // PDG 2025
  if (material_name == "StandardRock")      return 2.650;           // [g/cm3] PDG
  if (material_name == "ShieldingConcrete") return 2.300;           // [g/cm3] PDG
  if (material_name == "Polyvinyltoluene")  return 1.030;           // [g/cm3] PDG  (Scintillator)
  if (material_name == "Polyethylene")      return 0.890;           // [g/cm3] PDG  PE
  if (material_name == "Iron")              return 7.874;           // [g/cm3] PDG  Fe
  if (material_name == "Water")             return 1.000;           // [g/cm3] PDG  H2O
  if (material_name == "Air")               return 1.205 / 1000.0;  // [g/cm3] PDG  (dry, 1atm)
                                                                    // Engineering / Custom Values
  if (material_name == "NormalSoil")        return 2.000;           // [g/cm3] Civil Engineering (JIS A 1210)
  if (material_name == "Polyvinyl")         return 1.400;           // [g/cm3] Typical Rigid PVC
  if (material_name == "CastIron")          return 7.200;           // [g/cm3] Typical Ductile Iron
  if (material_name == "Steel")             return 7.850;           // [g/cm3] Industrial standard (JIS G 3192)
  if (material_name == "MuddyWater")        return 1.200;           // [g/cm3] Engineering value
  Error("GetDensity", "Unknown material name: %s", material_name.c_str());
  gSystem->Exit(EXIT_FAILURE);
}

static Color_t GetColor(const std::string& color_name) {
  return (Color_t)gInterpreter->Calc(color_name.c_str());
}

#endif

