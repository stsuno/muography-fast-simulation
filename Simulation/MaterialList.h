#ifndef MATERIALLIST_H
#define MATERIALLIST_H

#include <map>
#include <string>

#include <TColor.h>
#include <TError.h>
#include <TSystem.h>

#include "Simulation/Shape.h"

// Static-only registry of material properties, keyed by material name.
// Adding a new material = adding one row to Table() below.
//
// Sources:
//   - PDG Atomic and Nuclear Properties of Materials (2025 edition)
//     https://pdg.lbl.gov/2025/AtomicNuclearProperties/
//     <Z/A> is quoted directly from the PDG tables for compounds/mixtures;
//     for elements it is Z/A (Fe: 26/55.845, H: 1/1.008).
//   - Sternheimer coefficients are quoted verbatim from the header line
//     "Sternheimer coef:  a     k=m_s   x_0    x_1    I[eV]   Cbar  delta0"
//     of the PDG muon energy-loss tables (MUE/*.txt) for each material.
//   - Engineering values (soil, PVC, cast iron, steel, muddy water) follow
//     civil-engineering standards (JIS) or typical handbook values.
//   - Radiation lengths X0 [g/cm2] are quoted from the same PDG pages
//     ("Radiation length" row). They are used only by the multiple-Coulomb-
//     scattering model in Simulation/MuonTransport.h, never by Bethe-Bloch.
//     NormalSoil reuses the standard-rock X0 and CastIron the iron X0: both are
//     engineering approximations, flagged in the table below.
//
// Note: "Steel" and "Infinity" reuse the PDG iron (Fe) coefficients.
// (PDG iron density is 7.874 g/cm3 vs. 7.850 g/cm3 used here for steel; the
// resulting error in the density-effect correction is negligible.)
class MaterialList {
  public:
    MaterialList() = delete;  // static-only class

    // Materials for which only the density is known keep kUnset in the other
    // fields of Shape::Material; accessing an unset field through the
    // per-parameter getters below is an error (Error + Exit).
    static const Shape::Material& GetMaterial(const std::string& material_name);

    // Per-parameter getters (kept for callers that need a single property).
    static double GetDensity(const std::string& material_name);
    static double GetZOverA(const std::string& material_name);
    static double GetMeanExcitationEnergy(const std::string& material_name);
    static Shape::SternheimerParameters GetSternheimerParameters(const std::string& material_name);
    // Radiation length X0 [g/cm2] (Error + Exit if not registered).
    static double GetRadiationLength(const std::string& material_name);
    // Radiation length expressed as a length [cm], i.e. X0 / density.
    static double GetRadiationLengthCM(const std::string& material_name);

  private:
    // Returns the value unchanged, or fails loudly (Error + Exit) if it is
    // still kUnset for the requested material.
    static double RequireDefined(double value, const char* where, const std::string& material_name);

    // Function-local static: safe against initialization-order issues when this
    // header is included from multiple translation units.
    static const std::map<std::string, Shape::Material>& Table();
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline const Shape::Material& MaterialList::GetMaterial(const std::string& material_name) {
  const auto& table = Table();
  auto it = table.find(material_name);
  if (it == table.end()) {
    Error("MaterialList::GetMaterial", "Unknown material name: %s", material_name.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }
  return it->second;
}

inline double MaterialList::GetDensity(const std::string& material_name) {
  return RequireDefined(GetMaterial(material_name).m_density, "GetDensity", material_name);
}

inline double MaterialList::GetZOverA(const std::string& material_name) {
  return RequireDefined(GetMaterial(material_name).m_z_over_a, "GetZOverA", material_name);
}

inline double MaterialList::GetMeanExcitationEnergy(const std::string& material_name) {
  return RequireDefined(GetMaterial(material_name).m_mean_excitation_energy,
                        "GetMeanExcitationEnergy", material_name);
}

inline Shape::SternheimerParameters MaterialList::GetSternheimerParameters(const std::string& material_name) {
  const Shape::Material& material = GetMaterial(material_name);
  // Default-constructed (all-zero) parameters mean "not defined": Cbar > 0
  // for every real material (see Shape::Material).
  if (material.m_sternheimer.m_cbar <= 0.0) {
    Error("MaterialList::GetSternheimerParameters",
          "Sternheimer parameters not defined for: %s", material_name.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }
  return material.m_sternheimer;
}

inline double MaterialList::GetRadiationLength(const std::string& material_name) {
  return RequireDefined(GetMaterial(material_name).m_radiation_length,
                        "GetRadiationLength", material_name);
}

inline double MaterialList::GetRadiationLengthCM(const std::string& material_name) {
  return GetRadiationLength(material_name) / GetDensity(material_name);
}

inline double MaterialList::RequireDefined(double value, const char* where, const std::string& material_name) {
  if (value == Shape::Material::kUnset) {
    Error(Form("MaterialList::%s", where),
          "Property not defined for: %s", material_name.c_str());
    gSystem->Exit(EXIT_FAILURE);
  }
  return value;
}

inline const std::map<std::string, Shape::Material>& MaterialList::Table() {
  constexpr double kUnset = Shape::Material::kUnset;
  static const std::map<std::string, Shape::Material> table = {
    //                    density    <Z/A>    I [MeV]    Sternheimer {a, k, x0, x1, Cbar, delta0}                X0 [g/cm2]
    // PDG 2025
    {"Infinity",          {1000.0,   0.46557, 2.86e-4,  {0.1468, 2.9632, -0.0012, 3.1531,  4.2911, 0.12},        13.84}},  // Fe props (MUE/muE_iron_Fe.txt)
    {"ShieldingConcrete", {2.300,    0.50274, 1.352e-4, {0.0751, 3.5467,  0.1301, 3.0466,  3.9464, 0.00},        26.57}},  // MUE/muE_shielding_concrete.txt
    {"Polyvinyltoluene",  {1.030,    0.54141, 6.47e-5,  {0.1610, 3.2393,  0.1464, 2.4855,  3.1997, 0.00},        43.90}},  // Scintillator, MUE/muE_polyvinyltoluene.txt
    {"Air",               {1.205e-3, 0.49919, 8.57e-5,  {0.1091, 3.3994,  1.7418, 4.2759, 10.5961, 0.00},        36.62}},  // dry, 1 atm, MUE/muE_air_dry_1_atm.txt
    {"HydrogenLiquid",    {0.0708,   0.99212, 2.18e-5,  {0.1348, 5.6249,  0.4400, 1.8856,  3.0977, 0.00},        63.04}},  // MUE/muE_hydrogen_liquid.txt
    {"StandardRock",      {2.650,    0.50000, 1.364e-4, {0.0830, 3.4120,  0.0492, 3.0549,  3.7738, 0.00},        26.54}},  // MUE/muE_standard_rock.txt
    {"Polyethylene",      {0.890,    kUnset,  kUnset,   {},                                                      44.77}},  // PDG PE (density + X0)
    {"Iron",              {7.874,    kUnset,  kUnset,   {},                                                      13.84}},  // PDG Fe (density + X0)
    {"Water",             {1.0000,   0.55509, 7.97e-5,  {0.09116, 3.4773,  0.2400, 2.8004, 3.5017, 0.00},        36.08}}, // PDG H2O (density + X0)
    // Engineering / Custom Values
    {"Steel",             {7.850,    0.46557, 2.86e-4,  {0.1468, 2.9632, -0.0012, 3.1531,  4.2911, 0.12},        13.84}},  // JIS G 3192 density, Fe props
    {"NormalSoil",        {2.000,    kUnset,  kUnset,   {},                                                      26.54}},  // Civil Engineering (JIS A 1210), X0 approximated by standard rock
    {"Polyvinyl",         {1.400}},                                                                                          // Typical Rigid PVC
    {"CastIron",          {7.200,    kUnset,  kUnset,   {},                                                      13.84}},  // Typical Ductile Iron, X0 approximated by Fe
    {"MuddyWater",        {1.200}},                                                                                          // Engineering value
  };
  return table;
}

static Color_t GetColor(const std::string& color_name) {
  return (Color_t)gInterpreter->Calc(color_name.c_str());
}

#endif
