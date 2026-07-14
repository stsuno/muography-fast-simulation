#ifndef PIPE_H
#define PIPE_H

#include "Simulation/CompositeShape.h"
#include "Simulation/Cylinder.h"

// Pipe class consisting of two concentric Cylinders (outer shell and inner hollow).
class Pipe : public CompositeShape {
  public:
    Pipe(double cx, double cy, double cz, double length, double inner_radius, double outer_radius);
    virtual ~Pipe() = default;

    // Accessors (Setters and Getters)
    void SetShellDensity(double density);
    double GetShellDensity() const { return m_shell_density; }
    void SetHollowDensity(double density);
    double GetHollowDensity() const { return m_hollow_density; }

    void SetShellAtomicNumber(double atomic_number);
    double GetShellAtomicNumber() const { return m_shell_atomic_number; }
    void SetHollowAtomicNumber(double atomic_number);
    double GetHollowAtomicNumber() const { return m_hollow_atomic_number; }

    void SetHollowAtomicMass(double atomic_mass);
    double GetShellAtomicMass() const { return m_shell_atomic_mass; }
    void SetShellAtomicMass(double atomic_mass);
    double GetHollowAtomicMass() const { return m_hollow_atomic_mass; }

    void SetShellMeanExcitationEnergy(double mean_excitation_energy);
    double GetShellMeanExcitationEnergy() const { return m_shell_mean_excitation_energy; }
    void SetHollowMeanExcitationEnergy(double mean_excitation_energy);
    double GetHollowMeanExcitationEnergy() const { return m_hollow_mean_excitation_energy; }

  private:
    // Member variables
    double m_shell_density = 0.0;
    double m_hollow_density = 0.0;

    double m_hollow_atomic_number = 0.0;
    double m_hollow_atomic_mass = 0.0;
    double m_hollow_mean_excitation_energy = 0.0;

    double m_shell_atomic_number = 0.0;
    double m_shell_atomic_mass = 0.0;
    double m_shell_mean_excitation_energy = 0.0;

    // Geometric components
    Cylinder* m_outer_cylinder;
    Cylinder* m_inner_cylinder;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

// Constructor
inline Pipe::Pipe(double cx, double cy, double cz, double length, double inner_radius, double outer_radius) {
  m_inner_cylinder = new Cylinder(cx, cy, cz, length, inner_radius);
  m_outer_cylinder = new Cylinder(cx, cy, cz, length, outer_radius);

  AddShape(m_inner_cylinder);
  AddShape(m_outer_cylinder);
}

inline void Pipe::SetShellDensity(double density) {
  m_shell_density = density;
  if (m_outer_cylinder) m_outer_cylinder->SetDensity(density);
}

inline void Pipe::SetHollowDensity(double density) {
  m_hollow_density = density;
  if (m_inner_cylinder) m_inner_cylinder->SetDensity(density);
}

inline void Pipe::SetShellAtomicNumber(double atomic_number) {
  m_shell_atomic_number = atomic_number;
  if (m_outer_cylinder) m_outer_cylinder->SetAtomicNumber(atomic_number);
}

inline void Pipe::SetHollowAtomicNumber(double atomic_number) {
  m_hollow_atomic_number = atomic_number;
  if (m_inner_cylinder) m_inner_cylinder->SetAtomicNumber(atomic_number);
}

inline void Pipe::SetShellAtomicMass(double atomic_mass) {
  m_shell_atomic_mass = atomic_mass;
  if (m_outer_cylinder) m_outer_cylinder->SetAtomicMass(atomic_mass);
}

inline void Pipe::SetHollowAtomicMass(double atomic_mass) {
  m_hollow_atomic_mass = atomic_mass;
  if (m_inner_cylinder) m_inner_cylinder->SetAtomicMass(atomic_mass);
}

inline void Pipe::SetShellMeanExcitationEnergy(double mean_excitation_energy) {
  m_shell_mean_excitation_energy = mean_excitation_energy;
  if (m_outer_cylinder) m_outer_cylinder->SetMeanExcitationEnergy(mean_excitation_energy);
}

inline void Pipe::SetHollowMeanExcitationEnergy(double mean_excitation_energy) {
  m_hollow_mean_excitation_energy = mean_excitation_energy;
  if (m_inner_cylinder) m_inner_cylinder->SetMeanExcitationEnergy(mean_excitation_energy);
}

#endif
