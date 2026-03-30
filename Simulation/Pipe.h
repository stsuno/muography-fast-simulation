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

  private:
    // Member variables
    double m_shell_density = 0.0;
    double m_hollow_density = 0.0;

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

#endif
