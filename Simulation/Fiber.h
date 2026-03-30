#ifndef FIBER_H
#define FIBER_H

#include <iostream>

#include "Simulation/CompositeShape.h"
#include "Simulation/RectangularPrism.h"

class Fiber : public CompositeShape {
  public:
    Fiber(double cx, double cy, double cz, double radius, double windingangle, int direction);
    virtual ~Fiber() = default;

    static constexpr double s_fiber_width    =   0.2;
    static constexpr double s_element_length =   5.0;
    static constexpr double s_fiber_length   = 230.0;
    static constexpr int    s_fiber_element  = s_fiber_length/s_element_length;
    static constexpr double s_timing_resolution = 1.0;     // [ns]

    // シンチレーションファイバー内での光の速度
    // v = c（光速）/n（屈折）
    // v ≈ (3.00×10⁸)/1.59 ≈ 1.89×10⁸m/s = 1.89×10⁻¹m/ns = 18.9 cm/ns
    static constexpr double s_signal_speed = 18.9;   // [cm/ns]

    // Geometry Operations
    void MoveXYZ(double dx, double dy, double dz) override;
    void RotatePhi(double phi) override;
    void RotateTheta(double theta) override;

    double GetFiberLength()   { return s_fiber_length; }
    double GetFiberWidth()    { return s_fiber_width; }
    double GetElementLength() { return s_element_length; }
    int    GetNumElement()    { return s_fiber_element; }
    double GetSignalSpeed()   { return s_signal_speed; }
    double GetTimingResolution() { return s_timing_resolution; }

    TVector3 GetHitIDCenter(int id) { return m_scintillator[id]->GetCenter(); }
    int GetHitID(const Muon& muon);
    TVector3 GetHitPosition(const Muon& muon);

    std::vector<int> IsHitList(const Muon& muon);

  private:
    std::vector<RectangularPrism*> m_scintillator;
    TVector3 m_detectorCenter;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline Fiber::Fiber(double cx, double cy, double cz, double radius, double windingangle, int direction) {

  if (windingangle<=0.0) {
    Info("Fiber::Fiber", "Winding angle must be positive [degree].");
    abort();
  }
  double pi = TMath::Pi();
  m_scintillator.resize(s_fiber_element,nullptr);
  m_detectorCenter = TVector3(cx,cy,cz);

  double wangle = windingangle*pi/180.0;    // radian
  double sin    = TMath::Sin(wangle);
  double cos    = TMath::Cos(wangle);
  double height = s_fiber_length*sin;

  for (int i=0; i<s_fiber_element; i++) {
//    double degree = direction*s_fiber_width*i*cos/radius;
    double degree = direction*s_element_length*i*cos/radius;
    double gx = radius*TMath::Cos(degree);
    double gy = radius*TMath::Sin(degree);
    double gz = height*(-0.5+1.0*i/s_fiber_element);
//    m_scintillator[i] = new RectangularPrism(gx, gy, gz, s_fiber_width, s_fiber_width, s_fiber_width);
    m_scintillator[i] = new RectangularPrism(gx, gy, gz, s_fiber_width, s_element_length, s_fiber_width);
    m_scintillator[i] -> RotateTheta(wangle);
    m_scintillator[i] -> RotatePhi(degree);
  }

  m_scintillator[0] -> SetLineColor(1);
  for (int i=0; i<s_fiber_element; i++) {
    m_scintillator[i] -> MoveXYZ(cx,cy,cz);
    m_scintillator[i] -> SetUniqueID(i);
    AddShape(m_scintillator[i]);
  }
}

inline void Fiber::MoveXYZ(double dx, double dy, double dz) {
  TVector3 displacement(dx, dy, dz);
  m_detectorCenter += displacement;
  for (int i=0; i<s_fiber_element; i++) {
    std::vector<TVector3> vtx = m_scintillator[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex += displacement;
      m_scintillator[i]->UpdateVertex(j, vertex);
    }
    m_scintillator[i]->ComputeCenter();
    m_scintillator[i]->ComputeSurfaceEquation();
  }
}

inline void Fiber::RotatePhi(double phi) {
  for (int i=0; i<s_fiber_element; i++) {
    std::vector<TVector3> vtx = m_scintillator[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex -= m_detectorCenter;
      vertex.RotateZ(phi);
      vertex += m_detectorCenter;
      m_scintillator[i]->UpdateVertex(j, vertex);
    }
    m_scintillator[i]->ComputeCenter();
    m_scintillator[i]->ComputeSurfaceEquation();
  }
}

inline void Fiber::RotateTheta(double theta) {
  for (int i=0; i<s_fiber_element; i++) {
    std::vector<TVector3> vtx = m_scintillator[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex -= m_detectorCenter;
      vertex.RotateX(theta);
      vertex += m_detectorCenter;
      m_scintillator[i]->UpdateVertex(j, vertex);
    }
    m_scintillator[i]->ComputeCenter();
    m_scintillator[i]->ComputeSurfaceEquation();
  }
}

inline std::vector<int> Fiber::IsHitList(const Muon& muon) {
  std::vector<int> hitlist;
  for (int i=0; i<s_fiber_element; i++) {
    if (m_scintillator[i]->IsHit(muon)) {
      hitlist.push_back(i);
    }
  }
  return hitlist;
}

inline int Fiber::GetHitID(const Muon& muon) {
  std::vector<int> hitlist = IsHitList(muon);
  int minID = 1000000;
  for (size_t i=0; i<hitlist.size(); i++) {
    minID = TMath::Min(hitlist.at(i),minID); 
  }
  return (minID>10000) ? 0 : minID;
}

inline TVector3 Fiber::GetHitPosition(const Muon& muon) {
  int id = GetHitID(muon);
  std::vector<TVector3> hitpoints = m_scintillator[id]->HitPoints(muon);
  TVector3 point(0.0,0.0,0.0);
  int nface = 0;
  for (size_t i=0; i<hitpoints.size(); i++) {
    point += hitpoints.at(i);
    nface++;
  }
  if (nface>0) { point *= 1.0/nface; }
  return point;
}

#endif
