#ifndef SCIFIBERCYLINDER_H
#define SCIFIBERCYLINDER_H

#include "Simulation/CompositeShape.h"
#include "Simulation/Fiber.h"

class SciFiberCylinder : public CompositeShape {
  public:
    SciFiberCylinder(double cx, double cy, double cz);
    virtual ~SciFiberCylinder() = default;

    static constexpr int    s_nlayers            = 3;
    static constexpr int    s_nchannel_per_layer = 16;

//    static constexpr int    s_nlayers            = 1;
//    static constexpr int    s_nchannel_per_layer = 16;

    static constexpr double s_radius             = 10.0;
    static constexpr double s_space[3]           = {0.0, 0.2, 0.2};
    static constexpr double s_winding_angle[3]   = {2.91931, 2.91931, 3.4};
    static constexpr double s_ribbon_rotation[3] = {0.0, 60.0, 120.0};
    static constexpr int    s_ribbon_color[3]    = {3, 6, 5};

    // Geometry Operations
    void MoveXYZ(double dx, double dy, double dz) override;
    void RotatePhi(double phi) override;
    void RotateTheta(double theta) override;

    double GetDetectorHeight();
    Fiber* GetFiber(int id) { return m_fiber.at(id); }

    std::vector<std::pair<int, std::vector<int>>> IsHitList(const Muon& muon);

  private:
    std::vector<Fiber*> m_fiber;
    TVector3 m_detectorCenter;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline SciFiberCylinder::SciFiberCylinder(double cx, double cy, double cz) {

  double pi = TMath::Pi();
  m_fiber.resize(s_nlayers*s_nchannel_per_layer,nullptr);
  m_detectorCenter = TVector3(cx,cy,cz);

  double s_fiber_radius[3] = {s_radius+s_space[0], s_radius+s_space[1], s_radius+2.0*s_space[2]};
  int direction[3]         = {1,-1,1};

  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_fiber[s_nchannel_per_layer*i+j] = new Fiber(cx, cy, cz, s_fiber_radius[i], s_winding_angle[i], direction[i]);

      double width = m_fiber[s_nchannel_per_layer*i+j]->GetFiberWidth();
      m_fiber[s_nchannel_per_layer*i+j] -> MoveXYZ(0.0,0.0,width*j);
      m_fiber[s_nchannel_per_layer*i+j] -> RotatePhi(s_ribbon_rotation[i]*pi/180.0);
      m_fiber[s_nchannel_per_layer*i+j] -> SetLineColor(s_ribbon_color[i]);
      AddShape(m_fiber[s_nchannel_per_layer*i+j]);
    }
  }
}

inline void SciFiberCylinder::MoveXYZ(double dx, double dy, double dz) {
  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_fiber[s_nchannel_per_layer*i+j] -> MoveXYZ(dx,dy,dz);
    }
  }
}

inline void SciFiberCylinder::RotatePhi(double phi) {
  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_fiber[s_nchannel_per_layer*i+j] -> RotatePhi(phi);
    }
  }
}

inline void SciFiberCylinder::RotateTheta(double theta) {
  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_fiber[s_nchannel_per_layer*i+j] -> RotateTheta(theta);
    }
  }
}

inline double SciFiberCylinder::GetDetectorHeight() {
  double wid = m_fiber[0]->GetFiberWidth();
  double len = m_fiber[0]->GetFiberLength();
  double sin = TMath::Max(TMath::Sin(s_winding_angle[0]), TMath::Sin(s_winding_angle[1]));
  sin = TMath::Max(TMath::Sin(s_winding_angle[2]), sin);
  return len*sin+(s_nchannel_per_layer-1)*wid;
}

inline std::vector<std::pair<int, std::vector<int>>> SciFiberCylinder::IsHitList(const Muon& muon) {
  std::vector<std::pair<int, std::vector<int>>> hitlist;
  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      std::vector<int> hitlocal = m_fiber[s_nchannel_per_layer*i+j]->IsHitList(muon);
      if (hitlocal.size()>0) {
        hitlist.push_back(std::make_pair(s_nchannel_per_layer*i+j,hitlocal));
      }
    }
  }
  return hitlist;
}

#endif
