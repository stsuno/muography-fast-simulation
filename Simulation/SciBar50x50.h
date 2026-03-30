#ifndef SCIBAR50X50_H
#define SCIBAR50X50_H

#include "Simulation/CompositeShape.h"
#include "Simulation/RectangularPrism.h"

class SciBar50x50 : public CompositeShape {
  public:
    SciBar50x50(double cx, double cy, double cz);
    virtual ~SciBar50x50() = default;

    static constexpr int    s_nlayers            = 4;
    static constexpr int    s_nchannel_per_layer = 50;
    static constexpr int    s_total_channel      = s_nlayers*s_nchannel_per_layer*2;
    static constexpr double s_channel_width      = 1.0;
    static constexpr double s_offset             = -24.5;

    static constexpr double s_localzposition[4] = {-15.0, -5.0, 5.0, 15.0};

    int GetNumLayers() const          { return s_nlayers; }
    int GetNumIDs() const             { return s_total_channel; }
    int GetLayerIndex(int id) const   { return id/s_nchannel_per_layer; }
    int GetChannelIndex(int id) const { return id%s_nchannel_per_layer; }

    // Get local coordinate (X or Y) from strip index
    TVector3 GetLocalPosition(int id) const;
    // Get global coordinate from strip index
    TVector3 GetGlobalPosition(int id) const;
    TVector3 GetGlobalPositionError(int id) const;

    // Geometry Operations
    void MoveXYZ(double dx, double dy, double dz) override;
    void RotatePhi(double phi) override;
    void RotateTheta(double theta) override;

    std::vector<int> IsHitList(const Muon& muon);

  private:
    std::vector<RectangularPrism*> m_scintillator;
    TVector3 m_detectorCenter;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline SciBar50x50::SciBar50x50(double cx, double cy, double cz) {

  m_scintillator.resize(s_total_channel,nullptr);
  m_detectorCenter = TVector3(cx,cy,cz);

  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) { 
      m_scintillator[s_nchannel_per_layer*2*i+j] 
        = new RectangularPrism(0.0, 0.0,-0.5*s_channel_width, 
                               s_channel_width, s_nchannel_per_layer*s_channel_width, s_channel_width); 
    }
    for (int j=0; j<s_nchannel_per_layer; j++) { 
      m_scintillator[s_nchannel_per_layer*2*i+j+s_nchannel_per_layer] 
        = new RectangularPrism(0.0, 0.0, 0.5*s_channel_width, 
                               s_nchannel_per_layer*s_channel_width, s_channel_width, s_channel_width); 
    }
  }

  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_scintillator[s_nchannel_per_layer*2*i+j]                       -> MoveXYZ(s_offset+s_channel_width*j,0.0,s_localzposition[i]);
      m_scintillator[s_nchannel_per_layer*2*i+j+s_nchannel_per_layer]  -> MoveXYZ(0.0,s_offset+s_channel_width*j,s_localzposition[i]);
    }
  }

  for (int i=0; i<s_total_channel; i++) {
    m_scintillator[i] -> MoveXYZ(cx,cy,cz);
    m_scintillator[i] -> SetUniqueID(i);
    AddShape(m_scintillator[i]);
  }
}

inline void SciBar50x50::MoveXYZ(double dx, double dy, double dz) {
  TVector3 displacement(dx, dy, dz);
  m_detectorCenter += displacement;
  for (int i=0; i<s_total_channel; i++) {
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

inline void SciBar50x50::RotatePhi(double phi) {
  for (int i=0; i<s_total_channel; i++) {
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

inline void SciBar50x50::RotateTheta(double theta) {
  for (int i=0; i<s_total_channel; i++) {
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

inline std::vector<int> SciBar50x50::IsHitList(const Muon& muon) {
  std::vector<int> hitlist;
  for (int i=0; i<s_total_channel; i++) {
    if (m_scintillator[i]->IsHit(muon)) {
      hitlist.push_back(i);
    }
  }
  return hitlist;
}

inline TVector3 SciBar50x50::GetLocalPosition(int id) const {
  // Calculation logic based on the strip index
  int layer = GetLayerIndex(id)/2;
  return TVector3(s_offset+s_channel_width*GetChannelIndex(id),0.0,s_localzposition[layer]);
}

inline TVector3 SciBar50x50::GetGlobalPosition(int id) const {
  return m_scintillator[id]->GetCenter();
}

inline TVector3 SciBar50x50::GetGlobalPositionError(int id) const {
  int layer = GetLayerIndex(id)/2;
  double errx = s_channel_width;
  double erry = s_channel_width*s_nchannel_per_layer;
  if (id>49+100*layer && id<100*(layer+1)) {
    errx = s_channel_width*s_nchannel_per_layer;
    erry = s_channel_width;
  }
  return TVector3(errx,erry,0.0);
}

#endif
