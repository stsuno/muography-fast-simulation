#ifndef SAKURAJIMAMWPC_H
#define SAKURAJIMAMWPC_H

#include "Simulation/CompositeShape.h"
#include "Simulation/RectangularPrism.h"

class SakurajimaMWPC : public CompositeShape {
  public:
    SakurajimaMWPC() { };
    virtual ~SakurajimaMWPC() = default;

    static constexpr int    s_nlayers            = 8;
    static constexpr int    s_nchannel_per_layer = 64;
    static constexpr int    s_total_channel      = s_nlayers*s_nchannel_per_layer*2;
    static constexpr double s_channel_width      = 1.25;
    static constexpr double s_offset             = -0.5*s_channel_width*(s_nchannel_per_layer-1.0);

    int GetNumLayers() const          { return s_nlayers; }
    int GetNumIDs() const             { return s_total_channel; }
    int GetLayerIndex(int id) const   { return id/s_nchannel_per_layer; }
    int GetChannelIndex(int id) const { return id%s_nchannel_per_layer; }

    void SetMMOS(int mmos) {
      if (mmos==1) {
        m_localzposition[0] =  0.0;
        m_localzposition[1] =  321.0/12.0;
        m_localzposition[2] =  576.0/12.0;
        m_localzposition[3] =  866.0/12.0;
        m_localzposition[4] = 1173.0/12.0;
        m_localzposition[5] = 1471.0/12.0;
        m_localzposition[6] = 1724.0/12.0;
        m_localzposition[7] = 2059.0/12.0;
      }
      else if (mmos==2) {
        m_localzposition[0] =  0.0;
        m_localzposition[1] =  309.0/12.0;
        m_localzposition[2] =  579.0/12.0;
        m_localzposition[3] =  871.0/12.0;
        m_localzposition[4] = 1169.0/12.0;
        m_localzposition[5] = 1470.0/12.0;
        m_localzposition[6] = 1731.0/12.0;
        m_localzposition[7] = 2063.0/12.0;
      }
      else if (mmos==8) {
        m_localzposition[0] =  0.0;
        m_localzposition[1] =  263.0/12.0;
        m_localzposition[2] =  524.0/12.0;
        m_localzposition[3] =  805.0/12.0;
        m_localzposition[4] = 1123.0/12.0;
        m_localzposition[5] = 1399.0/12.0;
        m_localzposition[6] = 1689.0/12.0;
        m_localzposition[7] = 1952.0/12.0;
      }
      else {
        std::cout << "ERROR::SakurajimaMWPC:SetMMOS incorrect mmos: " << mmos << std::endl;
        abort();
      }
    }
    void BuildMWPC(double cx, double cy, double cz);

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
    std::vector<RectangularPrism*> m_sensewire;
    TVector3 m_detectorCenter;
    double m_localzposition[8];
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline void SakurajimaMWPC::BuildMWPC(double cx, double cy, double cz) {

  if (m_localzposition[7]<1.0) {
    std::cout << "ERROR::SakurajimaMWPC::BuildMWPC SetMMOS should be called in advance." << std::endl;
    abort();
  }

  m_sensewire.resize(s_total_channel,nullptr);
  m_detectorCenter = TVector3(cx,cy,cz);

  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) { 
      m_sensewire[s_nchannel_per_layer*2*i+j] 
        = new RectangularPrism(0.0, 0.0, 0.0, 
                               s_channel_width, s_nchannel_per_layer*s_channel_width, s_channel_width); 
    }
    for (int j=0; j<s_nchannel_per_layer; j++) { 
      m_sensewire[s_nchannel_per_layer*2*i+j+s_nchannel_per_layer] 
        = new RectangularPrism(0.0, 0.0, 0.0, 
                               s_nchannel_per_layer*s_channel_width, s_channel_width, s_channel_width); 
    }
  }

  for (int i=0; i<s_nlayers; i++) {
    for (int j=0; j<s_nchannel_per_layer; j++) {
      m_sensewire[s_nchannel_per_layer*2*i+j]                       -> MoveXYZ(s_offset+s_channel_width*j,0.0,m_localzposition[i]);
      m_sensewire[s_nchannel_per_layer*2*i+j+s_nchannel_per_layer]  -> MoveXYZ(0.0,s_offset+s_channel_width*j,m_localzposition[i]);
    }
  }

  for (int i=0; i<s_total_channel; i++) {
    m_sensewire[i] -> MoveXYZ(cx,cy,cz);
    m_sensewire[i] -> SetUniqueID(i);
    AddShape(m_sensewire[i]);
  }
}

inline void SakurajimaMWPC::MoveXYZ(double dx, double dy, double dz) {
  TVector3 displacement(dx, dy, dz);
  m_detectorCenter += displacement;
  for (int i=0; i<s_total_channel; i++) {
    std::vector<TVector3> vtx = m_sensewire[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex += displacement;
      m_sensewire[i]->UpdateVertex(j, vertex);
    }
    m_sensewire[i]->ComputeCenter();
    m_sensewire[i]->ComputeSurfaceEquation();
  }
}

inline void SakurajimaMWPC::RotatePhi(double phi) {
  for (int i=0; i<s_total_channel; i++) {
    std::vector<TVector3> vtx = m_sensewire[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex -= m_detectorCenter;
      vertex.RotateZ(phi);
      vertex += m_detectorCenter;
      m_sensewire[i]->UpdateVertex(j, vertex);
    }
    m_sensewire[i]->ComputeCenter();
    m_sensewire[i]->ComputeSurfaceEquation();
  }
}

inline void SakurajimaMWPC::RotateTheta(double theta) {
  for (int i=0; i<s_total_channel; i++) {
    std::vector<TVector3> vtx = m_sensewire[i]->GetVertices();
    for (size_t j=0; j<vtx.size(); j++) {
      TVector3 vertex = vtx[j];
      vertex -= m_detectorCenter;
      vertex.RotateX(theta);
      vertex += m_detectorCenter;
      m_sensewire[i]->UpdateVertex(j, vertex);
    }
    m_sensewire[i]->ComputeCenter();
    m_sensewire[i]->ComputeSurfaceEquation();
  }
}

inline std::vector<int> SakurajimaMWPC::IsHitList(const Muon& muon) {
  std::vector<int> hitlist;
  for (int i=0; i<s_total_channel; i++) {
    if (m_sensewire[i]->IsHit(muon)) {
      hitlist.push_back(i);
    }
  }
  return hitlist;
}

inline TVector3 SakurajimaMWPC::GetLocalPosition(int id) const {
  // Calculation logic based on the strip index
  int layer = GetLayerIndex(id)/2;
  return TVector3(s_offset+s_channel_width*GetChannelIndex(id),0.0,m_localzposition[layer]);
}

inline TVector3 SakurajimaMWPC::GetGlobalPosition(int id) const {
  return m_sensewire[id]->GetCenter();
}

inline TVector3 SakurajimaMWPC::GetGlobalPositionError(int id) const {
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
