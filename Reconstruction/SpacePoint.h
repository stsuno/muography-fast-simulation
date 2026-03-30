#ifndef SPACE_POINT_H
#define SPACE_POINT_H

#include <vector>

#include <TVector3.h>

class SpacePoint : public TVector3 {
  public:
    SpacePoint(double x, double y, double z)
      : TVector3(x, y, z) {}
    virtual ~SpacePoint() {}

    // Position management
    void SetPosition(double x, double y, double z) { this->SetXYZ(x, y, z); }

    // Error management
    double GetErrorX() const { return m_error_x; }
    double GetErrorY() const { return m_error_y; }
    double GetErrorZ() const { return m_error_z; }
    void SetErrorX(double error_x) { m_error_x = error_x; }
    void SetErrorY(double error_y) { m_error_y = error_y; }
    void SetErrorZ(double error_z) { m_error_z = error_z; }

    // Index and Hit management
    int GetUnitIndex() const { return m_unit_index; }
    void SetUnitIndex(int unit_index) { m_unit_index = unit_index; }
    void AddHitId(int hit_id) { m_hit_ids.push_back(hit_id); }
    const std::vector<int>&GetHitIds() const { return m_hit_ids; }

  private:
    int m_unit_index;
    double m_error_x = 0.5;
    double m_error_y = 0.5;
    double m_error_z = 0.5;

    // List of raw hit IDs forming this space point
    std::vector<int> m_hit_ids;
};

#endif
