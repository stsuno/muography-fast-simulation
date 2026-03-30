#ifndef HOUSE_H
#define HOUSE_H

#include "Simulation/CompositeShape.h"
#include "Simulation/RectangularPrism.h"
#include "Simulation/TriangularPrism.h"

class House : public CompositeShape {
  public:
    House(double cx, double cy, double cz,
        double width, double depth, double height);
    virtual ~House() = default;

    void SetRoofColor(Color_t color);

  private:
    RectangularPrism* m_body;
    RectangularPrism* m_base;
    RectangularPrism* m_door;
    TriangularPrism* m_roof;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline House::House(double cx, double cy, double cz,
    double width, double depth, double height) {

  double base_h = 30.0;
  double roof_h = height * 0.3;
  double body_h = height * 0.7;
  double overhang = 40.0;

  m_base = new RectangularPrism(cx, cy, cz + 0.5*base_h,
      width, depth, base_h);
  m_body = new RectangularPrism(cx, cy, cz + base_h + 0.5*body_h,
      0.95*width, 0.95*depth, body_h);
  m_roof = new TriangularPrism(cx, cy, cz + base_h + body_h + roof_h/3,
      width + overhang, roof_h, depth + 0.5*overhang, false);
  m_roof->RotateTheta(0.5*TMath::Pi());

  double door_w = 90.0;
  double door_h = 180.0;
  double door_t = 5.0;
  m_door = new RectangularPrism(cx, cy - 0.5*depth, cz + base_h + 0.5*door_h,
      door_w, door_t, door_h);

  Color_t wall_color = TColor::GetColor(253, 245, 230);

  m_base->SetFillColor(kGray);
  m_body->SetFillColor(wall_color);
  m_roof->SetFillColor(kOrange + 3);
  m_door->SetFillColor(kOrange - 8);

  m_base->SetLineColor(kGray + 2);
  m_body->SetLineColor(kGray + 2);
  m_roof->SetLineColor(kGray + 2);
  m_door->SetLineColor(kGray + 2);  

  AddShape(m_base);
  AddShape(m_body);
  AddShape(m_roof);
  AddShape(m_door);
}

inline void House::SetRoofColor(Color_t color) {
  m_roof->SetFillColor(color);
}

#endif
