#ifndef MUON_TRANSPORT_H
#define MUON_TRANSPORT_H

/**
 * @file MuonTransport.h
 * @brief 多重クーロン散乱つきのミュオン輸送エンジン(散乱法用)
 *
 * Muon::CalculateEnergyLoss() は「形状の入口/出口を結ぶ直線」に沿って
 * Bethe-Bloch を積分するだけで、飛跡は最後まで直進する。散乱法では飛跡の
 * 曲がりそのものが観測量なので、この輸送器では
 *
 *   1. 現在位置・現在方向から次の物質境界までの距離を求め、
 *   2. その均質領域をスラブに分割し、
 *   3. 各スラブで Bethe-Bloch のエネルギー損失と PDG の多重散乱
 *      (角度と横変位を相関付きでサンプル)を適用しながら位置と方向を更新する
 *
 * という素直なステップ輸送を行う。真空(密度 0)の領域は境界まで一気に飛ばすので、
 * 20 m 級のジオメトリでも 1 イベントあたり十数ステップで済む。
 *
 * [多重散乱の式]
 * PDG "Passage of Particles Through Matter" Eq. 34.16:
 *   theta0 = 13.6 MeV /(beta*c*p) * z * sqrt(x/X0) * [1 + 0.038*ln(x*z^2/(X0*beta^2))]
 * スラブ分割すると sqrt(x/X0) は自然に二乗和で復元されるが、対数項はスラブごとに
 * 評価すると全体の theta0 を過小評価する。そこで **sqrt はスラブ厚、対数は領域全体の
 * 厚み**で評価する。こうすると全スラブを二乗和した結果が、領域を一括で評価した
 * Highland 式と(beta, p が一定の極限で)厳密に一致する。
 *
 * 横変位は PDG Eq. 34.17-34.20 の相関付きサンプリング:
 *   y_plane     = z1*s*theta0/sqrt(12) + z2*s*theta0/2
 *   theta_plane = z2*theta0                      (z1, z2 は独立な標準正規乱数)
 * を、進行方向に直交する 2 軸それぞれに独立に適用する。
 *
 * [使い方]
 *   MuonTransport transport(seed);
 *   transport.SetShapes(absorber_shapes);   // 検出器は入れない(下記の注意参照)
 *   transport.SetStopZ(-500.0);
 *   MuonTransport::Result result = transport.Transport(muon);
 *   MuonTransport::State state;
 *   result.StateAtZ(detector_z, state);     // 検出器の高さでの位置と方向
 *
 * [注意]
 * SetShapes() に渡す形状は「物質」だけにすること。SciBar100x100 のような
 * CompositeShape は子形状 800 個ぶんの HitPoints を回すため、境界探索のたびに
 * 数千回の面交差計算が走って極端に遅くなる。検出器は密度 0 で扱い、飛跡が
 * 直線である区間の State から IsHitList() を呼ぶのが正しい使い方。
 */

#include <cmath>
#include <limits>
#include <vector>

#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TError.h>

#include "Simulation/Muon.h"
#include "Simulation/Shape.h"

class MuonTransport {
  public:
    // 飛跡上の 1 点。direction はこの点から「出ていく」方向の単位ベクトル。
    struct State {
      TVector3 position;
      TVector3 direction;
      double energy = 0.0;  // total energy [MeV]
    };

    struct Result {
      bool stopped = false;        // 物質中で静止した
      bool escaped = false;        // 下向きでなくなった等で輸送を打ち切った
      double initial_energy = 0.0; // [MeV] total
      double final_energy = 0.0;   // [MeV] total
      double path_length = 0.0;    // 実際に進んだ幾何学的距離 [cm]
      double absorber_path = 0.0;  // うち密度のある物質中の距離 [cm]
      double grammage = 0.0;       // integral rho ds [g/cm2]
      double x0_traversed = 0.0;   // integral ds/X0 [radiation lengths]
      double deflection = 0.0;     // 入射方向と射出方向のなす角 [rad]
      bool has_absorber = false;   // 密度のある物質を通ったか
      TVector3 absorber_entry;     // 最初に物質へ入った点
      TVector3 absorber_exit;      // 最後に物質から出た点
      std::vector<State> states;

      // z = z_target における位置・方向・エネルギーを線形補間で返す。
      // 飛跡がその z を通っていなければ false。
      bool StateAtZ(double z_target, State& state) const;

      const State& EntryState() const { return states.front(); }
      const State& ExitState() const { return states.back(); }
    };

    explicit MuonTransport(UInt_t random_seed = 1);

    // 輸送に使う物質形状。密度 0 の形状は素通りする。
    void SetShapes(const std::vector<Shape*>& shapes) { m_shapes = shapes; }
    void SetStopZ(double z) { m_stop_z = z; }
    void SetScattering(bool enable) { m_scattering = enable; }
    void SetEnergyLoss(bool enable) { m_energy_loss = enable; }
    // 均質領域あたりの最小スラブ数と、スラブ厚の上限 [cm]。
    // 実際のスラブ数は max(min_slabs, ceil(region/max_slab_length)) を
    // s_max_slabs_per_region で頭打ちにしたもの。
    void SetSlabsPerRegion(int n) { m_slabs_per_region = (n > 0) ? n : 1; }
    void SetMaxSlabLength(double length_cm) { m_max_slab_length = length_cm; }
    void SetMaxSteps(int n) { m_max_steps = n; }

    Result Transport(const TVector3& start, const TVector3& direction, double total_energy);
    // Muon の始点・終点・エネルギーから輸送する簡便版。
    Result Transport(const Muon& muon);

    // State を「その点を通る直線の Muon」に書き込む。検出器の IsHitList() に
    // 渡すための補助。half_length は始点/終点を State から前後に伸ばす距離 [cm]。
    static void SetMuonFromState(Muon& muon, const State& state, double half_length = 1.0e4);

    // Highland 式の theta0 [rad]。slab_x0 はスラブ厚 / X0、total_x0 は対数項に
    // 使う領域全体の厚み / X0(0 以下なら slab_x0 を使う)。
    static double HighlandTheta0(double total_energy, double slab_x0, double total_x0);

    static constexpr double s_density_floor = 1.0e-6;  // これ以下は真空扱い [g/cm3]
    static constexpr double s_push = 1.0e-4;           // 境界をまたぐための微小前進 [cm]
    static constexpr double s_tolerance = 1.0e-6;      // [cm]
    static constexpr int s_max_slabs_per_region = 2000;

  private:
    struct LocalMaterial {
      bool found = false;
      double density = 0.0;
      double z_over_a = 0.0;
      double mean_excitation_energy = 0.0;
      double radiation_length = 0.0;  // X0 [g/cm2]
      Shape::SternheimerParameters sternheimer{};
    };

    // Muon::CalculateEnergyLoss と同じ優先度規則(priority が小さいほど優先)で
    // 点における物質を決める。
    LocalMaterial MaterialAt(const TVector3& point) const;
    // 現在位置から進行方向にある最も近い物質境界までの距離 [cm]。
    // 境界が無ければ無限大を返す。
    double DistanceToBoundary(const TVector3& position, const TVector3& direction);
    // 1 スラブぶんの多重散乱を適用し、position と direction を更新する。
    void ApplyScattering(TVector3& position, TVector3& direction, double total_energy,
                         double slab_length, double x0_cm, double region_x0);
    // PDG Eq. 34.17-34.20 の相関付き (横変位, 角度) を 1 平面ぶんサンプルする。
    void SamplePlane(double theta0, double slab_length, double& y_plane, double& theta_plane);

    std::vector<Shape*> m_shapes;
    TRandom3 m_rnd;
    Muon m_probe;  // 境界探索用の使い回しトラック(乱数は使わない)

    double m_stop_z = -1.0e4;
    bool m_scattering = true;
    bool m_energy_loss = true;
    int m_slabs_per_region = 10;
    double m_max_slab_length = 50.0;  // [cm]
    int m_max_steps = 10000;

    bool m_warned_incomplete_material = false;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

inline MuonTransport::MuonTransport(UInt_t random_seed)
  : m_rnd(random_seed), m_probe(nullptr) {}

inline bool MuonTransport::Result::StateAtZ(double z_target, State& state) const {
  if (states.size() < 2) return false;
  // 折れ点ちょうどの z では上下 2 つのセグメントが該当する。散乱体の出口で
  // 「散乱前の方向」を返してしまわないよう、下側(=散乱後)のセグメントを優先する。
  // そのため後ろから走査して最初に見つかったものを返す。
  for (size_t j = states.size() - 1; j >= 1; --j) {
    size_t i = j - 1;
    double z_high = states[i].position.Z();
    double z_low = states[i + 1].position.Z();
    if (z_target > z_high + s_tolerance || z_target < z_low - s_tolerance) continue;

    double dz = z_high - z_low;
    double fraction = (dz > s_tolerance) ? (z_high - z_target) / dz : 0.0;
    state.position = states[i].position + fraction * (states[i + 1].position - states[i].position);
    state.direction = states[i].direction;
    state.energy = states[i].energy + fraction * (states[i + 1].energy - states[i].energy);
    return true;
  }
  return false;
}

inline void MuonTransport::SetMuonFromState(Muon& muon, const State& state, double half_length) {
  muon.SetTrackEndpoints(state.position - half_length * state.direction,
                         state.position + half_length * state.direction);
  muon.SetEnergy(state.energy);
}

inline double MuonTransport::HighlandTheta0(double total_energy, double slab_x0, double total_x0) {
  if (slab_x0 <= 0.0) return 0.0;
  double momentum2 = total_energy * total_energy - Muon::kMuonMass * Muon::kMuonMass;
  if (momentum2 <= 0.0) return 0.0;
  double momentum = TMath::Sqrt(momentum2);          // pc [MeV]
  double beta = momentum / total_energy;
  if (beta <= 0.0) return 0.0;

  // 対数項は領域全体の厚みで評価する(ファイル冒頭のコメント参照)。
  double log_thickness = (total_x0 > 0.0) ? total_x0 : slab_x0;
  double bracket = 1.0 + 0.038 * TMath::Log(log_thickness / (beta * beta));
  if (bracket < 0.0) bracket = 0.0;

  return 13.6 / (beta * momentum) * TMath::Sqrt(slab_x0) * bracket;
}

inline MuonTransport::LocalMaterial MuonTransport::MaterialAt(const TVector3& point) const {
  LocalMaterial material;
  int highest_priority = Shape::s_invalid_priority;
  for (auto shape : m_shapes) {
    if (!shape->IsInside(point)) continue;
    if (shape->GetPriority() >= highest_priority) continue;
    highest_priority = shape->GetPriority();
    material.found = true;
    material.density = shape->DensityAt(point);
    material.z_over_a = shape->ZOverAAt(point);
    material.mean_excitation_energy = shape->MeanExcitationEnergyAt(point);
    material.radiation_length = shape->RadiationLengthAt(point);
    material.sternheimer = shape->SternheimerParametersAt(point);
  }
  return material;
}

inline double MuonTransport::DistanceToBoundary(const TVector3& position, const TVector3& direction) {
  m_probe.SetTrackEndpoints(position, position + 100.0 * direction);

  double nearest = std::numeric_limits<double>::infinity();
  for (auto shape : m_shapes) {
    // HitPoints() は飛跡を無限直線として扱うので、進行方向前方だけを残す。
    for (const TVector3& point : shape->HitPoints(m_probe)) {
      double distance = (point - position).Dot(direction);
      if (distance > s_push && distance < nearest) nearest = distance;
    }
  }
  return nearest;
}

inline void MuonTransport::SamplePlane(double theta0, double slab_length,
                                       double& y_plane, double& theta_plane) {
  double z1 = m_rnd.Gaus(0.0, 1.0);
  double z2 = m_rnd.Gaus(0.0, 1.0);
  y_plane = z1 * slab_length * theta0 / TMath::Sqrt(12.0) + z2 * slab_length * theta0 / 2.0;
  theta_plane = z2 * theta0;
}

inline void MuonTransport::ApplyScattering(TVector3& position, TVector3& direction,
                                           double total_energy, double slab_length,
                                           double x0_cm, double region_x0) {
  double theta0 = HighlandTheta0(total_energy, slab_length / x0_cm, region_x0);
  if (theta0 <= 0.0) {
    position += direction * slab_length;
    return;
  }

  // 進行方向に直交する正規直交基底
  TVector3 axis_u = direction.Orthogonal().Unit();
  TVector3 axis_v = direction.Cross(axis_u).Unit();

  double offset_u = 0.0, angle_u = 0.0;
  double offset_v = 0.0, angle_v = 0.0;
  SamplePlane(theta0, slab_length, offset_u, angle_u);
  SamplePlane(theta0, slab_length, offset_v, angle_v);

  position += direction * slab_length + axis_u * offset_u + axis_v * offset_v;
  direction = (direction + axis_u * TMath::Tan(angle_u) + axis_v * TMath::Tan(angle_v)).Unit();
}

inline MuonTransport::Result MuonTransport::Transport(const Muon& muon) {
  TVector3 start = muon.GetStartPoint();
  TVector3 direction = muon.GetEndPoint() - start;
  return Transport(start, direction, muon.GetEnergyMuon());
}

inline MuonTransport::Result MuonTransport::Transport(const TVector3& start,
                                                      const TVector3& direction,
                                                      double total_energy) {
  Result result;
  result.initial_energy = total_energy;
  result.final_energy = total_energy;

  TVector3 position = start;
  TVector3 heading = direction.Unit();
  double energy = total_energy;

  const TVector3 entry_direction = heading;
  result.states.push_back({position, heading, energy});

  if (heading.Z() >= 0.0) {
    // 上向き/水平の飛跡は想定外。呼び出し側の座標系の取り違えを黙って通さない。
    Error("MuonTransport::Transport", "Track is not downward (dz = %.3f)", heading.Z());
    result.escaped = true;
    return result;
  }

  int steps = 0;
  while (steps < m_max_steps && position.Z() > m_stop_z) {
    ++steps;

    if (heading.Z() >= 0.0) { result.escaped = true; break; }

    double distance_to_boundary = DistanceToBoundary(position, heading);
    double distance_to_stop = (m_stop_z - position.Z()) / heading.Z();
    double region = TMath::Min(distance_to_boundary, distance_to_stop);
    if (!std::isfinite(region) || region <= s_push) break;

    LocalMaterial material = MaterialAt(position + heading * (0.5 * region));

    if (!material.found || material.density <= s_density_floor) {
      // 真空(または実質真空)の領域は境界まで一気に飛ばす。
      position += heading * (region + s_push);
      result.path_length += region;
      result.states.push_back({position, heading, energy});
      continue;
    }

    // Bethe-Bloch に必要な物性が揃っていない形状は、静かに間違った答えを返すより
    // 警告してエネルギー損失だけ止める(散乱は X0 があれば続行できる)。
    bool can_lose_energy = m_energy_loss &&
                           material.z_over_a > 0.0 && material.mean_excitation_energy > 0.0;
    if (m_energy_loss && !can_lose_energy && !m_warned_incomplete_material) {
      Warning("MuonTransport::Transport",
              "Shape without full Bethe-Bloch properties (<Z/A> or I missing): energy loss skipped");
      m_warned_incomplete_material = true;
    }

    double x0_cm = (material.radiation_length > 0.0) ? material.radiation_length / material.density : 0.0;
    double region_x0 = (x0_cm > 0.0) ? region / x0_cm : 0.0;

    int n_slabs = m_slabs_per_region;
    if (m_max_slab_length > 0.0) {
      int required = static_cast<int>(std::ceil(region / m_max_slab_length));
      if (required > n_slabs) n_slabs = required;
    }
    if (n_slabs > s_max_slabs_per_region) n_slabs = s_max_slabs_per_region;
    if (n_slabs < 1) n_slabs = 1;
    double slab_length = region / n_slabs;

    for (int i_slab = 0; i_slab < n_slabs; ++i_slab) {
      // 散乱で飛跡が曲がると、入口で測った領域から途中で抜けることがある。
      // その場合はここで抜けて外側のループで領域を測り直す。
      LocalMaterial local = MaterialAt(position + heading * (0.5 * slab_length));
      if (!local.found || local.density <= s_density_floor) break;

      if (!result.has_absorber) {
        result.has_absorber = true;
        result.absorber_entry = position;
      }

      if (can_lose_energy) {
        double dedx = Muon::BetheBlochdEdx(local.z_over_a, local.mean_excitation_energy,
                                           local.sternheimer, energy);
        double loss = dedx * local.density * slab_length;
        double kinetic_energy = energy - Muon::kMuonMass;
        if (loss >= kinetic_energy) {
          energy = Muon::kMuonMass;
          result.stopped = true;
        } else {
          energy -= loss;
        }
      }

      result.grammage += local.density * slab_length;
      result.absorber_path += slab_length;
      double local_x0_cm = (local.radiation_length > 0.0) ? local.radiation_length / local.density : 0.0;
      if (local_x0_cm > 0.0) result.x0_traversed += slab_length / local_x0_cm;

      if (m_scattering && local_x0_cm > 0.0 && !result.stopped) {
        ApplyScattering(position, heading, energy, slab_length, local_x0_cm, region_x0);
      } else {
        position += heading * slab_length;
      }
      result.path_length += slab_length;
      result.absorber_exit = position;
      result.states.push_back({position, heading, energy});

      if (result.stopped || position.Z() <= m_stop_z) break;
    }

    if (result.stopped) break;
    position += heading * s_push;
  }

  result.final_energy = energy;
  double cos_deflection = entry_direction.Dot(heading);
  if (cos_deflection > 1.0) cos_deflection = 1.0;
  if (cos_deflection < -1.0) cos_deflection = -1.0;
  result.deflection = TMath::ACos(cos_deflection);
  result.states.back().energy = energy;

  return result;
}

#endif
