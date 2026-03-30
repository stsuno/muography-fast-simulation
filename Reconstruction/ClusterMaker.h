#ifndef CLUSTER_MAKER_H
#define CLUSTER_MAKER_H

#include <map>
#include <vector>
#include <utility>
#include <algorithm>

#include "Reconstruction/SpacePoint.h"

template <typename VariousDetector>

class ClusterMaker {
  public:
    ClusterMaker() {}
    ClusterMaker(VariousDetector detector):m_detector(detector),m_model(CombinedLocal) {}
    virtual ~ClusterMaker() {}

    enum ClusterModel{Global,CombinedLocal,GlobalAndTime};

    void SetClusterModel(ClusterMaker::ClusterModel model) { m_model=model; }

    // Processes a list of raw hit IDs and returns a collection of SpacePoints.
    std::vector<SpacePoint> Execute(const std::vector<int>& hit_ids) const {
      if (hit_ids.empty()) return {};

      if (m_model==Global) {
        return SetGlobal(hit_ids);
      }
      else if (m_model==CombinedLocal) {
        return SetCombinedLocal(hit_ids);
      }
      return {};
    }

    std::vector<SpacePoint> Execute(const std::vector<std::pair<int,double>>& hit_ids) const {
      if (hit_ids.empty()) return {};

      if (m_model==GlobalAndTime) {
        return SetGlobalAndTime(hit_ids);
      }
      return {};
    }

  private:
    VariousDetector m_detector;
    ClusterModel m_model;

    std::vector<SpacePoint> SetGlobal(const std::vector<int>& hit_ids) const {
      std::vector<SpacePoint> space_points;
      for (auto& [layer, clusters] : Clustering(hit_ids)) {
        for (const auto& cluster : clusters) {
          // Calculate average position (Center of Gravity)
          TVector3 pos(0.0,0.0,0.0),err(0.0,0.0,0.0);
          for (int id : cluster) { 
            pos += m_detector->GetGlobalPosition(id); 
            err += m_detector->GetGlobalPositionError(id);
          }
          if (cluster.size()) { 
            pos *= 1.0/cluster.size();
            err *= 1.0/cluster.size();
          }
          SpacePoint space_point(pos.X(),pos.Y(),pos.Z());
          space_point.SetUnitIndex(layer);
          space_point.SetErrorX(err.X());
          space_point.SetErrorY(err.Y());
          space_point.SetErrorZ(err.Z());

          // TODO: Implement ghost rejection logic
          space_points.push_back(space_point);
        }
      }
      return space_points;
    }

    std::vector<SpacePoint> SetCombinedLocal(const std::vector<int>& hit_ids) const {
      // Clustering for sequential IDs
      std::map<int, std::vector<std::vector<int>>> clusters_by_layer = Clustering(hit_ids);

      // Combine X and Y layers to form 3D SpacePoints
      int nlayer = m_detector->GetNumLayers();
      int total_layer_channel = m_detector->GetNumIDs()/nlayer;
      std::vector<SpacePoint> space_points;
      for (int i=0; i<nlayer; i++) {
        const auto& x_clusters = clusters_by_layer[2*i];    // Even: X layer
        const auto& y_clusters = clusters_by_layer[2*i+1];  //  Odd: Y layer

        // Representative Z from the first ID of the X layer
        double z = m_detector->GetLocalPosition(i*total_layer_channel).Z();
        for (const auto& x_cluster : x_clusters) {
          for (const auto& y_cluster : y_clusters) {

            // Calculate average position (Center of Gravity)
            double xpos(0.0),ypos(0.0);
            for (int id : x_cluster) { xpos+=m_detector->GetLocalPosition(id).X(); }
            for (int id : y_cluster) { ypos+=m_detector->GetLocalPosition(id).X(); }
            if (x_cluster.size()) { xpos*=1.0/x_cluster.size(); }
            if (y_cluster.size()) { ypos*=1.0/y_cluster.size(); }

            SpacePoint space_point(xpos, ypos, z);
            space_point.SetUnitIndex(i);
            space_point.SetErrorX(x_cluster.size());
            space_point.SetErrorY(y_cluster.size());
            space_point.SetErrorZ(0.1);

            // TODO: Implement ghost rejection logic
            space_points.push_back(space_point);
          }
        }
      }
      return space_points;
    }

    std::map<int, std::vector<std::vector<int>>> Clustering(const std::vector<int>& hit_ids) const {
      // Group hits by layer: map<layer_index, vector<hit_id>>
      int nlayer = m_detector->GetNumLayers();
      int total_layer_channel = m_detector->GetNumIDs()/nlayer;
      std::map<int, std::vector<int>> hits_by_layer;
      for (int id : hit_ids) {
        hits_by_layer[m_detector->GetLayerIndex(id)].push_back(id);
      }

      // Perform clustering for each layer: map<layer_index, vector<contributing_ids>>
      std::map<int, std::vector<std::vector<int>>> clusters_by_layer;
      for (auto& [layer, ids] : hits_by_layer) {
        std::sort(ids.begin(), ids.end());  // Ensure IDs are sequential

        std::vector<int> current_cluster_ids;
        for (size_t i = 0; i < ids.size(); ++i) {
          current_cluster_ids.push_back(ids[i]);

          // Check if next ID is not continuous or it's the last ID
          if (i+1==ids.size() || ids[i+1]!=ids[i]+1) {
            clusters_by_layer[layer].push_back(current_cluster_ids);
            current_cluster_ids.clear();
          }
        }
      }
      return clusters_by_layer;
    }

    std::vector<SpacePoint> SetGlobalAndTime(const std::vector<std::pair<int,double>>& hit_ids) const {
      std::vector<SpacePoint> space_points;

//      std::cout << "START " << hit_ids.size() << std::endl;

      std::vector<std::pair<TVector3,TVector3>> points;
      for (const auto& hit : hit_ids) {
        int fiberid = hit.first;
        double time = hit.second;
        double speed = m_detector->GetFiber(fiberid)->GetSignalSpeed();
        double element = m_detector->GetFiber(fiberid)->GetElementLength();
        int hitID = speed*time/element;

        double resolution = m_detector->GetFiber(fiberid)->GetTimingResolution();
        int minusID    = speed*(time-resolution)/element;
        int positiveID = speed*(time+resolution)/element;
        if (minusID<0) { minusID=0; }
        int maxelement = m_detector->GetFiber(fiberid)->GetNumElement();
        if (positiveID>maxelement-1) { positiveID=maxelement-1; }

        TVector3 pos  = m_detector->GetFiber(fiberid)->GetHitIDCenter(hitID);
        TVector3 errm = m_detector->GetFiber(fiberid)->GetHitIDCenter(minusID);
        TVector3 errp = m_detector->GetFiber(fiberid)->GetHitIDCenter(positiveID);
        TVector3 err  = 0.5*(errm+errp);

        points.push_back(std::make_pair(pos,err));
      }

      // grouping
      int nlayer = 0;
      int chk[40]={0};
      for (size_t i=0; i<points.size(); i++) {
        if (chk[i]==0) {
          int icnt = 1;
          TVector3 pos = points[i].first;
          TVector3 err((points[i].second.X())*(points[i].second.X()),
                       (points[i].second.Y())*(points[i].second.Y()),
                       (points[i].second.Z())*(points[i].second.Z()));

          for (size_t j=i+1; j<points.size(); j++) {
            double distance = (points[i].first-points[j].first).Mag();
            if (distance<1.0) {
              pos += points[j].first;
              err += TVector3((points[j].second.X())*(points[j].second.X()),
                              (points[j].second.Y())*(points[j].second.Y()),
                              (points[j].second.Z())*(points[j].second.Z()));
              icnt++;
              chk[j] = 1;
            }
          }
          pos *= 1.0/icnt;

          SpacePoint space_point(pos.X(),pos.Y(),pos.Z());
          space_point.SetUnitIndex(nlayer);
          space_point.SetErrorX(TMath::Sqrt(err.X()));
          space_point.SetErrorY(TMath::Sqrt(err.Y()));
          space_point.SetErrorZ(TMath::Sqrt(err.Z()));
          nlayer++;

          // TODO: Implement ghost rejection logic
          space_points.push_back(space_point);
        }
      }
      return space_points;
    }

};

#endif
