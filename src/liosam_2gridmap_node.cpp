#include <math.h>
#include <nav_msgs/OccupancyGrid.h>
#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/features/normal_3d.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_ros/transforms.h>
#include <ros/ros.h>
#include <stdint.h>
#include <tf/transform_listener.h>

#include <fstream>
#include <iostream>

#include "lio_sam/cloud_info.h"

typedef pcl::PointXYZ PointT;
typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;

using namespace std;
using namespace lio_sam;

/* 全局变量 */
PointCloud::Ptr cloud_filterd(new PointCloud);
PointCloud::Ptr world_cloud(new PointCloud);
PointCloud::Ptr key_frame_cloud(new PointCloud);
PointCloud::Ptr key_frame_poses(new PointCloud);

bool newPointCloud = false;
double cellResolution = 0.05;
double maxRange = 50.00;
bool m_filterGroundPlane = false;
double m_groundFilterDistance = 0.04;
double m_groundFilterPlaneDistance = 0.07;
double m_groundFilterAngle = 0.15;

ros::Publisher pub_map_;
ros::Publisher pub_points_;
std::string mapFrame = "map";
std::string worldFrameId = "map";
tf::TransformListener *tfListener;

template <typename T>
void operator+=(std::vector<T> &v1, const std::vector<T> &v2) {
    v1.insert(v1.end(), v2.begin(), v2.end());
}

class GridMap {
   public:
    double resolution_;
    double max_range_;
    int width_;
    int height_;
    int size_;
    double mapOriginPosX_;
    double mapOriginPosY_;
    std::string mapFrameId_;
    std::vector<signed char> Grid_;
    // maxrange check

   public:
    GridMap() {
        resolution_ = 0.01;
    }
    ~GridMap() {}
    inline bool isValidGridIndex(const int &x, const int &y) const { return ((x < width_) && (y < height_) && (x >= 0) && (y >= 0)); }
    void SetMap(int width, int height, double mapOriginPosX, double mapOriginPosY, double resolution, std::string mapFrameId) {
        resolution_ = resolution;
        width_ = width;
        height_ = height;
        mapOriginPosX_ = mapOriginPosX;
        mapOriginPosY_ = mapOriginPosY;
        size_ = width_ * height_;
        mapFrameId_ = mapFrameId;
        Grid_.resize(size_, -1);
    }

    void rayCasting(int x0, int y0, int x1, int y1);
    void cloud_to_grid(PointCloud cloud);
};

void GridMap::cloud_to_grid(PointCloud cloud) {
    for (const auto &p : cloud) {
        int x = (int)((p.x - mapOriginPosX_) / resolution_);
        int y = (int)((p.y - mapOriginPosY_) / resolution_);

        int x0 = (int)((mapOriginPosX_) / resolution_);
        int y0 = (int)((-mapOriginPosY_) / resolution_);
        int x1 = (int)((p.x - mapOriginPosX_) / resolution_);
        int y1 = (int)((p.y - mapOriginPosY_) / resolution_);

        Grid_[y * width_ + x] = 100;  //-1:unknow;100:occupied;0:free;
        rayCasting(x0, y0, x1, y1);
    }
}
/*
all other points: free on ray, occupied on endpoint
ray trace is done by Bresenham's line algorithm.
*/
void GridMap::rayCasting(int x0, int y0, int x1, int y1) {
    int dx, dy, x, y, sx, sy, n, e;

    if (x0 > x1)
        dx = x0 - x1;
    else
        dx = x1 - x0;

    if (y0 > y1)
        dy = y0 - y1;
    else
        dy = y1 - y0;

    x = x0;
    y = y0;
    n = 1 + dx + dy;
    e = static_cast<long int>(dx) - static_cast<long int>(dy);
    sx = x0 < x1 ? 1 : -1;
    sy = y0 < y1 ? 1 : -1;
    dx *= 2;
    dy *= 2;

    for (; n > 0; --n) {
        if (isValidGridIndex(x, y) && Grid_[y * width_ + x] != 100)
            Grid_[y * width_ + x] = 0;

        if (e > 0) {
            x += sx;
            e -= dy;
        } else {
            y += sy;
            e += dx;
        }
    }
}

void pub_gridmsg(GridMap map) {
    static int seq = 0;
    seq++;
    nav_msgs::OccupancyGridPtr grid(new nav_msgs::OccupancyGrid);

    grid->header.frame_id = map.mapFrameId_;
    grid->header.seq = key_frame_cloud->header.seq;
    grid->header.stamp.sec = ros::Time::now().sec;
    grid->header.stamp.nsec = ros::Time::now().nsec;
    grid->info.map_load_time = ros::Time::now();

    grid->info.resolution = cellResolution;
    grid->info.width = map.width_;
    grid->info.height = map.height_;
    grid->info.origin.position.x = map.mapOriginPosX_;
    grid->info.origin.position.y = map.mapOriginPosY_;
    grid->info.origin.position.z = 0;
    grid->info.origin.orientation.w = 1;
    grid->info.origin.orientation.x = 0;
    grid->info.origin.orientation.y = 0;
    grid->info.origin.orientation.z = 0;
    grid->data = map.Grid_;

    pub_map_.publish(grid);
    ROS_INFO_STREAM("发布栅格地图");
}

void calcSize(PointCloud::Ptr &input_cloud, double *xMax, double *yMax, double *xMin, double *yMin) {
    pcl::PointXYZ minPt, maxPt;
    pcl::getMinMax3D(*input_cloud, minPt, maxPt);
    *xMax = maxPt.x;
    *yMax = maxPt.y;
    *xMin = minPt.x;
    *yMin = minPt.y;
}

void filterGroundPlane(const PointCloud::Ptr &pc, PointCloud::Ptr &ground, PointCloud::Ptr &nonground) {
    ground->header = pc->header;
    nonground->header = pc->header;

    if (pc->size() < 50) {
        ROS_WARN("Pointcloud in OctomapServer too small, skipping ground plane extraction");
        nonground = pc;
    } else {
        // plane detection for ground plane removal:
        pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
        pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

        // Create the segmentation object and set up:
        pcl::SACSegmentation<PointT> seg;
        seg.setOptimizeCoefficients(true);
        // TODO: maybe a filtering based on the surface normals might be more robust / accurate?
        seg.setModelType(pcl::SACMODEL_PERPENDICULAR_PLANE);
        seg.setMethodType(pcl::SAC_RANSAC);
        seg.setMaxIterations(200);
        seg.setDistanceThreshold(m_groundFilterDistance);
        seg.setAxis(Eigen::Vector3f(0, 0, 1));
        seg.setEpsAngle(m_groundFilterAngle);

        PointCloud cloud_filtered(*pc);
        // Create the filtering object
        pcl::ExtractIndices<PointT> extract;
        bool groundPlaneFound = false;

        while (cloud_filtered.size() > 10 && !groundPlaneFound) {
            seg.setInputCloud(cloud_filtered.makeShared());
            seg.segment(*inliers, *coefficients);
            if (inliers->indices.size() == 0) {
                ROS_INFO("PCL segmentation did not find any plane.");

                break;
            }

            extract.setInputCloud(cloud_filtered.makeShared());
            extract.setIndices(inliers);

            if (std::abs(coefficients->values.at(3)) < m_groundFilterPlaneDistance) {
                ROS_DEBUG("Ground plane found: %zu/%zu inliers. Coeff: %f %f %f %f", inliers->indices.size(), cloud_filtered.size(),
                          coefficients->values.at(0), coefficients->values.at(1), coefficients->values.at(2), coefficients->values.at(3));
                extract.setNegative(false);
                extract.filter(*ground);

                // remove ground points from full pointcloud:
                // workaround for PCL bug:
                if (inliers->indices.size() != cloud_filtered.size()) {
                    extract.setNegative(true);
                    PointCloud cloud_out;
                    extract.filter(cloud_out);
                    *nonground += cloud_out;
                    cloud_filtered = cloud_out;
                }

                groundPlaneFound = true;
            } else {
                ROS_DEBUG("Horizontal plane (not ground) found: %zu/%zu inliers. Coeff: %f %f %f %f", inliers->indices.size(), cloud_filtered.size(),
                          coefficients->values.at(0), coefficients->values.at(1), coefficients->values.at(2), coefficients->values.at(3));
                pcl::PointCloud<PointT> cloud_out;
                extract.setNegative(false);
                extract.filter(cloud_out);
                *nonground += cloud_out;

                if (inliers->indices.size() != cloud_filtered.size()) {
                    extract.setNegative(true);
                    cloud_out.points.clear();
                    extract.filter(cloud_out);
                    cloud_filtered = cloud_out;
                } else {
                    cloud_filtered.points.clear();
                }
            }
        }
        // TODO: also do this if overall starting pointcloud too small?
        if (!groundPlaneFound) {  // no plane found or remaining points too small
            ROS_WARN("No ground plane found in scan");

            // do a rough fitlering on height to prevent spurious obstacles
            pcl::PassThrough<PointT> second_pass;
            second_pass.setFilterFieldName("z");
            second_pass.setFilterLimits(-m_groundFilterPlaneDistance, m_groundFilterPlaneDistance);
            second_pass.setInputCloud(pc->makeShared());
            second_pass.filter(*ground);
            second_pass.setFilterLimitsNegative(true);
            second_pass.filter(*nonground);
        }
    }
}

void filter_cloud(PointCloud::Ptr &input_cloud) {
    PointCloud::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZ>);
    PointCloud::Ptr pc_ground(new pcl::PointCloud<pcl::PointXYZ>);
    PointCloud::Ptr pc_nonground(new pcl::PointCloud<pcl::PointXYZ>);

    /*点云降采样，减少处理时间*/
    pcl::VoxelGrid<pcl::PointXYZ> sor;
    sor.setInputCloud(input_cloud);
    sor.setLeafSize(0.3, 0.3, 0.3);
    sor.filter(*cloud_temp);

    /*将地面点云进行滤除*/
    if (m_filterGroundPlane)
        filterGroundPlane(cloud_temp, pc_ground, pc_nonground);

    cloud_filterd = pc_nonground;
}

void pub_cloud(PointCloud::Ptr &input_cloud, std::string frameId) {
    sensor_msgs::PointCloud2 pub_pc;
    pcl::toROSMsg(*input_cloud, pub_pc);
    pub_pc.header.frame_id = frameId;
    pub_pc.header.stamp = ros::Time::now();
    pub_points_.publish(pub_pc);
}

void callback(const lio_sam::cloud_infoConstPtr &msg) {
    ROS_INFO_STREAM("Convertor节点——接收到点云");
    pcl::fromROSMsg(msg->key_frame_cloud, *key_frame_cloud);
    pcl::fromROSMsg(msg->key_frame_poses, *key_frame_poses);
    key_frame_cloud->header.seq = msg->header.seq;
    key_frame_poses->header.seq = msg->header.seq;
    std::string baseFrameId = key_frame_cloud->header.frame_id;

    filter_cloud(key_frame_cloud);

    /*点云转成全局坐标系*/
    tf::StampedTransform baseToWorldTf, worldToBaseTf;

    try {
        tfListener->waitForTransform(worldFrameId, baseFrameId, msg->key_frame_cloud.header.stamp, ros::Duration(0.2));
        tfListener->lookupTransform(worldFrameId, baseFrameId, msg->key_frame_cloud.header.stamp, baseToWorldTf);
    } catch (tf::TransformException &ex) {
        ROS_DEBUG_STREAM("Transform error of baseToWorldTf: " << ex.what() << ", quitting callback");
    }

    Eigen::Matrix4f baseToWorld;
    pcl_ros::transformAsMatrix(baseToWorldTf, baseToWorld);
    pcl::transformPointCloud(*cloud_filterd, *world_cloud, baseToWorld);

    pub_cloud(world_cloud, worldFrameId);

    GridMap current_map;
    /*计算点云的最大和最小值*/
    double xMax = 0, yMax = 0, xMin = 0, yMin = 0;
    calcSize(cloud_filterd,&xMax, &yMax, &xMin, &yMin);
    /* 确定栅格地图的长和宽 */
    int xCells = ((int)((xMax) / cellResolution)) + 1;
    int yCells = ((int)((yMax - yMin) / cellResolution)) + 1;
    cout << "地图大小：" << xCells << " " << yCells << endl;
    current_map.SetMap(xCells, yCells, 0.0, yMin, cellResolution, baseFrameId);
    current_map.cloud_to_grid(*cloud_filterd);
    pub_gridmsg(current_map);

    GridMap Accumulated_map;
}

int main(int argc, char **argv) {
    setlocale(LC_ALL, "");
    ros::init(argc, argv, "convertor_node");
    ros::NodeHandle nh("~");
    tfListener = new tf::TransformListener();

    std::string subTopic = "/lio_sam/mapping/slam_info";
    std::string pubmapTopic = "/liosam_2gridmap/map";

    nh.param("sub_topic", subTopic, subTopic);
    nh.param("pub_map", pubmapTopic, pubmapTopic);
    nh.param("frame_id", mapFrame, mapFrame);
    nh.param("resolution", cellResolution, 0.05);
    nh.param("sensor_model/max_range", maxRange, 50.00);
    nh.param("filter_ground", m_filterGroundPlane, m_filterGroundPlane);
    nh.param("ground_filter/plane_distance", m_groundFilterPlaneDistance, m_groundFilterPlaneDistance);
    nh.param("ground_filter/distance", m_groundFilterDistance, m_groundFilterDistance);
    nh.param("ground_filter/angle", m_groundFilterAngle, m_groundFilterAngle);

    pub_map_ = nh.advertise<nav_msgs::OccupancyGrid>(pubmapTopic, 1, true);
    ros::Subscriber sub = nh.subscribe<lio_sam::cloud_info>("/lio_sam/mapping/slam_info", 1, callback);
    pub_points_ = nh.advertise<sensor_msgs::PointCloud2>("/filtered_points", 10);

    ros::spin();
    return 0;
}