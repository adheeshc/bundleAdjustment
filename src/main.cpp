#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>
#include <g2o/core/base_vertex.h>
#include <g2o/core/g2o_core_api.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/core/robust_kernel_impl.h>
#include <chrono>

#include "sophus/se3.hpp"
#include "sophus/so3.hpp"
#include "balReader.h"


struct poseAndIntrinsics {
    poseAndIntrinsics() {}
    Sophus::SO3d rotation;
    Eigen::Vector3d translation = Eigen::Vector3d::Zero();
    double focal = 0;
    double k1 = 0, k2 = 0;

    //set data
    explicit poseAndIntrinsics(double* data_address) {
        rotation = Sophus::SO3d::exp(Eigen::Vector3d(data_address[0], data_address[1], data_address[2]));
        translation = Eigen::Vector3d(data_address[3], data_address[4], data_address[5]);
        focal = data_address[6];
        k1 = data_address[7];
        k2 = data_address[8];
    }

    //set data address
    void setTo(double* data_address) {
        auto r = rotation.log();
        for (int i = 0; i < 3;i++)
            data_address[i] = r[i];
        for (int i = 0; i < 3;i++)
            data_address[i + 3] = translation[i];
        data_address[6] = focal;
        data_address[7] = k1;
        data_address[8] = k2;
    }
};

class VertexPoseAndIntrinsics : public g2o::BaseVertex<9, poseAndIntrinsics> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    VertexPoseAndIntrinsics() {};

    //override reset function
    virtual void setToOriginImpl() override {
        _estimate = poseAndIntrinsics();
    }

    //override plus operator
    virtual void oplusImpl(const double* update) override {
        _estimate.rotation = Sophus::SO3d::exp(Eigen::Vector3d(update[0], update[1], update[2])) * _estimate.rotation;
        _estimate.translation += Eigen::Vector3d(update[3], update[4], update[5]);
        _estimate.focal += update[6];
        _estimate.k1 += update[7];
        _estimate.k2 += update[8];
    }

    Eigen::Vector2d project(const Eigen::Vector3d& point) {
        Eigen::Vector3d pc = _estimate.rotation * point + _estimate.translation;
        pc = -pc / pc[2];
        double r2 = pc.squaredNorm();
        double distortion = 1.0 + r2 * (_estimate.k1 + _estimate.k2 * r2);
        return Eigen::Vector2d(_estimate.focal * distortion * pc[0], _estimate.focal * distortion * pc[1]);
    }

    //dummy read/write function
    virtual bool read(std::istream& in) {};
    virtual bool write(std::ostream& out) const {};
};

class VertexPoint : public g2o::BaseVertex<3, Eigen::Vector3d> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    VertexPoint() {};

    virtual void setToOriginImpl() override {
        _estimate = Eigen::Vector3d(0, 0, 0);
    }

    //override plus operator
    virtual void oplusImpl(const double* update) override {
        _estimate += Eigen::Vector3d(update[0], update[1], update[2]);
    }

    //dummy read/write function
    virtual bool read(std::istream& in) {};
    virtual bool write(std::ostream& out) const {};
};

class EdgeProjection : public g2o::BaseBinaryEdge<2, Eigen::Vector2d, VertexPoseAndIntrinsics, VertexPoint> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    //define error term computation
    virtual void computeError() override {
        auto v0 = (VertexPoseAndIntrinsics*)_vertices[0];
        auto v1 = (VertexPoint*)_vertices[1];
        auto proj = v0->project(v1->estimate());
        _error = proj - _measurement;
    }

    //dummy read/write function
    virtual bool read(std::istream& in) {};
    virtual bool write(std::ostream& out) const {};
};

void solveBA(BALProblem& bal_problem) {
    const int point_block_size = bal_problem.getPointBlockSize();
    const int camera_block_size = bal_problem.getCameraBlockSize();
    double* points = bal_problem.mutablePoints();
    double* cameras = bal_problem.mutableCameras();

    // pose dimension 9, landmark is 3
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<9, 3>> BlockSolverType;
    typedef g2o::LinearSolverCSparse<BlockSolverType::PoseMatrixType> LinearSolverType; // use LM
    auto solver = new g2o::OptimizationAlgorithmLevenberg(g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    optimizer.setVerbose(true);

    const double* observations = bal_problem.getObservations();
    //vertex
    std::vector<VertexPoseAndIntrinsics*> vertex_pose_intrinsics;
    std::vector<VertexPoint*> vertex_points;
    for (int i = 0; i < bal_problem.getNumCameras(); ++i) {
        VertexPoseAndIntrinsics* v = new VertexPoseAndIntrinsics();
        double* camera = cameras + camera_block_size * i;
        v->setId(i);
        v->setEstimate(poseAndIntrinsics(camera));
        optimizer.addVertex(v);
        vertex_pose_intrinsics.push_back(v);
    }
    for (int i = 0; i < bal_problem.getNumPoints(); ++i) {
        VertexPoint* v = new VertexPoint();
        double* point = points + point_block_size * i;
        v->setId(i + bal_problem.getNumCameras());
        v->setEstimate(Eigen::Vector3d(point[0], point[1], point[2]));
        v->setMarginalized(true);
        optimizer.addVertex(v);
        vertex_points.push_back(v);
    }

    //edges
    for (int i = 0; i < bal_problem.getNumObservations(); ++i) {
        EdgeProjection* edge = new EdgeProjection;
        edge->setVertex(0, vertex_pose_intrinsics[bal_problem.getCameraIndex()[i]]);
        edge->setVertex(1, vertex_points[bal_problem.getPointIndex()[i]]);
        edge->setMeasurement(Eigen::Vector2d(observations[2 * i + 0], observations[2 * i + 1]));
        edge->setInformation(Eigen::Matrix2d::Identity());
        edge->setRobustKernel(new g2o::RobustKernelHuber());
        optimizer.addEdge(edge);
    }
    optimizer.initializeOptimization();
    optimizer.optimize(20);

    // set to bal problem
    for (int i = 0; i < bal_problem.getNumCameras(); ++i) {
        double* camera = cameras + camera_block_size * i;
        auto vertex = vertex_pose_intrinsics[i];
        auto estimate = vertex->estimate();
        estimate.setTo(camera);
    }
    for (int i = 0; i < bal_problem.getNumPoints(); ++i) {
        double* point = points + point_block_size * i;
        auto vertex = vertex_points[i];
        for (int k = 0; k < 3; ++k) point[k] = vertex->estimate()[k];
    }

}

int main() {

    std::string filename = "../data/problem-16-22106-pre.txt";
    BALProblem balProblem(filename);
    balProblem.normalize();
    balProblem.perturb(0.1, 0.5, 0.5);
    balProblem.writeToPLYFile("../output/initial.ply");
    solveBA(balProblem);
    balProblem.writeToPLYFile("../output/final.ply");

    std::cout << "Done!" << std::endl;
    return 0;
}
