#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>


#include "balReader.h"
#include "rotation.h"
#include "random.h"

typedef Eigen::Map<Eigen::VectorXd> VectorRef;
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;


template<typename T>
void FscanfOrDie(FILE* fptr, const char* format, T* value) {
    int num_scanned = fscanf(fptr, format, value);
    if (num_scanned != 1) {
        std::cerr << "Invalid UW data file. ";
    }
}

void PerturbPoint3(const double sigma, double* point) {
    for (int i = 0; i < 3; ++i)
        point[i] += RandNormal() * sigma;
}

double Median(std::vector<double>* data) {
    int n = data->size();
    std::vector<double>::iterator midPoint = data->begin() + n / 2;
    std::nth_element(data->begin(), midPoint, data->end());
    return *midPoint;
}

BALProblem::BALProblem(const std::string& filename, bool useQuaternion) {
    FILE* fptr = fopen(filename.c_str(), "r");
    char buffer[255];
    if (fptr == NULL) {
        std::cerr << "Error: unable to open file " << filename;
        return;
    };


    // if (fgets(buffer, sizeof(buffer), fptr) != NULL) {
    //     printf("Line read from file: %s", buffer);
    // }
    // This wil die horribly on invalid files. Them's the breaks.
    FscanfOrDie(fptr, "%d", &_numCameras);
    FscanfOrDie(fptr, "%d", &_numPoints);
    FscanfOrDie(fptr, "%d", &_numObservations);

    _numParameters = 9 * _numCameras + 3 * _numPoints;

    std::cout << "Header: " << _numCameras << " " << _numPoints << " " << _numObservations << " " << _numParameters << "\n";

    _pointIndex = new int[_numObservations];
    _cameraIndex = new int[_numObservations];
    _observations = new double[2 * _numObservations];
    _parameters = new double[_numParameters];

    for (int i = 0; i < _numObservations; ++i) {
        FscanfOrDie(fptr, "%d", _cameraIndex + i);
        FscanfOrDie(fptr, "%d", _pointIndex + i);
        for (int j = 0; j < 2; ++j) {
            FscanfOrDie(fptr, "%lf", _observations + 2 * i + j);
        }
    }

    for (int i = 0; i < _numParameters; ++i) {
        FscanfOrDie(fptr, "%lf", _parameters + i);
    }


    fclose(fptr);

    _useQuaternions = useQuaternion;
    if (_useQuaternions) {
        // Switch the angle-axis rotations to quaternions.
        _numParameters = 10 * _numCameras + 3 * _numPoints;
        double* quaternionParameters = new double[_numParameters];
        double* originalCursor = _parameters;
        double* quaternionCursor = quaternionParameters;
        for (int i = 0; i < _numCameras; ++i) {
            AngleAxisToQuaternion(originalCursor, quaternionCursor);
            quaternionCursor += 4;
            originalCursor += 3;
            for (int j = 4; j < 10; ++j) {
                *quaternionCursor++ = *originalCursor++;
            }
        }
        // Copy the rest of the points.
        for (int i = 0; i < 3 * _numPoints; ++i) {
            *quaternionCursor++ = *originalCursor++;
        }
        // Swap in the quaternion parameters.
        delete[]_parameters;
        _parameters = quaternionParameters;
    }
}

void BALProblem::writeToFile(const std::string& filename) const {
    FILE* fptr = fopen(filename.c_str(), "w");

    if (fptr == NULL) {
        std::cerr << "Error: unable to open file " << filename;
        return;
    }

    fprintf(fptr, "%d %d %d %d\n", _numCameras, _numCameras, _numPoints, _numObservations);

    for (int i = 0; i < _numObservations; ++i) {
        fprintf(fptr, "%d %d", _cameraIndex[i], _pointIndex[i]);
        for (int j = 0; j < 2; ++j) {
            fprintf(fptr, " %g", _observations[2 * i + j]);
        }
        fprintf(fptr, "\n");
    }

    for (int i = 0; i < getNumCameras(); ++i) {
        double angleaxis[9];
        if (_useQuaternions) {
            //OutPut in angle-axis format.
            QuaternionToAngleAxis(_parameters + 10 * i, angleaxis);
            memcpy(angleaxis + 3, _parameters + 10 * i + 4, 6 * sizeof(double));
        }
        else {
            memcpy(angleaxis, _parameters + 9 * i, 9 * sizeof(double));
        }
        for (int j = 0; j < 9; ++j) {
            fprintf(fptr, "%.16g\n", angleaxis[j]);
        }
    }

    const double* points = _parameters + getCameraBlockSize() * _numCameras;
    for (int i = 0; i < getNumPoints(); ++i) {
        const double* point = points + i * getPointBlockSize();
        for (int j = 0; j < getPointBlockSize(); ++j) {
            fprintf(fptr, "%.16g\n", point[j]);
        }
    }

    fclose(fptr);
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::writeToPLYFile(const std::string& filename) const {
    std::ofstream of(filename.c_str());
    of << "ply"
        << '\n' << "format ascii 1.0"
        << '\n' << "element vertex " << _numCameras + _numPoints
        << '\n' << "property float x"
        << '\n' << "property float y"
        << '\n' << "property float z"
        << '\n' << "property uchar red"
        << '\n' << "property uchar green"
        << '\n' << "property uchar blue"
        << '\n' << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < getNumCameras(); ++i) {
        const double* camera = getCameras() + getCameraBlockSize() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        of << center[0] << ' ' << center[1] << ' ' << center[2]
            << " 0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double* points = _parameters + getCameraBlockSize() * _numCameras;
    for (int i = 0; i < getNumPoints(); ++i) {
        const double* point = points + i * getPointBlockSize();
        for (int j = 0; j < getPointBlockSize(); ++j) {
            of << point[j] << ' ';
        }
        of << " 255 255 255\n";
    }
    of.close();
}

void BALProblem::CameraToAngelAxisAndCenter(const double* camera, double* angleAxis, double* center) const {
    VectorRef angleAxisRef(angleAxis, 3);
    if (_useQuaternions) {
        QuaternionToAngleAxis(camera, angleAxis);
    }
    else {
        angleAxisRef = ConstVectorRef(camera, 3);
    }

    // c = -R't
    Eigen::VectorXd inverse_rotation = -angleAxisRef;
    AngleAxisRotatePoint(inverse_rotation.data(), camera + getCameraBlockSize() - 6, center);
    VectorRef(center, 3) *= -1.0;
}

void BALProblem::AngleAxisAndCenterToCamera(const double* angleAxis, const double* center, double* camera) const {
    ConstVectorRef angleAxisRef(angleAxis, 3);
    if (_useQuaternions) {
        AngleAxisToQuaternion(angleAxis, camera);
    }
    else {
        VectorRef(camera, 3) = angleAxisRef;
    }

    // t = -R * c
    AngleAxisRotatePoint(angleAxis, center, camera + getCameraBlockSize() - 6);
    VectorRef(camera + getCameraBlockSize() - 6, 3) *= -1.0;
}

void BALProblem::normalize() {
    // Compute the marginal median of the geometry
    std::vector<double> tmp(_numPoints);
    Eigen::Vector3d median;
    double* points = mutablePoints();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < _numPoints; ++j) {
            tmp[j] = points[3 * j + i];
        }
        median(i) = Median(&tmp);
    }

    for (int i = 0; i < _numPoints; ++i) {
        VectorRef point(points + 3 * i, 3);
        tmp[i] = (point - median).lpNorm<1>();
    }

    const double medianAbsoluteDeviation = Median(&tmp);

    // Scale so that the median absolute deviation of the resulting reconstruction is 100

    const double scale = 100.0 / medianAbsoluteDeviation;

    // X = scale * (X - median)
    for (int i = 0; i < _numPoints; ++i) {
        VectorRef point(points + 3 * i, 3);
        point = scale * (point - median);
    }

    double* cameras = mutableCameras();
    double angleAxis[3];
    double center[3];
    for (int i = 0; i < _numCameras; ++i) {
        double* camera = cameras + getCameraBlockSize() * i;
        CameraToAngelAxisAndCenter(camera, angleAxis, center);
        // center = scale * (center - median)
        VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
        AngleAxisAndCenterToCamera(angleAxis, center, camera);
    }
}

void BALProblem::perturb(const double rotationSigma, const double translationSigma, const double pointSigma) {
    assert(pointSigma >= 0.0);
    assert(rotationSigma >= 0.0);
    assert(translationSigma >= 0.0);

    double* points = mutablePoints();
    if (pointSigma > 0) {
        for (int i = 0; i < _numPoints; ++i) {
            PerturbPoint3(pointSigma, points + 3 * i);
        }
    }

    for (int i = 0; i < _numCameras; ++i) {
        double* camera = mutableCameras() + getCameraBlockSize() * i;

        double angleAxis[3];
        double center[3];
        // Perturb in the rotation of the camera in the angle-axis
        // representation
        CameraToAngelAxisAndCenter(camera, angleAxis, center);
        if (rotationSigma > 0.0) {
            PerturbPoint3(rotationSigma, angleAxis);
        }
        AngleAxisAndCenterToCamera(angleAxis, center, camera);

        if (translationSigma > 0.0)
            PerturbPoint3(translationSigma, camera + getCameraBlockSize() - 6);
    }
}