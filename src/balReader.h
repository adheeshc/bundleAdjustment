#ifndef BAL_READER_H
#define BAL_READER_H

#include <iostream>
#include <string>

class BALProblem {
private:
    int* _pointIndex;
    int* _cameraIndex;
    double* _observations;
    double* _parameters;

    int _numCameras;
    int _numPoints;
    int _numObservations;
    int _numParameters;
    int _useQuaternions;

    void CameraToAngelAxisAndCenter(const double* camera, double* angleAxis, double* center) const;

    void AngleAxisAndCenterToCamera(const double* angleAxis, const double* center, double* camera) const;

public:
    // load bal data from text file
    explicit BALProblem(const std::string& filename, bool useQuaternions = false);

    //destructor
    ~BALProblem() {
        delete[] _pointIndex;
        delete[] _cameraIndex;
        delete[] _observations;
        delete[] _parameters;
    }

    void writeToFile(const std::string& filename) const;
    void writeToPLYFile(const std::string& filename) const;
    void normalize();
    void perturb(const double rotationSigma, const double translationSigma, const double pointSigma);

    int getCameraBlockSize() const { return _useQuaternions ? 10 : 9; };
    int getPointBlockSize() const { return 3; }
    int getNumCameras() const { return _numCameras; }
    int getNumPoints() const { return _numPoints; }
    int getNumObservations() const { return _numObservations; }
    int getNumParameters() const { return _numParameters; }

    const int* getPointIndex() const { return _pointIndex; }
    const int* getCameraIndex() const { return _cameraIndex; }
    const double* getObservations() const { return _observations; }
    const double* getParameters() const { return _parameters; }
    const double* getCameras() const { return _parameters; }
    const double* getPoints() const { return _parameters + getCameraBlockSize() * _numCameras; }

    // address of camera parameter
    double* mutableCameras() { return _parameters; }
    double* mutablePoints() { return _parameters + getCameraBlockSize() * _numCameras; }

    double* mutableCameraForObservation(int i) { mutableCameras() + _cameraIndex[i] * getCameraBlockSize(); }
    double* mutablePointForObservation(int i) { mutablePoints() + _pointIndex[i] * getPointBlockSize(); }
    const double* cameraForObservation(int i) const { getCameras() + _cameraIndex[i] * getCameraBlockSize(); }
    const double* pointForObservation(int i) const { getPoints() + _pointIndex[i] * getPointBlockSize(); }
};

#endif //BAL_READER_H
