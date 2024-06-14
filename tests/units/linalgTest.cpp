#include <iostream>
#include <sstream>
#include <string>
#include <eigen3/Eigen/Dense>

int mxTest();


int main() {
    std::string vec2string(Eigen::Vector3d vector);

    Eigen::Vector3d vec(1.0, 0.0, 0.0);
    vec *= 3;

    double length = vec.norm();
    Eigen::Vector3d normalized = vec / length;

    std::cout << "The length of the vector is: " << length << std::endl;

    Eigen::Vector3d midpoint = (Eigen::Vector3d(1, 0, 0) + Eigen::Vector3d(0, 1, 0)) * 0.5;
    std::cout << "Midpoint: " << vec2string(midpoint) << std::endl;

    mxTest();

    return 0;
}

std::string vec2string(Eigen::Vector3d vector) {
    std::stringstream result;
    result << "[ " << vector.x() << ", "
        << vector.y() << ", "
        << vector.z() << " ]";
    return result.str();
}

int mxTest() {
    Eigen::Matrix4d mat;  // Define a 4x4 matrix of type double

    // Initialize the matrix
    mat << 1, 2, 3, 4,
           5, 6, 7, 8,
           9, 10, 11, 12,
           13, 14, 15, 16;

    // Print the matrix
    std::cout << "4x4 Matrix:\n" << mat << "\n\n";

    return 0;
}

