//Umit H. Coskun
//Last update : 07/15/2019

#ifndef SHAPE_HPP_INCLUDED_
#define SHAPE_HPP_INCLUDED_

#include "geometry.hpp"
#include <math.h>

class Shape
{
public:
    Shape(Vec3d newPosition, Vec3d newRotationAngles) : position(newPosition), rotationAngles(newRotationAngles)
    {
        calculateRotationMatrix();
        calculateAxisVector();
        shapeIndex = 0;
    }

    Vec3d getPosition()
    {
        return position;
    }

    Vec3d getRotationAngles()
    {
        return rotationAngles;
    }

    void calculateRotationMatrix()
    {
        double psi   = rotationAngles[0];
        double theta = rotationAngles[1];
        double phi   = rotationAngles[2];

        rotationMatrix[0][0] =  cos(theta) * cos(phi);
        rotationMatrix[0][1] =  sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
        rotationMatrix[0][2] =  cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
        rotationMatrix[0][3] =  0.;

        rotationMatrix[1][0] =  cos(theta) * sin(phi); 
        rotationMatrix[1][1] =  sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi);
        rotationMatrix[1][2] =  cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi);
        rotationMatrix[1][3] =  0.;

        rotationMatrix[2][0] = -sin(theta);
        rotationMatrix[2][1] =  sin(psi) * cos(theta);
        rotationMatrix[2][2] =  cos(psi) * cos(theta);
        rotationMatrix[2][3] =  0.;

        rotationMatrix[3][0] =  0.;
        rotationMatrix[3][1] =  0.;
        rotationMatrix[3][2] =  0.;
        rotationMatrix[3][3] =  1.;
    }

    void calculateAxisVector()
    {
        axisVector.x = rotationMatrix[0][2];
        axisVector.y = rotationMatrix[1][2];
        axisVector.z = rotationMatrix[2][2];
    }

    void setPosition(Vec3d newPosition)
    {
        position = newPosition;
    }

    void setRotationAngles(Vec3d newRotationAngles)
    {
        rotationAngles = newRotationAngles;
    }

    Vec3d getAxisVector()
    {
        return axisVector;
    }

    Matrix44d getRotationMatrix()
    {
        return rotationMatrix;
    }

    virtual double getRadius() = 0;
    virtual double getLengthA() = 0;
    virtual double getLengthB() = 0;
    virtual double getHeight() = 0;
    virtual std::vector<Vec3d> getVertices()=0;
    virtual Vec3d getNormalVector(const Vec3d &point)=0;
    unsigned int getShapeIndex()
    {
        return shapeIndex;
    }

protected:
    Vec3d position;
    Vec3d rotationAngles;
    Vec3d axisVector;
    Matrix44d rotationMatrix;
    unsigned int shapeIndex;
};

class Rectangle : public Shape
{
public:

    Rectangle(Vec3d newPosition, Vec3d newEulerAngle, double newLengthA, double newLengthB) : Shape(newPosition, newEulerAngle)
    {
        shapeIndex = 1;
        setLengthA(newLengthA);
        setLengthB(newLengthB);

        // Fix vertex positions before rotation
        vertex.push_back(Vec3d( 0.5 * lengthA, 0.5 * lengthB,0));
        vertex.push_back(Vec3d(-0.5 * lengthA, 0.5 * lengthB,0));
        vertex.push_back(Vec3d(-0.5 * lengthA, -0.5 * lengthB,0));
        vertex.push_back(Vec3d( 0.5 * lengthA, -0.5 * lengthB,0));



        //Rotate
        rotationMatrix.multMatrixDir(vertex.at(0));
        rotationMatrix.multMatrixDir(vertex.at(1));
        rotationMatrix.multMatrixDir(vertex.at(2));
        rotationMatrix.multMatrixDir(vertex.at(3));

        //Translate
        for(int i = 0; i<4; i++)
        {
            vertex.at(i) = vertex.at(i) + position;
        }
    }

    std::vector<Vec3d> getVertices()
    {
        return vertex;
    }
    double getRadius()
    {
        return 0;
    }

    void setLengthA(double newLengthA)
    {
        lengthA = newLengthA;
    }

    void setLengthB(double newLengthB)
    {
        lengthB = newLengthB;
    }

    double getLengthA()
    {
        return lengthA;
    }

    double getLengthB()
    {
        return lengthB;
    }

    double getHeight()
    {
        return 0;
    }

    Vec3d getNormalVector(const Vec3d &point)
    {
        return axisVector;
    }
private:
    double lengthA;
    double lengthB;
    std::vector<Vec3d> vertex;
};

class Disc : public Shape
{
public:

    Disc(Vec3d newPosition, Vec3d newEulerAngle, double newRadius) : Shape(newPosition, newEulerAngle)
    {
        setRadius(newRadius);
        shapeIndex = 2;
    }

    void setRadius(double newRadius)
    {
        radius = newRadius;
    }

    double getRadius()
    {
        return radius;
    }
    
    double getHeight()
    {
        return 0;
    }
    double getLength()
    {
        return 0;
    }
    
    double getLengthA()
    {
        return 0;
    }
    double getLengthB()
    {
        return 0;
    }

    std::vector<Vec3d> getVertices()
    {
        Vec3d ve(0,0,0);
        std::vector<Vec3d> unn;
        unn.push_back(ve);
        return unn;
    }

    Vec3d getNormalVector(const Vec3d &point)
    {
        return axisVector;
    }
private:
    double radius;
};

class Cylinder : public Shape
{
public:

    Cylinder(Vec3d newPosition, Vec3d newEulerAngle, double newRadius, double newHeight) : Shape(newPosition, newEulerAngle)
    {
        shapeIndex = 3;
        setRadius(newRadius);
        setHeight(newHeight);
    }

    void setRadius(double newRadius)
    {
        radius = newRadius;
    }

    double getRadius()
    {
        return radius;
    }

    void setHeight(double newHeight)
    {
        height = newHeight;
    }

    double getHeight()
    {
        return height;
    }

    double getLength()
    {
        return 0;
    }

    double getLengthA()
    {
        return 0;
    }

    double getLengthB()
    {
        return 0;
    }
    
    std::vector<Vec3d> getVertices()
    {
        Vec3d ve(0,0,0);
        std::vector<Vec3d> unn;
        unn.push_back(ve);
        return unn;
    }

    Vec3d getNormalVector(const Vec3d &point)
    {
        return (((point - position).dotProduct(axisVector)) * axisVector - (point - position)).normalize();
    }
private:
    double radius;
    double height;
};

class Sphere : public Shape
{
public:

    Sphere(Vec3d newPosition, Vec3d newEulerAngle, double newRadius) : Shape(newPosition, newEulerAngle)
    {
        setRadius(newRadius);
        shapeIndex = 4;
    }

    void setRadius(double newRadius)
    {
        radius = newRadius;
    }

    double getRadius()
    {
        return radius;
    }

private:
    double radius;
};

#endif // SHAPE_HH_INCLUDED_