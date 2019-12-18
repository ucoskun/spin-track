//Umit H. Coskun
//Last update : 07/15/2019
//To compile : g++ -lgsl -lgslcblas -lm -O3 spintrack.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <memory>
#include "shape.hpp"
#include "geometry.hpp"
#include "tracker.hpp"
#include "paramread.hpp"
#include <ctime>

int main(void) 
{
    //Create the rectangular box
    
    double dim1 = 0.102;
    double dim2 = 0.4;
    double dim3 = 0.075;

    Vec3d pos0(0.,0.,-dim3/2.);
    Vec3d dir0(0.,0.,0.);

    Vec3d pos1(0.,0.,dim3/2);
    Vec3d dir1(-M_PI,0.,0.);

    Vec3d pos2(0.,-dim2/2.,0.);
    Vec3d dir2(-M_PI/2.0,0.,0.);

    Vec3d pos3(0.,dim2/2.,0.);
    Vec3d dir3(M_PI/2.0,0.,0.);

    Vec3d pos4(dim1/2.,0.,0.);
    Vec3d dir4(0.,-M_PI/2.0,0.);

    Vec3d pos5(-dim1/2.,0.,0.);
    Vec3d dir5(0.,M_PI/2.0,0.);


    std::vector<Shape *> myObject;
    myObject.push_back(new Rectangle(pos0, dir0, dim1, dim2));
    myObject.push_back(new Rectangle(pos1, dir1, dim1, dim2));
    myObject.push_back(new Rectangle(pos2, dir2, dim1, dim3));
    myObject.push_back(new Rectangle(pos3, dir3, dim1, dim3));
    myObject.push_back(new Rectangle(pos4, dir4, dim3, dim2));
    myObject.push_back(new Rectangle(pos5, dir5, dim3, dim2));

    // Position and velocity vectors will be randomly assigned bt PathTracer
    Vec3d position;
    Vec3d velocity;

    //Initialize parameters
    Parameters par;
    par.printParametes();

    // Create a PathTracer to track the motion and spin
    PathTracer myTracker( myObject, par.spin, velocity, position, par.gravityZ,
                         par.initialTime, par.finalTime, par.gamma );
    
    myTracker.setParticleType(par.particleType);
    myTracker.setSeed(par.seed);

    // Set configurations
    myTracker.setBField(par.BField);
    myTracker.setEField(par.EField);

    // Set UCN speed
    myTracker.setUcnSpeed(par.ucnSpeed);

    // Set max speed of He3
    myTracker.setTemperature(par.temperature);
    myTracker.setMaxSpeed(par.maxSpeed);

    //Set output name
    myTracker.setTimeSpinPath(par.spinPath);

    // Run pathtracing numberOfParticles times 
    for ( int i = 0 ; i < par.numberOfParticles; i++ )
    {
        std::clock_t begin = clock();
        myTracker.run(i);
        std::clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        std::cout << "Computation time = " << elapsed_secs << std::endl;
    }
    return 0;
}


    /*
    //Create a cube
    double len = 0.5;
    Vec3d pos0(0.,0.,0.);
    Vec3d dir0(0.,0.,0.);

    Vec3d pos1(0.,0.,len);
    Vec3d dir1(-M_PI,0.,0.);

    Vec3d pos2(0.,-len/2.,len/2.);
    Vec3d dir2(-M_PI/2.0,0.,0.);

    Vec3d pos3(0.,len/2.,len/2.);
    Vec3d dir3(M_PI/2.0,0.,0.);

    Vec3d pos4(len/2.,0.,len/2.);
    Vec3d dir4(0.,-M_PI/2.0,0.);
    
    Vec3d pos5(-len/2.,0.,len/2.);
    Vec3d dir5(0.,M_PI/2.0,0.);

    std::vector<Shape *> myObject;
    myObject.push_back(new Rectangle(pos0, dir0, len, len));
    myObject.push_back(new Rectangle(pos1, dir1, len, len));
    myObject.push_back(new Rectangle(pos2, dir2, len, len));
    myObject.push_back(new Rectangle(pos3, dir3, len, len));
    myObject.push_back(new Rectangle(pos4, dir4, len, len));
    myObject.push_back(new Rectangle(pos5, dir5, len, len));
    */

    /*
    //Create a cylinder (ILL Chamber)
    double radius = 0.235;
    double height = 0.12;
    Vec3d posSide( 0., 0., -height*0.5 );
    Vec3d dirSide( 0., 0., 0. );

    Vec3d posTop( 0., 0., height*0.5 );
    Vec3d dirTop( -M_PI, 0., 0. );

    Vec3d posBottom( 0., 0., -height*0.5 );
    Vec3d dirBottom( 0., 0., 0. );

    // Initialize the object and push back the components
    std::vector< Shape* > myObject;
    myObject.push_back( new Cylinder ( posSide, dirSide, radius, height ) );
    myObject.push_back( new Disc ( posTop, dirTop, radius ) );
    myObject.push_back( new Disc ( posBottom, dirBottom, radius ) );
    */
