//Umit H. Coskun
//Last update : 07/15/2019

#ifndef TRACKER_HPP_INCLUDED_
#define TRACKER_HPP_INCLUDED_

#include "shape.hpp"
#include "geometry.hpp"
#include <gsl/gsl_poly.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <ctime>
#include <boost/numeric/odeint.hpp>

class write_cout
{
private:
    double time;
public:
    write_cout( const double &newTime ) : time(newTime) { }
    void operator() ( const std::vector<double> &x , const double t )
    {
        std::cout << std::setprecision(16) << t + time << " "<< x[0] << " " << x[1] << " " << x[2] << std::endl;
        //std::cout << std::setprecision(16) << t + time << " "<< x[0] << " " << x[1] << " " << x[2] << " " << std::cos((t+time)*1.83247179e2)<< " " << std::sin((t+time)*1.83247179e2)<< std::endl;    
    }
};

class writeTimeSpin
{
private:
    double time;
    std::string path;
    std::ofstream &streamm;
public:
    writeTimeSpin( const double &newTime, std::ofstream &stream,const std::string& newPath) : time(newTime), path(newPath), streamm(stream)
    {}
    void operator() ( const std::vector<double> &x , const double t )
    {
        streamm << std::setprecision(16) << t + time << " "<< x[0] << " " << x[1] << " " << x[2] << std::endl;
        //std::cout << std::setprecision(16) << t + time << " "<< x[0] << " " << x[1] << " " << x[2] << " " << std::cos((t+time)*1.83247179e2)<< " " << std::sin((t+time)*1.83247179e2)<< std::endl;    
    }
};

class EquationSystem
{
private:
    std::vector<double> G;
    std::vector<double> p;
    std::vector<double> v;
    std::vector<double> a;
    std::vector<double> EF;
    double gamma;
    double overCSq = 1.112650056053618432174E-17;

public:
    EquationSystem(const std::vector<double> &newG, const std::vector<double> &newP, const std::vector<double> &newV, const std::vector<double> &newA, const double &newGamma, const std::vector<double> &newEF) 
    : G(newG), p(newP), v(newV), a(newA), gamma(newGamma), EF(newEF)
    {
        
    }

    void operator() ( const std::vector<double> &s , std::vector<double> &dsdt , const double t )
    {
        dsdt[0] = gamma*(0.5*a[1]*G[4]*s[1] + 0.5*a[2]*G[5]*s[1] + 0.5*a[0]*G[6]*s[1] - 0.5*a[0]*G[3]*s[2] - 0.5*a[2]*G[4]*s[2] + 0.25*a[1]*G[5]*s[2] + 0.5*a[1]*G[7]*s[2])*t * t + 
                  gamma*t*(a[1]*EF[0]*overCSq*s[1] - a[0]*EF[1]*overCSq*s[1] + a[2]*EF[0]*overCSq*s[2] - a[0]*EF[2]*overCSq*s[2] + G[6]*s[1]*v[0] - G[3]*s[2]*v[0] + G[4]*s[1]*v[1] + 0.5*G[5]*s[2]*v[1] + G[7]*s[2]*v[1] + G[5]*s[1]*v[2] - G[4]*s[2]*v[2]) + 
                  gamma*(G[1]*s[1] - G[0]*s[2] - EF[1]*overCSq*s[1]*v[0] - EF[2]*overCSq*s[2]*v[0] + EF[0]*overCSq*s[1]*v[1] + EF[0]*overCSq*s[2]*v[2] + G[6]*s[1]*p[0] - G[3]*s[2]*p[0] + G[4]*s[1]*p[1] + 0.5*G[5]*s[2]*p[1] + G[7]*s[2]*p[1] + G[5]*s[1]*p[2] - G[4]*s[2]*p[2]);


        dsdt[1] = gamma*(-0.5*a[1]*G[4]*s[0] - 0.5*a[2]*G[5]*s[0] - 0.5*a[0]*G[6]*s[0] + 0.5*a[1]*G[3]*s[2] - 0.25*a[0]*G[5]*s[2] + 0.5*a[2]*G[6]*s[2] + 0.5*a[0]*G[7]*s[2])*t * t + 
                  gamma*t*(-(a[1]*EF[0]*overCSq*s[0]) + a[0]*EF[1]*overCSq*s[0] + a[2]*EF[1]*overCSq*s[2] - a[1]*EF[2]*overCSq*s[2] - G[6]*s[0]*v[0] - 0.5*G[5]*s[2]*v[0] + G[7]*s[2]*v[0] - G[4]*s[0]*v[1] + G[3]*s[2]*v[1] - G[5]*s[0]*v[2] + G[6]*s[2]*v[2]) + 
                  gamma*(-(G[1]*s[0]) + G[2]*s[2] + EF[1]*overCSq*s[0]*v[0] - EF[0]*overCSq*s[0]*v[1] - EF[2]*overCSq*s[2]*v[1] + EF[1]*overCSq*s[2]*v[2] - G[6]*s[0]*p[0] - 0.5*G[5]*s[2]*p[0] + G[7]*s[2]*p[0] - G[4]*s[0]*p[1] + G[3]*s[2]*p[1] - G[5]*s[0]*p[2] + G[6]*s[2]*p[2]);


        dsdt[2] = gamma*(0.5*a[0]*G[3]*s[0] + 0.5*a[2]*G[4]*s[0] - 0.25*a[1]*G[5]*s[0] - 0.5*a[1]*G[7]*s[0] - 0.5*a[1]*G[3]*s[1] + 0.25*a[0]*G[5]*s[1] - 0.5*a[2]*G[6]*s[1] - 0.5*a[0]*G[7]*s[1])*t * t + 
                  gamma*t*(-(a[2]*EF[0]*overCSq*s[0]) + a[0]*EF[2]*overCSq*s[0] - a[2]*EF[1]*overCSq*s[1] + a[1]*EF[2]*overCSq*s[1] + G[3]*s[0]*v[0] + 0.5*G[5]*s[1]*v[0] - G[7]*s[1]*v[0] - 0.5*G[5]*s[0]*v[1] - G[7]*s[0]*v[1] - G[3]*s[1]*v[1] + G[4]*s[0]*v[2] - G[6]*s[1]*v[2]) + 
                  gamma*(G[0]*s[0] - G[2]*s[1] + EF[2]*overCSq*s[0]*v[0] + EF[2]*overCSq*s[1]*v[1] - EF[0]*overCSq*s[0]*v[2] - EF[1]*overCSq*s[1]*v[2] + G[3]*s[0]*p[0] + 0.5*G[5]*s[1]*p[0] - G[7]*s[1]*p[0] - 0.5*G[5]*s[0]*p[1] - G[7]*s[0]*p[1] - G[3]*s[1]*p[1] + G[4]*s[0]*p[2] - 
                        G[6]*s[1]*p[2]);
                
    }
};


class PathTracer
{
private:
    const std::vector<Shape*> shape;

    Vec3d preVelocity;
    Vec3d velocity;
    Vec3d position;
    double speed;
    Vec3d prePosition;
    Vec3d gravity;
    Vec3d axisVector;
    Vec3d spin;

    std::vector<double> gList; // G0-1, G00, G01, G1-2, G1-1, G10, G11, G12
    std::vector<double> EF;    // Electric field

    double initTime;
    double time;
    double finalTime;
    double gamma;

    double sumZ;
    int count;

    const gsl_rng_type * rngType;
    gsl_rng * myRng;

    std::ofstream timeSpinOut;

    //flags
    bool isLastRun;

    //some constants
    double kBoltz = 1.38064852e-23; // m^2 kg s^-2 K^-1
    double mEff   = 1.1018115853052553036E-26; // kg
    double diffusion = 1.6e-4; // m^2 s^-1 K^-7
    
    // Most prob. He3 speed
    double vMP;

    double tauP;

    // Particle type 1 = UCN 2 = 3He
    int particleType;

    //Temperature
    double temp; // In K

    //Speed of the UCN
    double ucnSpeed;

    //Maximum speed of He3
    double maxSpeed;

    //Output path names
    std::string timeSpinPath;

public:
    PathTracer(const std::vector<Shape*> newShapeList, Vec3d &newSpin, Vec3d &newVelocity, Vec3d &newPosition, double newGravity, double &newInitTime, double newFinalTime, double newGamma)
    : shape(newShapeList), spin(newSpin), velocity(newVelocity), position(newPosition), gravity({0.,0.,newGravity}), initTime(newInitTime), finalTime(newFinalTime), gamma(newGamma)
    {
        isLastRun = false;

        gsl_rng_env_setup();
        rngType = gsl_rng_default;
        myRng = gsl_rng_alloc (rngType);
    }

    void setSeed(const long unsigned int &seed)
    {
        gsl_rng_set(myRng, seed);
    }

    Vec3d getPosition()
    {
        return position;
    }

    Vec3d getVelocity()
    {
        return velocity;
    }

    double solveQuadratic(double coefficients[3])
    {
        double epsilon = 1e-09;
        double a = coefficients[0];
        double b = coefficients[1];
        double c = coefficients[2];

        double firstSolution;
        double secondSolution;


        double solution;
        gsl_poly_solve_quadratic(a, b, c, &firstSolution, &secondSolution);

        if ( std::fabs( firstSolution ) < epsilon || firstSolution <- epsilon )
        {
            solution = secondSolution;
        }
        else
        {
            solution = firstSolution;
        }
        //std::cout << "Solution = " << solution << std::endl;
        return solution;
    }

    double solveCubic ( double coefficients[3] )
    {
        double epsilon = 1e-09;
        double a = coefficients[0];
        double b = coefficients[1];
        double c = coefficients[2];

        double firstSolution;
        double secondSolution;
        double thirdSolution;

        double solution;
        double rootNumber = gsl_poly_solve_cubic(a, b, c, &firstSolution, &secondSolution, &thirdSolution);

        if (rootNumber == 1)
        {
            solution = firstSolution;
        }
        else
        {
            if (firstSolution > epsilon)
            {
                solution = firstSolution;
            }
            else if(secondSolution > epsilon)
            {
                solution = secondSolution;
            }
            else
            {
                solution = thirdSolution;
            }            
        }
        return solution;
    }

    double solveQuartic ( std::vector<double> coefficientsVec )
    {
        double *coefficients = &coefficientsVec[0];
        double epsilon = 1e-13;
        
        double solution;

        if ( std::fabs( coefficients[0] ) < epsilon && std::fabs( coefficients[1] ) < epsilon )
        {   
            double newCoefficients[3] = { coefficients[2], coefficients[3], coefficients[4] };
            solution = solveQuadratic ( newCoefficients );
        }
        else if( std::fabs( coefficients[0] ) < epsilon && std::fabs( coefficients[1] ) > epsilon )
        {
            double newCoefficients[3] = {coefficients[2]/coefficients[1], coefficients[3]/coefficients[1], coefficients[4]/coefficients[1]};
            solution = solveCubic(newCoefficients);
        }
        else
        {
            int degree = 4;
            double z[degree*2];

            gsl_poly_complex_workspace * w
            = gsl_poly_complex_workspace_alloc (degree + 1);

            gsl_poly_complex_solve (coefficients, degree + 1, w, z);

            gsl_poly_complex_workspace_free (w);

            std::vector<double> tempRoots;

            for(int i = 0; i<degree; i++)
            {
                if( std::fabs(z[2*i+1] < epsilon) && z[2*i] > epsilon )
                {
                    tempRoots.push_back(z[2*i]);
                }
            }
            solution = (*std::min_element(tempRoots.begin(), tempRoots.end()));
        }

        return solution;
    }

    void cylinderIntersection(Shape *currentShape, std::vector<double> &timeList, std::vector<Vec3d> &normalList, std::vector<Vec3d> &positionList)
    {
        double epsilon = 1e-10;

        Vec3d  b       = currentShape->getPosition();
        Vec3d  n       = currentShape->getAxisVector();
        double radius  = currentShape->getRadius();
        double height  = currentShape->getHeight();

        Vec3d p0B = position - b;

        Vec3d p0BCrossN = p0B.crossProduct(n);
        Vec3d v0CrossN  = velocity.crossProduct(n);
        Vec3d aCrossN   = 0.5 * gravity.crossProduct(n);

        std::vector<double> coefficients;

        coefficients.push_back(aCrossN.dotProduct(aCrossN));
        coefficients.push_back(2. * v0CrossN.dotProduct(aCrossN));
        coefficients.push_back(v0CrossN.dotProduct(v0CrossN) + 2. * p0BCrossN.dotProduct(aCrossN));
        coefficients.push_back(2. * p0BCrossN.dotProduct(v0CrossN));
        coefficients.push_back(p0BCrossN.dotProduct(p0BCrossN) - radius * radius);

        double collisionTime = solveQuartic(coefficients);

        Vec3d temPosition = position + velocity * collisionTime + 0.5 * gravity * collisionTime * collisionTime;

        // Inside-outside test
        double distance = ((temPosition-b).crossProduct(n)).length();
        double distanceFromCenter = (temPosition-b).dotProduct(n);

        Vec3d normalVector = currentShape->getNormalVector(temPosition);
        if(distanceFromCenter < height && distanceFromCenter > epsilon && distance <= radius+epsilon && collisionTime > epsilon)
        {
            timeList.push_back(collisionTime);
            normalList.push_back(normalVector);
            positionList.push_back(b);
        }
    
    }

    void discIntersection(Shape *currentShape, std::vector<double> &timeList, std::vector<Vec3d> &normalList, std::vector<Vec3d> &positionList)
    {
        //For numbers which are really zero
        double epsilon = 1e-12;

        // Position and axis vector(normal vector for the square) of the shape
        Vec3d b  = currentShape->getPosition();
        Vec3d n  = currentShape->getAxisVector();
        double r = currentShape->getRadius();

        double coefficients[3];

        coefficients[0] = 0.5 * gravity.dotProduct(n);
        coefficients[1] = velocity.dotProduct(n);
        coefficients[2] = (position-b).dotProduct(n);

        double collisionTime = solveQuadratic(coefficients);

        Vec3d temPosition = position + velocity * collisionTime + 0.5 * gravity * collisionTime * collisionTime;

        // Inside-outside test
        double distance = std::fabs((temPosition-b).dotProduct(n));
        double distanceFromCenter = (temPosition-b).length();
        Vec3d normalVector = currentShape->getNormalVector(temPosition);
        if(distanceFromCenter < r && distance < epsilon && collisionTime > epsilon)
        {
            timeList.push_back(collisionTime);
            normalList.push_back(normalVector);
            positionList.push_back(b);
        }
    }

    void rectangleIntersection(Shape *currentShape, std::vector<double> &timeList, std::vector<Vec3d> &normalList, std::vector<Vec3d> &positionList)
    {
        //For numbers which are really zero
        double epsilon = 1e-14;

        // Position and axis vector(normal vector for the square) of the shape
        Vec3d b = currentShape->getPosition();
        Vec3d n = currentShape->getAxisVector();

        double coefficients[3];

        coefficients[0] = 0.5 * gravity.dotProduct(n);
        coefficients[1] = velocity.dotProduct(n);
        coefficients[2] = (position-b).dotProduct(n);

        double collisionTime = solveQuadratic(coefficients);

        Vec3d temPosition = position + velocity * collisionTime + 0.5 * gravity * collisionTime * collisionTime;

        // Inside-outside test
        std::vector<Vec3d> vertices = currentShape->getVertices();
        Vec3d AB = vertices.at(0) - vertices.at(1);
        Vec3d AM = vertices.at(0) - temPosition;
        Vec3d AD = vertices.at(0) - vertices.at(3);

        double distance = std::fabs((temPosition-b).dotProduct(n));

        Vec3d normalVector = currentShape->getNormalVector(temPosition);
        if(AM.dotProduct(AB) > 0 && AB.dotProduct(AB) > AM.dotProduct(AB) &&  AM.dotProduct(AD) > 0 && AD.dotProduct(AD) >AM.dotProduct(AD) && distance < epsilon && collisionTime > epsilon)
        {
            timeList.push_back(collisionTime);
            normalList.push_back(normalVector);
            positionList.push_back(b);
        }
    
    }

    void setBField(const std::vector<double> &bField)
    {
        gList = bField;
    }

    void setEField(const std::vector<double> &eField)
    {
        EF = eField;
    }

    void setRandomPosition()
    {
        Vec3d temp(gsl_rng_uniform (myRng) * 0.102 - 0.051,gsl_rng_uniform (myRng) * 0.4 - 0.2 , 0/*gsl_rng_uniform (myRng) * 0.075 - 0.0375*/);
        position = temp;
    }

    void setParticleType(int type)
    {
        particleType = type;
    }

    void setSpeed(int type)
    {
        if(type==1)
        {
            speed = ucnSpeed;
        }
        else if(type==2)
        {
            bool accept = false;

            while(!accept)
            {
                double testRand = gsl_rng_uniform (myRng) * He3Distrubition(vMP);
                double vCandidate = gsl_rng_uniform (myRng) * maxSpeed;
                if(He3Distrubition(vCandidate) > testRand)
                {
                    accept = true;
                }
                speed = vCandidate;
            }
        }

    }

    void setUcnSpeed(const double& newUcnSpeed)
    {
        ucnSpeed = newUcnSpeed;
    }

    void setMaxSpeed(const double& newMaxSpeed)
    {
        maxSpeed = newMaxSpeed;
    }

    void setRandomVelocity()
    {
        Vec3d temp(gsl_rng_uniform (myRng),gsl_rng_uniform (myRng),0.);
        temp.normalize();
        temp = speed * temp;
        velocity = temp;
    }

    void run(int runNumber)
    {
        bool phononCollision = false;
        auto spinResult = spin.toStdVec();
        time = initTime;

        // Set random position and velocity
        setSpeed(particleType);
        std::cout << speed << std::endl;
        setRandomPosition();
        setRandomVelocity();

        timeSpinOut << std::setprecision(16) << velocity;
        while(isLastRun == false)
        {
            double timeEpsilon = 1e-14;
            std::vector<double> timeList;
            std::vector<Vec3d>  normalList;
            std::vector<Vec3d>  positionList;

            for(int i=0; i<shape.size();i++)
            {
                auto currentShape = shape.at(i);
                auto currentShapeIndex = currentShape->getShapeIndex();

                if(currentShapeIndex==1)
                {
                    rectangleIntersection(currentShape, timeList, normalList, positionList);
                }
                if(currentShapeIndex==2)
                {
                    discIntersection(currentShape, timeList, normalList, positionList);
                }
                if(currentShapeIndex==3)
                {
                    cylinderIntersection(currentShape, timeList, normalList, positionList);
                }
            }

            double collisionTime = (*std::min_element(timeList.begin(), timeList.end()));
            int    minPosition   = std::distance(timeList.begin(), std::min_element(timeList.begin(), timeList.end()));

            Vec3d currentNormalVector  = normalList.at(minPosition);
            Vec3d currentPositonVector = positionList.at(minPosition);

            double beginTime    = time;
            Vec3d beginPosition = position;
            Vec3d beginVelocity = velocity;

            auto vPosition = position.toStdVec();
            auto vVelocity = velocity.toStdVec();
            auto vGravity  = gravity.toStdVec();

            if(particleType == 2)
            {
                double He3colTime = generatePhononHe3CollisionTime();
                
                if(He3colTime < collisionTime && collisionTime > 1e6)
                {
                    collisionTime = He3colTime;
                    phononCollision = true;
                }
            }

            if(time + collisionTime > finalTime)
            {
                collisionTime = finalTime - time;
                isLastRun = true;
            }
            
            EquationSystem eqs(gList, vPosition, vVelocity, vGravity, gamma, EF);

            writeTimeSpin woo(time, timeSpinOut, timeSpinPath);

            boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled( 1E-16 , 1E-16 , boost::numeric::odeint::runge_kutta_fehlberg78< std::vector<double> >() ) ,
                        eqs , spinResult , 0.0 , collisionTime , 1E-15);

            time += collisionTime;

            position = position + velocity * collisionTime + 0.5 * gravity * collisionTime * collisionTime;
            velocity = velocity + collisionTime * gravity;
            
            //std::cout << position <<std::endl;
            if(phononCollision)
            {
                setRandomVelocity();
                phononCollision = false;
            }
            else
            {
                // Specular reflection
                velocity = velocity - 2. * (velocity.dotProduct(currentNormalVector)) * currentNormalVector;
            }
            
            // Set any number close to 0 to 0.
            numericalCorrect(position);
            numericalCorrect(velocity);

        }
        isLastRun = false;
        timeSpinOut << std::setprecision(16) << " "<< spinResult[0] << " " << spinResult[1] << " " << spinResult[2] << std::endl;
    }

/*
    void printPath(double collisionTime, Vec3d beginPosition, Vec3d beginVelocity, double beginTime)
    {
        double dt = 0.0001;
        double t = 0;
        while(t<=collisionTime)
        {                
            trajectoryOut << std::setprecision(6) << beginTime + t <<" " <<beginPosition + t * beginVelocity + 0.5 * t * t * gravity << std::endl;
            t += dt;
        }
    }
*/
    void numericalCorrect(Vec3d &vector)
    {
        double epsilon = 1e-12;
        for(int i=0; i<3; i++)
        {
            if(std::fabs(vector[i]) < epsilon)
            {
                vector[i] = 0.;
            }
        }
    }

    // He3 Specific Functions
    void setTemperature(const double& newTemp)
    {
        temp = newTemp;
        vMP = std::sqrt(2. * kBoltz * temp / mEff);
        tauP = 2. * diffusion / (vMP * vMP); 
    }

    double He3Distrubition(const double &v)
    {
        // note that 1/sqrt(pi) = 0.5641895835477562869481
        // mEff is the effective mass of 3He diffusing in liquid Helium (In Riccardo Schmidt's thesis the notation is M*)
        return 2. / (vMP * vMP * vMP * vMP) * v * v * v * std::exp(-v * v / (vMP * vMP));
    }

    double phononHe3Distrubition(const double &t)
    {
        return std::exp(- t / tauP) / tauP;
    }

    double generatePhononHe3CollisionTime()
    {
        bool accept = false;
        double tCol;
        while(!accept)
        {
            double testRand = gsl_rng_uniform (myRng) / tauP;
            double tCandidate = gsl_rng_uniform (myRng) * 0.001;
            if(phononHe3Distrubition(tCandidate) > testRand)
            {
                accept = true;
            }
            tCol = tCandidate;
        }
        return tCol;
    }

    void setTimeSpinPath(const std::string& pathname)
    {
        timeSpinPath = pathname;
        if(timeSpinOut.is_open())
        {
            timeSpinOut.close();
            timeSpinOut.open(timeSpinPath + ".out");
        }
        else
        {
            timeSpinOut.open(timeSpinPath + ".out");
        }        
    }
};

#endif // TRACKER_HPP_INCLUDED_
