//Umit H. Coskun
//Last update : 07/15/2019
//To read the parameter files and initialize the prgram

#ifndef PARAMREAD_HPP_INCLUDED_
#define PARAMREAD_HPP_INCLUDED_

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "geometry.hpp"

class Parameters
{
public:
    double gamma;
    double initialTime;
    double finalTime;
    double gravityZ;
    unsigned long int seed;
    unsigned int numberOfParticles;
    double maxSpeed;
    int particleType;
    double temperature;
    double ucnSpeed;
    std::vector<double> BField{0,0,0,0,0,0,0,0};
    std::vector<double> EField{0,0,0};
    Vec3d spin{0,0,0};
    std::string spinPath;

    Parameters()
    {
        // Initializes parameters and initial conditions etc
        std::ifstream pfin;

        pfin.open("config/parameters.dat");  
        if (pfin.is_open())
        { 
            pfin >> gamma;
            pfin >> initialTime;
            pfin >> finalTime;
            pfin >> gravityZ;
            pfin >> seed;
            pfin >> numberOfParticles;
            pfin >> maxSpeed;
            pfin >> particleType;
            pfin >> temperature;
            pfin >> ucnSpeed;
        }
        else
        {
            std::cout << "parameters.dat cannot be read!" << std::endl;
            exit(1);
        }
        pfin.close();

        pfin.open("config/bfield.dat");  
        if (pfin.is_open())
        { 
            pfin >> BField.at(0);
            pfin >> BField.at(1);
            pfin >> BField.at(2);
            pfin >> BField.at(3);
            pfin >> BField.at(4);
            pfin >> BField.at(5);
            pfin >> BField.at(6);
            pfin >> BField.at(7);
        }
        else
        {
            std::cout << "efield.dat cannot be read!" << std::endl;
            exit(1);
        }
        pfin.close();

        pfin.open("config/efield.dat");  
        if (pfin.is_open())
        { 
            pfin >> EField.at(0);
            pfin >> EField.at(1);
            pfin >> EField.at(2);
        }
        else
        {
            std::cout << "efield.dat cannot be read!" << std::endl;
            exit(1);
        }
        pfin.close();

        pfin.open("config/spin.dat");  
        if (pfin.is_open())
        { 
            pfin >> spin[0];
            pfin >> spin[1];
            pfin >> spin[2];
        }
        else
        {
            std::cout << "spin.dat cannot be read!" << std::endl;
            exit(1);
        }
        pfin.close();

        pfin.open("config/outconfig.dat");  
        if (pfin.is_open())
        { 
            pfin >> spinPath;
        }
        else
        {
            std::cout << "outconfig.dat cannot be read!" << std::endl;
            exit(1);
        }
        pfin.close();
    }

    void printParametes()
    {
        std::cout << "||from config/parameters.dat :" << std::endl;
        std::cout << "||gamma        = " << gamma << std::endl;
        std::cout << "||Initial Time = " << initialTime << std::endl;
        std::cout << "||Final Time   = " << finalTime << std::endl;
        std::cout << "||gravityZ     = " << gravityZ << std::endl;
        std::cout << "||SEED         = " << seed << std::endl;
        std::cout << "||NoP          = " << numberOfParticles << std::endl;
        std::cout << "||vMAX of He3s = " << maxSpeed<< " m/s" << std::endl;
        std::cout << "||Part. Type   = " << particleType<< " (1 = UCN, 2 = He3)" << std::endl;
        std::cout << "||Temperature  = " << temperature << " K" << std::endl;
        std::cout << "||vUCN         = " << ucnSpeed << " m/s" << std::endl;


        std::cout << "==========================================================" << std::endl;

        std::cout << "||from config/bfield.dat" << std::endl;
        std::cout << "||"<<BField.at(0) << " " << BField.at(1) << " " << BField.at(2) << " " << BField.at(3) << " " << BField.at(4) << " " << BField.at(5) << " " << BField.at(6) << " " << BField.at(7) << std::endl;

        std::cout << "==========================================================" << std::endl;

        std::cout << "||from config/efield.dat" << std::endl;
        std::cout << "||"<<EField.at(0) << " " << EField.at(1) << " " << EField.at(2) << std::endl;

        std::cout << "==========================================================" << std::endl;

        std::cout << "||from config/spin.dat"  << std::endl;
        std::cout << "||"<<spin << std::endl;

        std::cout << "==========================================================" << std::endl;
    }

};

#endif
