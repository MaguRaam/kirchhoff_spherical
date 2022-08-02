#include "../include/Headers.h"
#include "../include/Exceptions.h"

void ThrowOrderError(unsigned int Order) {
    std::cerr << "Order = " << Order << " not implemented yet." << std::endl; 
    std::exit(EXIT_FAILURE);
}

void ThrowNegativePressureDensityError(Vector<double> WL, Vector<double> WR, Point<2> P, bool boundary) {
    std::cerr << "Negative Pressure/Density" << std::endl; 
    std::cerr << "rhoL = " << WL[0] << ", rhoR = " << WR[0] << std::endl;
    std::cerr << "PL = " << WL[3] << ", PR = " << WR[3] << std::endl; 
    std::cerr << "Location: " << P << std::endl; 

    if (boundary == true) {
        std::cerr << "The face is at boundary"  << std::endl;
    }

    else {
       std::cerr << "The face is in interior"  << std::endl; 
    }

    std::exit(EXIT_FAILURE);
}
 
