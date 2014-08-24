#ifndef SPHSolver_H
#define SPHSolver_H

#include <set>
#include <vector>
#include <iostream>
#include <math.h>

#include "Particle.h"
#include "PDCFileWriter.h"

/// \file SPHSolver.h
/// \brief Contains the SPHSolver class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class SPHSolver
/// \brief The implementation of SPH method. Finds new position and velocity values for each particle

class SPHSolver
{

public:
	/// \brief Ctor
    /// @param[in] _pdcWriter the reference to the PDCFileWriter object
    /// @param[in] _particle the reference to the particle vector to be simulated
    /// @param[in] _startFrame the start frame of the simulation
    /// @param[in] _frameCount the number of frames to be cached out	
	SPHSolver( PDCFileWriter& pdcWriter, std::vector<Particle*>& particle, int startFrame, int frameCount );	

	//~SPHSolver();

	/// \brief Smoothing Kernel used to calculate density for each particle
    /// @param[in] _magnitudeSquared the disctance between the particle and its neighbour
	double WKernelPoly( double magnitudeSquared );

	/// \brief Smoothing Kernel used to calculate surface tension force and interface tension force
    /// @param[in] _ magnitudeVectorCoor the magnitude vector value (could be X,Y or Z value)
	/// @param[in] _magnitudeSquared the disctance between the particle and its neighbour
	double WKernelPolyGradient( double magnitudeVectorCoor, double magnitudeSquared );

	/// \brief Smoothing Kernel used to calculate surface tension force and interface tension force
    /// @param[in] _magnitudeSquared the disctance between the particle and its neighbour
	double WKernelPolyLaplacian( double magnitudeSquared );

	/// \brief Smoothing Kernel used to calculate Pressure Force
	/// @param[in] _ magnitudeVectorCoor the magnitude vector value (could be X,Y or Z value)
	/// @param[in] _magnitude the disctance between the particle and its neighbour
	double WKernelPressureGradient( double magnitudeVectorCoor, double magnitude );

	/// \brief Smoothing Kernel used to calculate Viscosity (Friction) Force
   /// @param[in] _magnitude the disctance between the particle and its neighbour
	double WKernelViscosityLaplacian( double magnitude );

	/// \brief Starts simulation
	void runSimulation();

	/// \brief Updates density and pressure values for each particle
	/// @param[in] _particle the particle to be updated
    /// @param[in] _neighbourParticle the reference to the particle neighbour vector    
	void updateDensityPressure( Particle* particle, std::vector<Particle*>& neighbourParticle );

	/// \brief Updates position and velocity for each particle
	/// @param[in] _particle the particle to be updated
    /// @param[in] _neighbourParticle the reference to the particle neighbour vector
	void updatePositionVelocity( Particle * particle, std::vector<Particle*>& neighbourParticle);
	
	
private:

	/// \brief The vector of the particles to be simulated
	std::vector<Particle*>& m_particle;

	/// \brief PDCFileWriter object
	PDCFileWriter& m_PDCFileWriter;

	/// \brief The radius for the neighbour search
	double m_neighbourRadius;

	//the values below can be precalculated in advance

	/// \brief The radius for the neighbour search (squared)
	double m_radiusSquared;

	/// \brief Used to calculate WKernelPoly
	double m_weightPoly;
	
	/// \brief Used to calculate WKernelPolyGradient
	double m_weightPolyGradient;

	/// \brief Used to calculate WKernelPolyLaplacian
    double m_weightPolyLaplacian;

	/// \brief Used to calculate WKernelPressureGradient
    double m_weightPressureGradient;

	/// \brief Used to calculate WKernelViscosityLaplacian
    double m_weightViscosityLaplacian;	

};

#endif
