#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

#include "Cell.h"
#include "CellManager.h"

/// \file Particle.h
/// \brief Contains the Particle class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class Particle
/// \brief Creates a particle class. Runs the Neighbour search algorithm.

class Particle
{
public:

	/// \brief Default ctor
	Particle();
	//~Particle();
	
	/// \brief Sets the id of the particle
	void setId( void* input );		
		
	/// \brief Sets the coordinates of the particle position
	void setPosition( void* input );		
		
	/// \brief Sets the values of the particle velocity
	void setVelocity( void* input );		
		
	/// \brief Sets the mass of the particle
	void setMass( void* input );

	/// \brief Sets the rest density of the substance
	void setRestDensity( void* input );

	/// \brief Sets the viscosity constant of the substance
	void setViscosityConstant( void* input );

	/// \brief Sets the radius used for the neighbour search
	void setNeighbourRadius( void* input ); 

	/// \brief Sets the gas stiffness constant
	void setGasConstant( void* input );
	
	/// \brief Sets the surface tension value of the substance
	void setSurfaceTensionConstant( void* input );

	/// \brief setter typedefs
	typedef void (Particle::*Setter)( void* );
		
	/// \brief According to a string input calls the right method
	/// @param[in] _name the string to allocate the method to be called
	static Setter getSetterMethod( std::string _name );	
	
	/// \brief The pointer to the cell manager object created in the main method
	static CellManager* m_cellManager;
	
	/// \brief Updates particle position and velocity
	/// @param[in] _pX the x coordinate of the particle position
	/// @param[in] _pY the y coordinate of the particle position
	/// @param[in] _pZ the z coordinate of the particle position
	/// @param[in] _vX the x value of the particle velocity vector
	/// @param[in] _vY the y value of the particle velocity vector
	/// @param[in] _vZ the z value of the particle velosity vector	
	void setNewPositionVelocity( double _pX , double _pY , double _pZ, double _vX, double _vY, double _vZ );	

	/// \brief Sets new value of density
	/// @param[in] _density the new value of density
	void setDensity( double _density );

	/// \brief Calculates pressure
	void calculatePressure();
			
	/// \brief Sets cell for the particle in which the particle is
	/// @param[in] _cell the pointer to the cell
	void setCell( Cell* _cell ); 

	/// \brief Gets particle neighbours for the particle
	/// @param[in] _particleNeighbour the reference to the particle neighbour vector to be filled with particle neighbours
	void getNeighbours( std::vector<Particle*>& _particleNeighbour );
	
	/// \brief Gets the particle id
	unsigned getId(); 

	/// \brief Gets the particle position coordinate
	/// @param[in] _i the index of the particle position coordinate
	double getPosition( unsigned _i) const;	

	/// \brief Gets the particle velocity value
	/// @param[in] _i the index of the particle velocity value
	double getVelocity( unsigned _i ) const;	

	/// \brief Gets the particle mass
	double getMass() const;
	
	/// \brief Gets the particle density
	double getDensity() const;
	
	/// \brief Gets pressure acting on the particle
	double getPressure() const;
	
	/// \brief Gets the rest density of the substance
	double getSubstanceRestDensity() const;
	
	/// \brief Gets the viscosity constant of the substance
	double getViscosityConstant() const;

	/// \brief Gets the radius used for the neighbour search
	double getNeighbourRadius() const;	
		
	// surface tension
	/// \brief Gets the surface color coefficient
	double getSurfaceColorCoefficient() const;
	
	/// \brief Gets the surface tension threshold
	double getSurfaceTensionThreshold() const;
	
	/// \brief Gets the surface tension coefficient
	double getSurfaceTensionCoefficient() const;
	
	// interface tension
	/// \brief Gets the interface color coefficient
	double getInterfaceColorCoefficient() const;
	
	/// \brief Gets the interface tension threshold
	double getInterfaceTensionThreshold() const;
	
	/// \brief Gets the interface tension coefficient
	double getInterfaceTensionCoefficient() const;
	
private:

	/// \brief Holds the particle id
	unsigned 	m_id;				// get from .pda file

	/// \brief Holds the coordinates of the particle position
	double 		m_position[3];		// get from .pda file and later calculate in the program

	/// \brief Holds the values of the particle velocity
	double 		m_velocity[3];		// get from .pda file and later calculate in the program

	/// \brief Holds the particle mass
	double		m_mass;				// get from .pda file

	/// \brief Holds the particle density
	double		m_density;			// calculate in the program

	/// \brief Holds the value of the pressure acting on the particle
	double		m_pressure;			// calculate in the program

	/// \brief Holds the value of the radius used for the neighbour search
	double		m_neighbourRadius;	// get from .pda file

	/// \brief Holds the value of the gass stiffnes constant
	double		m_gasConstant;			// get from .pda file

	/// \brief Holds the value of the substance rest density
	double		m_substanceRestDensity;	// get from .pda file

	/// \brief Holds the value of the substance viscosity constant
	double		m_viscosityConstant;	// get from .pda file	

	/// \brief Holds the value of the surface color coefficient
	double		m_surfaceColorCoefficient;

	/// \brief Holds the value of the surface tension threshold
	double		m_surfaceTensionThreshold;

	/// \brief Holds the value of the surface tension coefficient
	double		m_surfaceTensionCoefficient; // get from .pda file

	/// \brief Holds the value of the interface color coefficient
	double		m_interfaceColorCoefficient;

	/// \brief Holds the value of the interface tension threshold
	double		m_interfaceTensionThreshold;

	/// \brief Holds the value of the interface tension coefficient
	double		m_interfaceTensionCoefficient;

	/// \brief Holds the pointer to the cell in which the particle is
	Cell*		m_cell;

	/// \brief The private function to check if a particle is close enough to be a neighbour
	/// @param[in] _particle the particle to be checked
	bool		isNear( Particle* _particle );
	
};

#endif