#ifndef CELL_H
#define CELL_H

#include <set>
#include <vector>
#include <iostream>

/// \file Cell.h
/// \brief Contains the Cell class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class Cell
/// \brief Holds the list of particles which are in a cell

class Particle;
class CellManager;

class Cell
{

public:
	 /// \brief Default ctor
	Cell();

	/// \brief Copy ctor
    /// @param[in] _cell the copied object
	Cell( const Cell& _cell );

	/// \brief Ctor
    /// @param[in] _x the array holding x coordinates of cell's lower and upper points
    /// @param[in] _y the array holding y coordinates of cell's lower and upper points
    /// @param[in] _z the array holding z coordinates of cell's lower and upper points
    /// @param[in] _i the cell id along x axis
	/// @param[in] _j the cell id along y axis
	/// @param[in] _k the cell id along z axis
	Cell( double* _x , double* _y , double* _z ,  int _i , int _j , int _k );	  

	/// \brief Dtor
	~Cell();

	/// \brief Adds a particle to the cell's particle list/set
    /// @param[in] _particle the particle to be added
	void addParticle( Particle* _particle );
	
	/// \brief Removes a particle from the cell's particle list/set
    /// @param[in] _particle the particle to be removed
	void removeParticle( Particle* _particle );	

	/// \brief Gets the id of the cell
    /// @param[in] id the reference to id array to be filled with cell's id values
	void getId (int (&id)[3] ) const;	

	/// \brief Returns the set of particles in the cell
    /// @returns The set of particles in the cell
 	const std::set<Particle* >& getParticles();
	

private:

	/// \brief Holds x coordinates of the cell's upper and lower points
	double m_x[2];

	/// \brief Holds y coordinates of the cell's upper and lower points
	double m_y[2];

	/// \brief Holds z coordinates of the cell's upper and lower points
	double m_z[2];
	
	/// \brief The set of particles in the cell
	std::set<Particle* > m_particle;

	/// \brief Holds the cells id
	int m_id[3];

};

#endif