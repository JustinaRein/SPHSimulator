#ifndef CELL_MANAGER_H
#define CELL_MANAGER_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>

#include "Cell.h"

/// \file CellManager.h
/// \brief Contains the CellManager class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class CellManager
/// \brief Creates and manages cells

class Particle;

class CellManager
{
public:
	
	// \brief Ctor
    /// @param[in] _cellSize the radius used to set the number of cells 
    /// @param[in] _minPoint the coordinates of the lower point of the container
    /// @param[in] _maxPoint the coordinates of the upper point of the container    
	CellManager( double _cellSize , double* _minPoint, double* _maxPoint);

	/// \brief Dtor
	~CellManager();

	// \brief Gets the cell in which the particle is
    /// @param[in] _particlePosition the coordinates of the particle position
	Cell* getCell( double* _particlePosition );
	
	/// \brief Finds the neighbour cells for a given cell
    /// @param[in] cellId the reference to the cell id 
	/// @param[in] neighbourCell the reference to the neighbour cell vector
	void getNeighbourCells( const int (&cellId)[3] , std::vector<Cell*>& neighbourCell );

private:
	
	/// \brief Holds the number of cells along x, y and z axis
	int			m_numCells[3];

	/// \brief Holds the lengths of a cell along x, y and z axis
	double		m_delta[3];	

	/// \brief Holds the coordinates of the lower point of the container
	double		m_minPoint[3];

	/// \brief Holds all cells of the container
	std::vector<std::vector<std::vector<Cell*> > > m_cell;
};

#endif