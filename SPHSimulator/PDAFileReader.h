#ifndef PDA_FILE_READER_H
#define PDA_FILE_READER_H

#include "stdlib.h"
#include "stdio.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Particle.h"
#include "Cell.h"
#include "CellManager.h"

/// \file PDAFileReader.h
/// \brief Contains the PDAFileReader class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class PDAFileReader
/// \brief Reads a .pda file and gets particle and scene data from Maya


class PDAFileReader
{
public:
	/// \brief Default ctor
	PDAFileReader();

	/// \brief Ctor
    /// @param[in] _inputfile the file to be read
    /// @param[in] _particle the reference to the particle vector to be filled
    /// @param[in] _cellManager the reference to the existing cell manager to find a cell for the particle     
	PDAFileReader(char* _inputfile , std::vector<Particle*>& _particle , CellManager &cellManager ); 

	/// \brief Dtor
	~PDAFileReader();
		
	/// \brief Gets the start frame for the simulation
	int getStartFrame() const;
	
	/// \brief Gets the number of cache files to be writen
	int getCacheWidth() const;
		
private:

	/// \brief The number of particles read from the .pda file
	int m_particleCount;

	/// \brief The start frame for the simulation	
	int m_startFrame;

	/// \brief The number of cache files to be writen
	int m_cacheWidth;
	
	/// \brief The cell manager object
	CellManager m_cellManager;

	/// \brief The vector of particles
	std::vector<Particle*>& m_particle;
		
};

#endif
