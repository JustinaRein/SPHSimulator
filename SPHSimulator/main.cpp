#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include "PDCFileWriter.h"
#include "PDAFileReader.h"

#include "Particle.h"
#include "CellManager.h"
#include "SPHSolver.h"

int main()
{
	char* inputFile = "C:/Programming/fluid_simulator/nParticleShape1.1.pda";
	
	double minPoint[3] = { -2, 0, -2};
	double maxPoint[3] = {2, 4, 2};
	
	double cellSize = 0.05;
	
	CellManager cellManger(cellSize, minPoint, maxPoint);
	Particle::m_cellManager = &cellManger;
	
	std::vector<Particle*> particle;

	PDAFileReader PDAFileReader( inputFile , particle , cellManger );

	int startFrame = (int) PDAFileReader.getStartFrame();
	int frameCount = 150; //(int) PDAFileReader.getCacheWidth();

	PDCFileWriter pdcWriter( inputFile );

	SPHSolver SPHSolver( pdcWriter , particle , startFrame , frameCount );
		
	return  EXIT_SUCCESS;
}