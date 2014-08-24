#include "SPHSolver.h"
//-----------------------------------------------------------------------------
SPHSolver::SPHSolver( PDCFileWriter& pdcWriter, std::vector<Particle*>& particle, int startFrame, int frameCount ) : m_particle ( particle ) , m_PDCFileWriter( pdcWriter )
{

	//gets the radius for the neighbour search
	m_neighbourRadius = particle[0]->getNeighbourRadius();	
	
	double pi = 3.14;

	m_radiusSquared = m_neighbourRadius * m_neighbourRadius;

	// precalculates some values used to calculate Smoothing Kernels
	m_weightPoly = (315.0 / (64.0 * pi * pow(m_neighbourRadius, 9)));		

	m_weightPolyGradient = (-945.0 / (32.0 * pi * pow(m_neighbourRadius, 9)));

    m_weightPolyLaplacian = (-945.0 / (32.0 * pi * pow(m_neighbourRadius, 9)));

	m_weightPressureGradient = (-45.0 / (pi * pow(m_neighbourRadius, 6)));

    m_weightViscosityLaplacian = (45.0 / (pi * pow(m_neighbourRadius, 6)));

	// for each frame, we run simulation and result export into .pdc file
	for (int i = startFrame; i <=frameCount; i++ )
	{
		std::cout <<"frame count = "<< i <<std::endl;
		runSimulation();
		m_PDCFileWriter.writePDCFile( m_particle,  i );		
	}

}
//-----------------------------------------------------------------------------
double SPHSolver::WKernelPoly( double magnitudeSquared )
{
	if( magnitudeSquared <= m_radiusSquared ) 
	{
		return m_weightPoly * pow( m_radiusSquared - magnitudeSquared , 3 );
	}
	else
	{
		return 0.0;
	}
}
//-----------------------------------------------------------------------------
double SPHSolver::WKernelPolyGradient( double magnitudeVectorCoor, double magnitudeSquared ) 
{
    return m_weightPolyGradient * magnitudeVectorCoor * pow( m_radiusSquared - magnitudeSquared, 2 );   
}
//-----------------------------------------------------------------------------
double SPHSolver::WKernelPolyLaplacian( double magnitudeSquared )
{
    return m_weightPolyLaplacian * ( m_radiusSquared - magnitudeSquared ) * ( ( 3.0 * m_radiusSquared ) - ( 7.0 * magnitudeSquared ) );
}
//-----------------------------------------------------------------------------------------------------
double SPHSolver::WKernelPressureGradient( double magnitudeVectorCoor, double magnitude )
{
    return m_weightPressureGradient * pow( m_neighbourRadius - magnitude, 2 ) * ( magnitudeVectorCoor / magnitude );
}
//-----------------------------------------------------------------------------
double SPHSolver::WKernelViscosityLaplacian( double magnitude )
{
   return m_weightViscosityLaplacian * (m_neighbourRadius - magnitude);
}
//-----------------------------------------------------------------------------
void SPHSolver::runSimulation()
{
	std::vector<Particle*> neighbourParticle;
	
	// updates density and pressure values for each particle
	for( unsigned i = 0; i < m_particle.size(); ++i )
	{
		neighbourParticle.clear();
		//double r = m_particle[i]->getNeighbourRadius();
		m_particle[i]->getNeighbours( neighbourParticle );
		//std::cout<< "particle = " << i << " = neighbour count = " << neighbourParticle.size() << std::endl;
		//std::cout <<" particle = " << i <<std::endl;
		if (neighbourParticle.size() ) 
		{
			//std::cout<< "particle = " << i << " = neighbour count = " << neighbourParticle.size() << std::endl;
			updateDensityPressure( m_particle[i], neighbourParticle );
		}
		
	}
	
	// updates particle position and velocity
	for( unsigned i = 0; i < m_particle.size(); ++i )
	{
		neighbourParticle.clear();
		//std::cout <<" particle = " << i <<std::endl;
		m_particle[i]->getNeighbours( neighbourParticle );
		
		if (neighbourParticle.size() )
		{
			updatePositionVelocity( m_particle[i], neighbourParticle );
		}				
	}
}
//-----------------------------------------------------------------------------
void SPHSolver::updateDensityPressure( Particle* particle, std::vector<Particle*>& neighbourParticle )
{
	double magnitudeVector[3];
	double magnitudeSquared;  
	double density = 0;

	// it is going through the particle neighbours and calculates density for the particle
	for (unsigned n = 0; n < neighbourParticle.size(); n++)
	{
		magnitudeVector[0] = particle->getPosition(0) - neighbourParticle[n]->getPosition(0);
		magnitudeVector[1] = particle->getPosition(1) - neighbourParticle[n]->getPosition(1);
		magnitudeVector[2] = particle->getPosition(2) - neighbourParticle[n]->getPosition(2);

		magnitudeSquared = pow( magnitudeVector[0], 2 ) + pow( magnitudeVector[1], 2 ) + pow( magnitudeVector[2], 2 );
		
	
		density = density + neighbourParticle[n]->getMass() * WKernelPoly( magnitudeSquared );
		
	}
	
	// memorises new density value
	particle->setDensity( density );
			
	// calculates preassure value for the particle
	particle->calculatePressure();
	
}
//-----------------------------------------------------------------------------
void SPHSolver::updatePositionVelocity( Particle * particle, std::vector<Particle*>& neighbourParticle )
{
	double density = 0;
	
	double pressure[3] = { 0.0 , 0.0 , 0.0};
	double viscosity[3] = { 0.0 , 0.0 , 0.0};

	double magnitudeVector[3] = { 0.0 , 0.0 , 0.0};
	double magnitude;

	double massPerDensity = 0;

	double viscosityLaplacian = 0;

	double surfaceGradient[3] = { 0.0 , 0.0 , 0.0};
    double surfaceLaplacian = 0;

    double surfaceTensionGradient[3] = { 0.0 , 0.0 , 0.0};
    double surfaceTensionLaplacian = 0;

    double interfaceTensionGradient[3] = { 0.0 , 0.0 , 0.0};
    double interfaceTensionLaplacian = 0;
	
	// 
	for (unsigned n = 0; n < neighbourParticle.size(); n++)
	{	
		magnitudeVector[0] = particle->getPosition(0) - neighbourParticle[n]->getPosition(0);
		magnitudeVector[1] = particle->getPosition(1) - neighbourParticle[n]->getPosition(1);
		magnitudeVector[2] = particle->getPosition(2) - neighbourParticle[n]->getPosition(2);
		
		magnitude = sqrt( pow( magnitudeVector[0], 2 ) + pow( magnitudeVector[1], 2 ) + pow( magnitudeVector[2], 2 ) );

		massPerDensity = neighbourParticle[n]->getMass() / neighbourParticle[n]->getDensity();

		viscosityLaplacian = WKernelViscosityLaplacian( magnitude );
		
		// it is going through the particle neighbours and calculates pressure force, viscosity force, surface tension force and interface tension force
		for (unsigned i = 0; i < 3; i++)
		{	
			
			pressure[i] = pressure[i] + ( ( particle->getPressure() / pow( particle->getDensity() , 2 )) + 
												( neighbourParticle[n]->getPressure() / pow( neighbourParticle[n]->getDensity(), 2 ))  ) *
										     neighbourParticle[n]->getMass() * WKernelPressureGradient( magnitudeVector[i], magnitude );

			viscosity[i] = viscosity[i] + ( massPerDensity *                                       
										  ( neighbourParticle[n]->getVelocity(i) - particle->getVelocity(i)) *
											viscosityLaplacian );
			/*
			//accumulate surface tension-interface force component
            surfaceGradient[i] = WKernelPolyGradient( magnitudeVector[i], magnitude );            

            //accumulate surface gradient
            surfaceTensionGradient[i] = surfaceTensionGradient[i] + 
										( massPerDensity * particle->getSurfaceColorCoefficient() * surfaceGradient[i] );

            

            //accumulate interface gradient
            interfaceTensionGradient[i] = interfaceTensionGradient[i] + 
										  ( massPerDensity * particle->getInterfaceColorCoefficient() * surfaceGradient[i] );

            */
		}
		/*
		surfaceLaplacian = WKernelPolyLaplacian( magnitude );
		//accumulate surface laplacian
        surfaceTensionLaplacian = surfaceTensionLaplacian +
                                  ( massPerDensity * particle->getSurfaceColorCoefficient() * surfaceLaplacian );
		//accumulate interface laplacian
        interfaceTensionLaplacian = surfaceTensionLaplacian +
                                    ( massPerDensity * particle->getInterfaceColorCoefficient() * surfaceLaplacian );
									*/

	}
	
	double gravity[3]= { 0.0 , 0.0 , 0.0};
	double gravityAccel[3] = { 0.0, -9.8, 0.0 };
	/*
	double surfaceTension[3];
	double interfaceTension[3];

	double surfaceTensionGradientMagnitude = sqrt( pow ( surfaceTensionGradient[0], 2 ) + 
												   pow ( surfaceTensionGradient[1], 2 ) + 
												   pow ( surfaceTensionGradient[2], 2 ) );

	double interfaceTensionGradientMagnitude = sqrt( pow ( interfaceTensionGradient[0], 2 ) + 
													 pow ( interfaceTensionGradient[1], 2 ) + 
												     pow ( interfaceTensionGradient[2], 2 ) );
	*/
	double newAcceleration[3] = { 0.0 , 0.0 , 0.0};
	double newVelocity[3] = { 0.0 , 0.0 , 0.0};
	double newPosition[3] = { 0.0 , 0.0 , 0.0};
	
	// finalizes the calculations of the pressure force, viscosity force and calculates gravity force
	// calculates the new acceleration, velosity and position values
	for (unsigned i = 0; i < 3; i++) 
	{
		pressure[i] = (-1.0) * particle->getDensity() * pressure[i];
		viscosity[i] = particle->getViscosityConstant() * viscosity[i];
		gravity[i] = particle->getDensity() * gravityAccel[i];
		//std::cout << "preasure = " << i << " = " << pressure[i] << std::endl;
		/*		    
        if (surfaceTensionGradientMagnitude > particle->getSurfaceTensionThreshold())
        {
            surfaceTension[i] = (-1.0) * particle->getSurfaceTensionCoefficient() * 
								(surfaceTensionGradient[i]/surfaceTensionGradientMagnitude) * surfaceTensionLaplacian;
        }		
       
        if (interfaceTensionGradientMagnitude > particle->getInterfaceTensionThreshold())
        {
            interfaceTension[i] = (-1.0) * particle->getInterfaceTensionCoefficient() * 
								(interfaceTensionGradient[i]/interfaceTensionGradientMagnitude) * interfaceTensionLaplacian;
        }
		*/		
		newAcceleration[i] = ( pressure[i] + viscosity[i] + gravity[i] ) / particle->getDensity(); //+ surfaceTension[i] + interfaceTension[i]
		//why ir does not stop here? I was about to ask the same thing.
		//Semi-implicit_Euler
		newVelocity[i] = particle->getVelocity(i) + newAcceleration[i] * 0.004; //* timeStep = 1 second/25 frames = 0.04

		newPosition[i] = particle->getPosition(i) + newVelocity[i] * 0.004; //* timeStep = 1 second/25 frames = 0.04		
	}
	/*
	std::cout << "density = " << particle->getDensity() << std::endl;
	std::cout << "pressure = [ " << pressure[0] << " , " << pressure[1] << " , " << pressure[2] << " ]" << std::endl;
	std::cout << "viscosity = [ " << viscosity[0] << " , " << viscosity[1] << " , " << viscosity[2] << " ]" << std::endl;
	std::cout << "gravity = [ " << gravity[0] << " , " << gravity[1] << " , " << gravity[2] << " ]" << std::endl;
	std::cout << "newAcceleration = [ " << newAcceleration[0] << " , " << newAcceleration[1] << " , " << newAcceleration[2] << " ]" << std::endl;
	std::cout << "newVelocity = [ " << newVelocity[0] << " , " << newVelocity[1] << " , " << newVelocity[2] << " ]" << std::endl;
	std::cout << "oldPosition = [ " << particle->getPosition(0) << " , " << particle->getPosition(1) << " , " << particle->getPosition(2) << " ]" << std::endl;
	std::cout << "newPosition = [ " << newPosition[0] << " , " << newPosition[1] << " , " << newPosition[2] << " ]" << std::endl;
	std::cout << "-------------------------------------------------------------------------------------------------------------" << std::endl;
	*/
	//std::cout << "new position = " << newPosition[0] << " " << newPosition[1] << " " << newPosition[2] <<std::endl;
	particle->setNewPositionVelocity(newPosition[0], newPosition[1], newPosition[2],
									 newVelocity[0], newVelocity[1], newVelocity[2]);

}



//-----------------------------------------------------------------------------
