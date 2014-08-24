#include "Particle.h"

CellManager* Particle::m_cellManager;

Particle::Particle() : m_cell( NULL )
{
	m_position[0] = m_position[1] = m_position[2] = 0.0;
}
//-----------------------------------------------------------------------------
void Particle::setId( void* input )
{
	m_id = *((unsigned*)input); //
	//std::cout << "Id = " << m_id << std::endl;
}

//-----------------------------------------------------------------------------
void Particle::setPosition( void* input )
{	
	m_position[0] = ((double*)input)[0];
	m_position[1] = ((double*)input)[1];
	m_position[2] = ((double*)input)[2];

	//std::cout << m_id << " Position = [" << m_position[0] << " , " << m_position[1] << " , " << m_position[2] << "]" << std::endl; 
}

//-----------------------------------------------------------------------------
void Particle::setVelocity( void* input )
{
	m_velocity[0] = ((double*)input)[0];
	m_velocity[1] = ((double*)input)[1];
	m_velocity[2] = ((double*)input)[2];
	
	//std::cout << "Velocity = [" << m_velocity[0] << " , " << m_velocity[1] << " , " << m_velocity[2] << "]" << std::endl;
}
//-----------------------------------------------------------------------------
void Particle::setMass( void* input )
{
	m_mass = *((double*)input);
	m_mass = 0.3;
	//std::cout << "mass = " << m_mass << std::endl;
}
//-----------------------------------------------------------------------------
void Particle::setRestDensity( void* input )
{
	m_substanceRestDensity = *((double*)input);	
	//std::cout << " RestDensity = " << m_substanceRestDensity << std::endl;
	m_substanceRestDensity = 998.29;	
}
//-----------------------------------------------------------------------------
void Particle::setViscosityConstant( void* input )
{
	m_viscosityConstant = *((double*)input);
	//std::cout << " ViscosityConstant = " << m_viscosityConstant << std::endl;
	m_viscosityConstant = 0.0;
}
//-----------------------------------------------------------------------------
void Particle::setNeighbourRadius( void* input )
{
	m_neighbourRadius = *((double*)input);
	m_neighbourRadius = 0.04;
	//std::cout << " NeighbourRadius = " << m_neighbourRadius << std::endl;
}
//-----------------------------------------------------------------------------
void Particle::setGasConstant( void* input )
{
	m_gasConstant = *((double*)input); 
	//std::cout << " GasConstant = " << m_gasConstant << std::endl;
	m_gasConstant = 0.1;
}

//-----------------------------------------------------------------------------
void Particle::setSurfaceTensionConstant( void* input )
{
	m_surfaceTensionCoefficient = *((double*)input);
}
//-----------------------------------------------------------------------------
Particle::Setter Particle::getSetterMethod( std::string _name )
{
	//calls the right method according to an attribute name read from .pda file
		if( _name == "id" ) { return &Particle::setId; }
		
		if( _name == "position" ) { return &Particle::setPosition; }

		if( _name == "velocity" ) { return &Particle::setVelocity; }

		if( _name == "mass" ) { return &Particle::setMass; }	

		if( _name == "viscosity" ) { return &Particle::setViscosityConstant; }

		if( _name == "restDensity" ) { return &Particle::setRestDensity; }

		if( _name == "radius" ) { return &Particle::setNeighbourRadius; }

		if( _name == "stickiness" ) { return &Particle::setGasConstant; }

		if( _name == "surfaceTension" ) { return &Particle::setSurfaceTensionConstant; }
		
		//std::cout << "No setter method for " << name << std::endl;
		return NULL;
}
//-----------------------------------------------------------------------------
void Particle::setNewPositionVelocity( double _pX , double _pY , double _pZ, double _vX, double _vY, double _vZ )
{
	double position[3];

	position[0] = _pX;
	position[1] = _pY;
	position[2] = _pZ;


	// as the particle has moves, we need to check, if it is still in the same cell
	Cell *c = m_cellManager->getCell( position );

	if ( c != NULL)
	{
		m_position[0] = _pX;
		m_position[1] = _pY;
		m_position[2] = _pZ;

		m_velocity[0] = _vX;
		m_velocity[1] = _vY;
		m_velocity[2] = _vZ;

		if( c != m_cell )
		{
			// if the particle is not in the cell anymore, it is removed from the cell's particle list
			m_cell->removeParticle( this );
			// the particle is added to the new cell in which it is now
			c->addParticle( this );
			// the particle memorised the cell
			m_cell = c;
		}
	}
	else
	{
		// if the new calculated cell does not exist, that means, the particle tries to escape the container
		// reverce the direction of velocity values
		m_velocity[0] =  m_velocity[0] * (-0.01);
		m_velocity[1] =  m_velocity[1] * (-0.01);
		m_velocity[2] =  m_velocity[2] * (-0.01);

	}

}

//-----------------------------------------------------------------------------
void Particle::setCell( Cell* _cell )  
{
	m_cell = _cell;
	
	//std::cout << m_id << " Particle is in Cell = [" << m_cell->getIdX() << " , " << m_cell->getIdY() << " , " << m_cell->getIdZ() << "]" << std::endl;
}
//-----------------------------------------------------------------------------
void Particle::setDensity( double _density )
{
	m_density = _density;
}
//-----------------------------------------------------------------------------
void Particle::calculatePressure()
{
	//modified ideal ga state equation
	m_pressure = m_gasConstant * (m_density - m_substanceRestDensity);
}
//-----------------------------------------------------------------------------
void Particle::getNeighbours( std::vector<Particle*>& _particleNeighbour )
{
	std::vector<Cell*> neighbourCell;
	
	int id[3];

	m_cell->getId( id ); 
	//std::cout<< "particle "<< m_id << "cell = [ " << id[0] << " , " << id[1] << " , " << id[2] << " ] =  "; //<<std::endl;
	// gets neighbour cells of the cell in which particle is
	m_cellManager->getNeighbourCells( id , neighbourCell );
	unsigned len = neighbourCell.size();
	
	for(unsigned i = 0; i < len; i++)
	{
		//gets the particle list from each neigbour cell
		const std::set<Particle* >& p = neighbourCell[i]->getParticles();	

		std::set<Particle*>::const_iterator pIt = p.begin();
		std::set<Particle*>::const_iterator pItEnd = p.end();

		while( pIt != pItEnd )
		{
			Particle *pt = *pIt;
			// checks if a particle from the list is near enough to be memorised as a neighbour
			if( this != pt && isNear( pt ) )
			{
				_particleNeighbour.push_back( pt );
			}
			pIt++;
		}
	}

	// checks the particle list of the cell in which the particle is
	const std::set<Particle* >& p = m_cell->getParticles();	

	std::set<Particle*>::const_iterator pIt = p.begin();
	std::set<Particle*>::const_iterator pItEnd = p.end();

	while( pIt != pItEnd )
	{
		Particle *pt = *pIt;
		// checks if a particle from the list is near enough to be memorised as a neighbour
		if( this != pt && isNear( pt ) )
		{
			_particleNeighbour.push_back( pt );
		}
		pIt++;
	}

}
//-----------------------------------------------------------------------------
bool Particle::isNear( Particle* _particle )
{
	//calculates the distance between the particle and the potential neighbour particle
	// sqtr is not used, because it is expencieve to calculate
	double L2 = pow( (_particle->m_position[0] - m_position[0]), 2 ) + 
				pow( (_particle->m_position[1] - m_position[1]), 2 ) + 
				pow( (_particle->m_position[2] - m_position[2]), 2 );

	// the distance has to be smaller than the given radius for the neighbour search
	return L2 < (m_neighbourRadius * m_neighbourRadius);
}

//-----------------------------------------------------------------------------
unsigned Particle::getId()
{
	//std::cout << " id = " << m_id << std::endl;
	return m_id;
}
//-----------------------------------------------------------------------------
double Particle::getPosition( unsigned _i) const
{
	//std::cout << " position = " << i << " = " << m_position[i] << std::endl;
	return m_position[_i];
}
//-----------------------------------------------------------------------------
double Particle::getVelocity( unsigned _i ) const
{
	//std::cout << " velocity = " << i << " = " << m_velocity[i] << std::endl;
	return m_velocity[_i];
}
//-----------------------------------------------------------------------------
double Particle::getMass() const
{
	//std::cout << " mass = " << m_mass << std::endl;
	return m_mass;
}
//-----------------------------------------------------------------------------
double Particle::getDensity() const
{
	//std::cout << " Density = " << m_density << std::endl;
	return m_density;
}
//-----------------------------------------------------------------------------
double Particle::getPressure() const
{
	//std::cout << " Pressure = " << m_pressure << std::endl;
	return m_pressure;
}
//-----------------------------------------------------------------------------
double Particle::getSubstanceRestDensity() const
{
	//std::cout << " RestDensity = " << m_substanceRestDensity << std::endl;
	return m_substanceRestDensity;
}	
//-----------------------------------------------------------------------------
double Particle::getViscosityConstant() const
{
	//std::cout << " viscosity = " << m_viscosityConstant << std::endl;
	return m_viscosityConstant;
}
//-----------------------------------------------------------------------------
double Particle::getNeighbourRadius() const
{		
	return m_neighbourRadius;
}
//-----------------------------------------------------------------------------
// surface tension
double Particle::getSurfaceColorCoefficient() const
{
	return m_surfaceColorCoefficient;
}
//-----------------------------------------------------------------------------
double Particle::getSurfaceTensionThreshold() const
{
	return m_surfaceTensionThreshold;
}
//-----------------------------------------------------------------------------
double Particle::getSurfaceTensionCoefficient() const
{
	return m_surfaceTensionCoefficient;
}
//-----------------------------------------------------------------------------
// interface tension
double Particle::getInterfaceColorCoefficient() const
{
	return m_interfaceColorCoefficient;
}
//-----------------------------------------------------------------------------
double Particle::getInterfaceTensionThreshold() const
{
	return m_interfaceTensionThreshold;
}
//-----------------------------------------------------------------------------
double Particle::getInterfaceTensionCoefficient() const
{
	return m_interfaceTensionCoefficient;
}
//-----------------------------------------------------------------------------