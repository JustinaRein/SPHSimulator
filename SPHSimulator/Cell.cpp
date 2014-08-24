#include "Cell.h"

//-----------------------------------------------------------------------------
Cell::Cell()
{
	m_id[0] = m_id[1] = m_id[2] = 0;
	//std::cout << "Cell::Cell()" << std::endl;
}

//-----------------------------------------------------------------------------
Cell::Cell( const Cell&  _cell )
{
	m_id[0] = _cell.m_id[0];
	m_id[1] = _cell.m_id[1];
	m_id[2] = _cell.m_id[2];

	//std::cout << "Copying cell [" << m_id[0] << "][" << m_id[1] << "][" << m_id[2] << "]" <<std::endl;
}

//-----------------------------------------------------------------------------
Cell::Cell( double* _x , double* _y , double* _z ,  int _i , int _j , int _k )
{
	m_id[0] = _i;
	m_id[1] = _j;
	m_id[2] = _k;

	m_x[0] = _x[0];
	m_x[1] = _x[1];

	m_y[0] = _y[0];
	m_y[1] = _y[1];

	m_z[0] = _z[0];
	m_z[1] = _z[1];

	//std::cout << "Creating cell [" << m_id[0] << "][" << m_id[1] << "][" << m_id[2] << "] "; // <<std::endl;
	//std::cout	<< "[" << m_x[0] << " , " << m_y[0] << " , " << m_z[0] << " ] ["
	//			<< m_x[1] << " , " << m_y[1] << " , " << m_z[1] << " ]" << std::endl;
}
//-----------------------------------------------------------------------------
void Cell::addParticle( Particle* _particle )
{
	m_particle.insert( _particle );

	//std::cout << "cell [ " << m_id[0] << " " << m_id[1] << " " << m_id[2] << " ] = " << m_particle.size() << std::endl;
}

//-----------------------------------------------------------------------------
void Cell::removeParticle( Particle* _particle )
{
	m_particle.erase( _particle );
}

//-----------------------------------------------------------------------------
void Cell::getId (int (&id)[3] ) const
	{
		id[0] = m_id[0];
		id[1] = m_id[1];
		id[2] = m_id[2];
	}
//-----------------------------------------------------------------------------
const std::set<Particle* >& Cell::getParticles()
{
	return m_particle;
}
//-----------------------------------------------------------------------------
