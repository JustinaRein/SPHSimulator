#ifndef FLUID_PARTICLE_H
#define FLUID_PARTICLE_H


//#include <vector>

//using namespace std;
//#include "Cell.h"

class Cell;

class FluidParticle
{

public:

	FluidParticle();
	FluidParticle(int id, double vec_x, double vec_y, double vec_z);

	~FluidParticle();

	int get_fluid_particle_id() const
	{
		return m_id;
	}

	double get_fluid_particle_vecx() const
	{
		return m_position[0];
	}

	double get_fluid_particle_vecy() const
	{
		return m_position[1];
	}

	double get_fluid_particle_vecz() const
	{
		return m_position[2];
	}

	void set_cell( Cell *cell )
	{
		m_cell = cell;
	}
		
private:

	int m_id;
	double m_position[3];
	Cell*	m_cell;
	//Cell m_cell;
};

#endif