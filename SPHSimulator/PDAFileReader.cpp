
#include "PDAFileReader.h"
//-----------------------------------------------------------------------------
PDAFileReader::PDAFileReader(char * _inputfile , std::vector<Particle*>& _particle , CellManager &cellManager ) : m_particle( _particle ) , m_cellManager( cellManager )
{
	std::vector <std::string> attributeName;
	std::vector <char> attributeType;
	//int particleCount = 0;
	int attributeCount = 0;
	std::vector <std::string> particleDataString;	
	
	//Open and read input file-----------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------
	std::ifstream inFile(_inputfile);
	std::string line;
	int n = 0;

	while (getline(inFile, line))
	{
		n++;
		//the second line is the list of the attribute names
		if (n == 2)
		{			
			std::string name;
			std::istringstream iss(line);
			while (iss >> name)
			{
				attributeName.push_back(name);
				attributeCount++;
				//cout<<"attribute_name = "<< a_name<< endl;
			}			
		}
		//the fourth line is the list of attribute types (could be I - int, R - real, V - vector)
		if (n == 4)
		{			
			char aType;
			std::istringstream iss(line);			
			while (iss >> aType)
			{				
				attributeType.push_back(aType);				
				//cout << "attribute_type = " << a_type << endl;
			}
		}
		// the number of particles is got from the fifth line
		if (n == 5)
		{
			char * wordSimbols;
			std::string word;
			std::istringstream iss(line);
			while (iss.good())
			{
				iss >> word;
			}
			wordSimbols = new char [word.size()];			
			strcpy (wordSimbols, word.c_str()); 			
			m_particleCount = atoi(wordSimbols);			
			//cout<< "number of particles = "<< particle_count<<endl;
		}
		// after the fifth line every line represents a paticle data
		if (n > 6)
		{
			particleDataString.push_back(line);		
		}
	}	
	
	inFile.close(); 

	//---------------------------------------------------------------------------------------------------
	// now we will sort all the data
	//---------------------------------------------------------------------------------------------------
	
	Cell *c;
	Particle::Setter s;
	Particle *particle;

	int intData;
	double realData;
	double vecData[3];

	//-------GET ALL ATTRIBUTES FOR EACH PARTICLE-----------
	// (not all attributes will be used for the simulation, as Maya has its own value system, which is not fully understood)
	// -----------------------------------------------------
	// id 				particle id
	// casheWidth		the number of frames to be cashed
	// -----------------------------------------------------
	// mass 			particle mass		
	// viscosity		substance viscosity
	// rest density		substance rest density
	// radius			radius for neigbours
	// stickiness		in another words it is a gas constant
	// start frame		the start frame of the simulation 
	// surface tension  substance surface tension
	//-------------------------------------------------------
	// position			position (x y z) of a particle
	// velosity			velocity (x y z) of a particle

	//fills particle vector with particles and gets all data from the particle_data_string vector (see above)
	m_particle.resize( m_particleCount );

	for (int i = 0; i < m_particleCount; i++)
	{
		particle = new Particle();

		std::istringstream iss(particleDataString[i]);
				
		for (int j = 0; j < attributeCount; j++)
		{		
			if (attributeType[j] == 'I')
			{				
				iss >> intData;
				if (attributeName[j] == "id")
				{					
					s = Particle::getSetterMethod("id");
					if( s != NULL )
					{
						unsigned id = intData;
						(particle->*s)(&id);
					}						
				}

				if (attributeName[j] == "cacheWidth")
				{					
					m_cacheWidth = intData;
					//std::cout<<m_cacheWidth<<std::endl;
				}
			}			
			if (attributeType[j] == 'R') 
			{
				iss >> realData;
				
				if (attributeName[j] == "mass")
				{					
					s = Particle::getSetterMethod("mass");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}
				else if (attributeName[j] == "viscosity")
				{					
					s = Particle::getSetterMethod("viscosity");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}
				else if (attributeName[j] == "restDensity")
				{					
					s = Particle::getSetterMethod("restDensity");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}
				else if (attributeName[j] == "radius")
				{	
					s = Particle::getSetterMethod("radius");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}
				else if (attributeName[j] == "stickiness")
				{					
					s = Particle::getSetterMethod("stickiness");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}
				else if (attributeName[j] == "startFrame")
				{					
					m_startFrame = (int)realData; 				
				}
				else if (attributeName[j] == "surfaceTension")
				{					
					s = Particle::getSetterMethod("surfaceTention");
					if( s != NULL )
					{
						(particle->*s)(&realData);
					}					
				}

			}
			if (attributeType[j] == 'V')
			{
				iss >> vecData[0];
				iss >> vecData[1];
				iss >> vecData[2];

				if (attributeName[j] == "position")
				{					
					s = Particle::getSetterMethod("position");
					if( s != NULL )
					{
						(particle->*s)(vecData);
					}
					// as we know the position of each particle, 
					// we could find the cell in which a particle is at the start of simulation
					c = m_cellManager.getCell( vecData );
					//std::cout<< i << " particle = cell [ " << c->getIdX() << " , " << c->getIdY() << " , " << c->getIdZ() << " ] " << std::endl;
					particle->setCell( c );						
					c->addParticle( particle );  
				}
				else if (attributeName[j] == "velocity") 
				{					
					s = Particle::getSetterMethod("velocity");
					if( s != NULL )
					{
						(particle->*s)(vecData);
					}					
				}
			} 			
		}		
		m_particle[i] = particle;
	}
}
//-----------------------------------------------------------------------------
PDAFileReader::~PDAFileReader()
{
	for( int i = 0; i < m_particleCount; ++i )
	{
		//std::cout << "Destroying Particle " << m_particle[i]->getId() << std::endl;
		delete m_particle[i];
	}

}
//-----------------------------------------------------------------------------
int PDAFileReader::getStartFrame() const
{
	return m_startFrame;
}
//-----------------------------------------------------------------------------
int PDAFileReader::getCacheWidth() const
{
	return m_cacheWidth;
}
//-----------------------------------------------------------------------------