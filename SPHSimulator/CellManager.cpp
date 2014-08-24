#include "CellManager.h"
//-----------------------------------------------------------------------------
CellManager::CellManager(double _cellSize , double* _minPoint, double* _maxPoint )
{
	//calculates thenumbers of cells along x, y and z axis
	m_numCells[0] = (int) floor((_maxPoint[0] - _minPoint[0])/ _cellSize);
	m_numCells[1] = (int) floor((_maxPoint[1] - _minPoint[1])/ _cellSize);
	m_numCells[2] = (int) floor((_maxPoint[2] - _minPoint[2])/ _cellSize);

	//calculates the lengths of cells along x, y and z axis
	m_delta[0] = (double) ((_maxPoint[0] - _minPoint[0])/m_numCells[0]);
	m_delta[1] = (double) ((_maxPoint[1] - _minPoint[1])/m_numCells[1]);
	m_delta[2] = (double) ((_maxPoint[2] - _minPoint[2])/m_numCells[2]);

	//memorises the coordinates of the lower point of the container
	m_minPoint[0] = _minPoint[0];
	m_minPoint[1] = _minPoint[1];
	m_minPoint[2] = _minPoint[2];


	double x[2] , y[2] , z[2];

	m_cell.resize( m_numCells[0] );

	//creates cells
	for (int i = 0; i < m_numCells[0]; ++i )
	{
		m_cell[i].resize( m_numCells[1] );

		for (int j = 0; j < m_numCells[1]; ++j )
		{
			m_cell[i][j].resize( m_numCells[2] );

			for (int k = 0; k < m_numCells[2] ; ++k )
			{
				x[0] = _minPoint[0] + _cellSize * i;
				y[0] = _minPoint[1] + _cellSize * j;
				z[0] = _minPoint[2] + _cellSize * k;

				x[1] = x[0] + _cellSize;
				y[1] = y[0] + _cellSize;
				z[1] = z[0] + _cellSize;

				m_cell[i][j][k] = new Cell( x , y , z , i , j , k );
				//std::cout << i << " " << j << " " << k << std::endl;				
			}
		}
	}
}

//-----------------------------------------------------------------------------
CellManager::~CellManager()
{
	// for some reason the program stops running when this section is active
	/*for (int i = 0; i < m_numCells[0]; ++i )
	{
		for (int j = 0; j < m_numCells[1]; ++j )
		{
			for (int k = 0; k < m_numCells[2] ; ++k )
			{
					delete m_cell[i][j][k];			
							
			}
		}
	}*/
}

//-----------------------------------------------------------------------------
void CellManager::getNeighbourCells( const int (&cellId)[3] , std::vector<Cell*>& neighbourCell )
{
	
	//Middle plane (z coordinate the same)
	if ((cellId[0] - 1) >= 0) 
	{
		//(x_left) = true;
		neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1]][cellId[2]]);		
	}
	
	if ((cellId[0] + 1) < m_numCells[0])
	{
		//x_right = true;
		neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1]][cellId[2]]);
	}
	
	if ((cellId[1] + 1) < m_numCells[1])		
	{
		//y_up = true;
		neighbourCell.push_back(m_cell[cellId[0]][cellId[1] + 1][cellId[2]]);
		if ((cellId[0] - 1) >= 0)				//(x_left)
		{
			neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] + 1][cellId[2]]);
		}
		if ((cellId[0] + 1) < m_numCells[0])		//(x_right)
		{
			neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] + 1][cellId[2]]);
		}
	}

	if ( (cellId[1] - 1) >= 0 )		
	{
		//y_down = true;
		neighbourCell.push_back(m_cell[cellId[0]][cellId[1] - 1][cellId[2]]);
		if ((cellId[0] - 1) >= 0)				//(x_left)
		{
			neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] - 1][cellId[2]]);
		}
		if ((cellId[0] + 1) < m_numCells[0])		//(x_right)
		{
			neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] - 1][cellId[2]]);
		}
	}

	//Front plane ---------------------------------
	if ((cellId[2] - 1) >= 0)
	{
		//z_close = true;
		neighbourCell.push_back(m_cell[cellId[0]][cellId[1]][cellId[2] - 1]);
		if ((cellId[0] - 1) >= 0)				//(x_left)
		{
			neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1]][cellId[2] - 1]);
			if ((cellId[1] + 1) < m_numCells[1])			//(y_up)
			{
				neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] + 1][cellId[2] - 1]);
			}
			if ((cellId[1] - 1) >= 0 )	//(y_down)
			{
				neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] - 1][cellId[2] - 1]);
			}
		} 
		if ((cellId[0] + 1) < m_numCells[0])		//(x_right)
		{
			neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1]][cellId[2] - 1]);
			if ((cellId[1] + 1) < m_numCells[1])			//(y_up)
			{
				neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] + 1][cellId[2] - 1]);
			}
			if ((cellId[1] - 1) >= 0 )	//(y_down)
			{
				neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] - 1][cellId[2] - 1]);
			}
		}

		if ((cellId[1] + 1) < m_numCells[1])				//(y_up)
		{
			neighbourCell.push_back(m_cell[cellId[0]][cellId[1] + 1][cellId[2] - 1]);
		}

		if ((cellId[1] - 1) >= 0 )		//(y_down)
		{
			neighbourCell.push_back(m_cell[cellId[0]][cellId[1] - 1][cellId[2] - 1]);
		}
	}	

	//Back plain -------------------------------------	
	if ((cellId[2] + 1) < m_numCells[2])
	{
		//z_far = true;
		neighbourCell.push_back(m_cell[cellId[0]][cellId[1]][cellId[2] + 1]);
		if ((cellId[0] - 1) >= 0)				//(x_left)
		{
			neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1]][cellId[2] + 1]);
			if ((cellId[1] + 1) < m_numCells[1])			//(y_up)
			{
				neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] + 1][cellId[2] + 1]);
			}
			if ((cellId[1] - 1) >= 0 )	//(y_down)
			{
				neighbourCell.push_back(m_cell[cellId[0] - 1][cellId[1] - 1][cellId[2] + 1]);
			}
		}

		if ((cellId[0] + 1) < m_numCells[0])		//(x_right)
		{
			neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1]][cellId[2] + 1]);
			if ((cellId[1] + 1) < m_numCells[1])			//(y_up)
			{
				neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] + 1][cellId[2] + 1]);
			}
			if ((cellId[1] - 1) >= 0 )	//(y_down)
			{
				neighbourCell.push_back(m_cell[cellId[0] + 1][cellId[1] - 1][cellId[2] + 1]);
			}
		}

		if ((cellId[1] + 1) < m_numCells[1])				//(y_up)
		{
			neighbourCell.push_back(m_cell[cellId[0]][cellId[1] + 1][cellId[2] + 1]);
		}

		if ((cellId[1] - 1) >= 0 )		//(y_down)
		{
			neighbourCell.push_back(m_cell[cellId[0]][cellId[1] - 1][cellId[2] + 1]);
		}
	}

	
}

//-----------------------------------------------------------------------------
Cell* CellManager::getCell( double* _particlePosition )
{
	//gets the cell id for the particle
	int idX = (int)floor((_particlePosition[0] - m_minPoint[0])/m_delta[0]);
	int idY = (int)floor((_particlePosition[1] - m_minPoint[1])/m_delta[1]);
	int idZ = (int)floor((_particlePosition[2] - m_minPoint[2])/m_delta[2]);
	
	//checks if the cell id exists, if it returns NULL, that means the particle tries to leave the container
	if( idX < 0 || idY < 0 || idZ < 0 || idX >= m_numCells[0] || idY >= m_numCells[1] || idZ >= m_numCells[2] )
	{
		return NULL;
	}

	return m_cell[idX][idY][idZ];
}
//-----------------------------------------------------------------------------

