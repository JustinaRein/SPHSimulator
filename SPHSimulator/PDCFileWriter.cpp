/***********************************
This program is written by Peter J. Lu (Harvard University) and
Amit Chourasia with help from John Moreland (San Diego Supercomputer Center)
For bugs and other info contact the authors.
This information should be left unmodified with the source.
The program is provided "AS-IS". Use it at your own risk. Authors are not liable for any damages resulting from the use of the program
Last Modified: January 29, 2005

*///////////////////////////////////

#include "PDCFileWriter.h"

// Notes:
// Uses all C++ filestreams for outputs (C-style fwrite operations tended to get mixed up with byte alignment issues, writing junk to the file.
// The C++ filestreams do not seem to suffer from this problem with the Intel Compiler under WindowsXP.
// Byte order is assumed to be Big Endian for Maya, so that all numerical data must have its byte order reversed using the swap_int and swap_double functions.
// The swap functions reverse byte order in-place (pass variables by reference), so make sure to pass copies of data that should not be modified
// PDCHeader class keeps track of the header information only, and has a member function to properly write it to the disk file. This must be
// declared/constructed first, and the write_PDCHeader member function called first. The data is stored as private variables, which cannot be modified (though
// various get_ functions return the actual data without offering the possibility of it being modified. In particular, the header class holds both the number
// of attributes, and the number of particles, which can be useful for sizing loops.
// After that, for each attribute, the write_attrib function must be called, followed immediately a list of the appropriate data called by the correct
// function for that data type. Available from: http://users.sdsc.edu/~amit/web/book/export/html/21
//-----------------------------------------------------------------------------
PDCFileWriter::PDCFileWriter( char* fileName ): m_formatVersion(1), m_byteOrder(1), m_extra1(0), m_extra2(0)
{
	m_format[0]='P';
	m_format[1]='D';
	m_format[2]='C';
	m_format[3]=' ';	
	
	//std::cout<<"particle count = " << m_particleCount << std::endl;

	m_attributeCount = 2; 
	
	//preapare output file name	
	string outputfileName;	
	istringstream str(fileName);
	getline (str, outputfileName, '.');
	
	m_outputFile = new char [outputfileName.size()];
	strcpy (m_outputFile, outputfileName.c_str());
}

//-----------------------------------------------------------------------------
PDCFileWriter::~PDCFileWriter()
{
	//delete [] m_outputFile;
}

//-----------------------------------------------------------------------------
void PDCFileWriter::writePDCFile( std::vector<Particle*>& particle, int frameNumber)
{
	m_particleCount = particle.size();
	 //Append .pdc to the filename provided
	ostringstream outfilenametext;
	outfilenametext << m_outputFile << "." << frameNumber * 250 <<".pdc";
	string outfilename = outfilenametext.str();

	//create C++ file stream to handle output
	ofstream outFile;
	outFile.open(outfilename.c_str(),ios::out | ios::binary);
	
	writePDCHeader(outFile);

	
	int attributeNameLength;

	// --- write id ---------------------

	char attributeName1[] = "id";

	attributeNameLength = (int) strlen( "id" );

	writeAttribute(outFile, attributeNameLength, attributeName1, 1);
	
	for (int i = 0; i < m_particleCount; i++)
	{
		write1intscalar(outFile,particle[i]->getId());	
	}
	
	// --- write position ---------------------

	char attributeName2[] = "position";
	//cout << attributeName <<endl;
	attributeNameLength = (int) strlen( "position" );
	//cout << attributeNameLength <<endl;
	writeAttribute(outFile, attributeNameLength, attributeName2, 5);
	
	for (int i = 0; i < m_particleCount; i++)
	{
		write3vector(outFile, particle[i]->getPosition(0), particle[i]->getPosition(1), particle[i]->getPosition(2));		
	}
	
	//for (int i = 0; i < attributeCount; i++)
	//{
		/*int attribute_name_length = (int) strlen (pda_data.m_particle[0][i].attribute_name);
		//int attribute_name_length = (int) strlen (pda_data.m_particle[0][i].attribute_name);
		
		
		write_attrib(outFile, attribute_name_length, pda_data.m_particle[0][i].attribute_name, pda_data.m_particle[0][i].attribute_type);
		
		for (int j = 0; j < particle_count; j++)
		{
			if (pda_data.m_particle[j][i].attribute_type == 1)
				write_1intscalar(outFile,pda_data.m_particle[j][i].attribute_intType);
			if (pda_data.m_particle[j][i].attribute_type == 3)
				write_1scalar(outFile,pda_data.m_particle[j][i].attribute_realType);
			if (pda_data.m_particle[j][i].attribute_type == 5)
				write_3vector(outFile, pda_data.m_particle[j][i].attribute_vectorType[0], pda_data.m_particle[j][i].attribute_vectorType[1], pda_data.m_particle[j][i].attribute_vectorType[2]);		
		}
		*/
	//}
	outFile.close();
}

//-----------------------------------------------------------------------------
void PDCFileWriter::writePDCHeader(ofstream &outfile)
{
	//create temporary versions of the variables to byte-swap

	int temp_formatVersion = getFormatVersion();
	int temp_byteOrder = getByteOrder();
	int temp_extra1 = getExtra1();
	int temp_extra2 = getExtra2();
	int temp_numParticles = getParticleCount();
	int temp_numAttributes = getAttributeCount();

	//do swap of byte order

	swapInt((char*) &temp_formatVersion);
	swapInt((char*) &temp_byteOrder);
	swapInt((char*) &temp_extra1);
	swapInt((char*) &temp_extra2);
	swapInt((char*) &temp_numParticles);
	swapInt((char*) &temp_numAttributes);

	//write out to file

	for(int i=0;i<4;i++) 
	{
		outfile.put( m_format[i]);
	}

	outfile.write((char*) &temp_formatVersion, sizeof(int));
	outfile.write((char*) &temp_byteOrder, sizeof(int));
	outfile.write((char*) &temp_extra1, sizeof(int));
	outfile.write((char*) &temp_extra2, sizeof(int));
	outfile.write((char*) &temp_numParticles, sizeof(int));
	outfile.write((char*) &temp_numAttributes, sizeof(int));
}

//-----------------------------------------------------------------------------
char PDCFileWriter::getFormat0() const
{
	return m_format[0];
}

//-----------------------------------------------------------------------------
char PDCFileWriter::getFormat1() const
{
	return m_format[1];
}

//-----------------------------------------------------------------------------
char PDCFileWriter::getFormat2() const
{
	return m_format[2];
}

//-----------------------------------------------------------------------------
char PDCFileWriter::getFormat3() const
{
	return m_format[3];
}

//-----------------------------------------------------------------------------
int PDCFileWriter::getFormatVersion() const
{
	return m_formatVersion;
}
//-----------------------------------------------------------------------------
int PDCFileWriter::getByteOrder() const
{
	return m_byteOrder;
}

//-----------------------------------------------------------------------------
int PDCFileWriter::getExtra1() const
{
	return m_extra1;
}

//-----------------------------------------------------------------------------
int PDCFileWriter::getExtra2() const
{
	return m_extra2;
}

//-----------------------------------------------------------------------------
int PDCFileWriter::getParticleCount() const
{
	return m_particleCount;
}
//-----------------------------------------------------------------------------
int PDCFileWriter::getAttributeCount() const
{
	return m_attributeCount;
}

//-----------------------------------------------------------------------------
void PDCFileWriter::swapInt(char* int_array)
{
	// Swaps byte order of 4-element character array.
	// Pass to this function a number which is an integer cast to a char*, allowing byte order to be reversed.
	// note that by passing a reference to this function, the original variable itself is modified and is no longer
	// useful as an integer for any numerical or comparison operation
	char swap;

	for(int i=0;i<2;i++) 
	{
		swap = int_array[3-i];
		int_array[3-i]=int_array[i];
		int_array[i]=swap;
	}
}

//-----------------------------------------------------------------------------
void PDCFileWriter::swapDouble(char* double_array)
{
	// Swaps byte order of 8-element character array.
	// Pass to this function a number which is a double cast to a char*, allowing byte order to be reversed.
	// Note that by passing a reference to this function, the original variable itself is modified and is no longer
	// useful as an double for any numerical or comparison operation
	char swap;

	for(int i=0;i<4;i++) 
	{
		swap = double_array[7-i];
		double_array[7-i]=double_array[i];
		double_array[i]=swap;
	}
}

//-----------------------------------------------------------------------------
void PDCFileWriter::writeAttribute(ofstream &outfile, const int attrib_length, const char* attrib_name, const int attrib_type)
{
	// writes the text for each new attribute to the data file
	// 'attrib_length' is the number of characters in the attribute name 'attrib_name'
	// attrib type is the integer data type, following the convention:
	// 0=int, 1=intArray, 2=intArray, 3=double, 4= doublearray, 5=vector(3doubles), 6=vectorArray(array of 3doubles)

	int attribLen = attrib_length;
	swapInt((char*) &attribLen);
	outfile.write((char*) &attribLen, sizeof(int));

	for(int i=0;i<attrib_length;i++)
	{
		outfile.write((char*) &attrib_name[i], sizeof(char));
	}

	int attribType=attrib_type;
	swapInt((char*) &attribType);
	outfile.write((char*) &attribType, sizeof(int));
}

//-----------------------------------------------------------------------------
void PDCFileWriter::write3vector(ofstream &outfile, const double x1, const double y1, const double z1)
{
	//writes a 3-column data vector to the file: use this for data type 5
	double x = x1;
	double y = y1;
	double z = z1;

	swapDouble((char*) &x);
	swapDouble((char*) &y);
	swapDouble((char*) &z);

	outfile.write((char*) &x, sizeof(double));
	outfile.write((char*) &y, sizeof(double));
	outfile.write((char*) &z, sizeof(double));
}
//-----------------------------------------------------------------------------
void PDCFileWriter::write1scalar(ofstream &outfile, const double x1)
{
	//writes a 1-column data scalar/number to the file: use this for data type 3
	double x = x1;
	swapDouble((char*) &x);
	outfile.write((char*) &x,sizeof(double));
}

//-----------------------------------------------------------------------------
void PDCFileWriter::write1intscalar(ofstream &outfile, const int ix)
{
	//writes a 1-column data scalar/number to the file: use this for data type 0
	int x = ix;
	swapInt((char*) &x);
	outfile.write((char*) &x,sizeof(int));
}

//-----------------------------------------------------------------------------