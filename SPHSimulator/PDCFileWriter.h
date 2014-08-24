#ifndef PDC_FILE_WRITER_H
#define PDC_FILE_WRITER_H


#include <string>

using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>

#include "Particle.h"

/// \file PDCFileWriter.h
/// \brief Contains the PDCFileWriter class
/// \author Justina Reingardtaite
/// \version 1.0
/// \date 29/03/2012

/// \class PDCFileWriter
/// \brief Writes particle data into a .pdc file, which is a standard format of Maya cache file.

/// The following section is from :-
/// Amit Chourasia (San Diego Supercomputer Center) and Peter J. Lu (Harward University). Maya Particle Disk Cache (PDC) Utils [online]. [Accessed 2012]. 
/// Available from: http://users.sdsc.edu/~amit/web/book/export/html/21

class PDCFileWriter
{
public:
	
	// \brief Ctor
    /// @param[in] fileName the file to be used to name cache files    
	PDCFileWriter( char* fileName);

	// \brief Exports particle data into .pdc file
    /// @param[in] particle the reference to the particle vector to be written into the file
	/// @param[in] frameNumber the frame numer to be used in creating the cache file name
	void writePDCFile( std::vector<Particle*>& particle, int frameNumber);

	// \brief Writes the .pdc file header
    /// @param[in] outfile the reference to the opened file where the particle data will be written
	void writePDCHeader(ofstream &outfile);

	/// \brief Dtor
	~PDCFileWriter();

	/// \brief Gets the first simbol of the file format
	char getFormat0() const;

	// \brief Gets the second simbol of the file format
	char getFormat1() const;	

	// \brief Gets the third simbol of the file format
	char getFormat2() const;	

	// \brief Gets the forth simbol of the file format
	char getFormat3() const;

	// \brief Gets the format version
	int getFormatVersion() const;
	
	// \brief Gets the byte order
	int getByteOrder() const;

	// \brief Gets the extra symbol 1
	int getExtra1() const;	

	// \brief Gets the extra symbol 2
	int getExtra2() const;
	
	// \brief Gets the particle count
	int getParticleCount() const;	

	// \brief Gets the attribute count
	int getAttributeCount() const;

	// \brief Swaps byte order of 4-element character array.
	void swapInt(char* int_array);

	// \brief Swaps byte order of 8-element character array.
	void swapDouble(char* double_array);

	// \brief writes the text for each new attribute to the data file
	void writeAttribute(ofstream &outfile, const int attrib_length, const char* attrib_name, const int attrib_type);

	// \brief Writes a 3-column data vector to the file: use this for data type 5
	void write3vector(ofstream &outfile, const double x1, const double y1, const double z1);
	
	// \brief Writes a 1-column data scalar/number to the file: use this for data type 3
	void write1scalar(ofstream &outfile, const double x1);

	// \brief Writes a 1-column data scalar/number to the file: use this for data type 0
	void write1intscalar(ofstream &outfile, const int x1);

private:

	// \brief The name of output file
	char * m_outputFile;

	// \brief Holds four format symbols
	char m_format[4];

	// \brief The format format
	int m_formatVersion;

	// \brief The byte order
	int m_byteOrder;

	// \brief The extra symbol 1
	int m_extra1;

	// \brief The extra symbol 2
	int m_extra2;

	// \brief The particle count
	int m_particleCount;

	// \brief The attribute count
	int m_attributeCount;
};

#endif
