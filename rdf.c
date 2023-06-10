#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

typedef struct trajectory
{
	int atomID, atomType, molType, ix, iy, iz;
	float x, y, z;
	int isEndGroup;
} TRAJECTORY;

typedef struct vector
{
	float x1, y1, z1;
	float x2, y2, z2;
	float xc, yc, zc;
} VECTOR;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct rdf
{
	float rlo, rhi, gofr;
} RDF;

typedef struct stats
{
	float average, standardDeviation;
} STATS;

typedef struct orderParameterBins
{
	float orderParameter, rlo, rhi, count;
} ORDERPARAMETER_BINS;

int countNAtoms (FILE *file_dump, int *nAtomEntries)
{
	int nAtoms, currentAtomID, nAtomsFixed;
	char lineString[2000];
	rewind (file_dump);

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, file_dump); }

	sscanf (lineString, "%d\n", &nAtoms);
	(*nAtomEntries) = nAtoms;
	rewind (file_dump);
	nAtomsFixed = nAtoms;

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > nAtoms) {
			nAtomsFixed = currentAtomID; }
	}

	return nAtomsFixed;
}

SIMULATION_BOUNDARY readDumpBoundary (FILE *file_dump, SIMULATION_BOUNDARY boundary)
{
	rewind (file_dump);
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}

	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.xlo, &boundary.xhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.ylo, &boundary.yhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.zlo, &boundary.zhi);
	rewind (file_dump);

	boundary.xLength = boundary.xhi - boundary.xlo;
	boundary.yLength = boundary.yhi - boundary.ylo;
	boundary.zLength = boundary.zhi - boundary.zlo;

	return boundary;
}

TRAJECTORY *initializeAtoms (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].molType = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].isEndGroup = 0;
	}

	return atoms;
}

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtomEntries, SIMULATION_BOUNDARY *boundary)
{
	char lineString[2000];
	int currentAtomID = 1;

	for (int i = 0; i < 5; ++i) {
		fgets (lineString, 2000, file_dump); }

	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).xlo, &(*boundary).xhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).ylo, &(*boundary).yhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).zlo, &(*boundary).zhi);
	fgets (lineString, 2000, file_dump);

	for (int i = 0; i < nAtomEntries; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);
		sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
		atoms[currentAtomID - 1].isEndGroup = 0;
	}

	return atoms;
}

int countNAtoms_byType (TRAJECTORY *atoms, int nAtoms, int nAtoms_byType, int requiredAtomType)
{
	nAtoms_byType = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == requiredAtomType) {
			nAtoms_byType++; }
	}

	return nAtoms_byType;
}

RDF *initializeRDF (RDF *rdf_atomType2, int rdf_nBins, float rdf_binWidth)
{
	for (int i = 0; i < rdf_nBins; ++i)
	{
		if (i == 0) {
			rdf_atomType2[0].rlo = 0;
			rdf_atomType2[0].rhi = rdf_binWidth; }
		else {
			rdf_atomType2[i].rlo = rdf_atomType2[i - 1].rhi;
			rdf_atomType2[i].rhi = rdf_atomType2[i].rlo + rdf_binWidth; }

		rdf_atomType2[i].gofr = 0;
	}

	return rdf_atomType2;
}

int main(int argc, char const *argv[])
{
	if (argc != 6)
	{
		printf("REQUIRED ARGUMENTS\n~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = program\n {~} argv[1] = (char *)input dump filename\n {~} argv[2] = (int)atom type 1\n {~} argv[3] = (int)atom type 2\n {~} argv[4] = (float)maximum rdf distance to compute\n {~} argv[5] = (float)bin width\n\n");
		exit (1);
	}

	FILE *file_dump;
	file_dump = fopen (argv[1], "r");

	int nAtomEntries, nAtoms = countNAtoms (file_dump, &nAtomEntries), atomType1 = atoi (argv[2]), atomType2 = atoi (argv[3]), file_status;
	SIMULATION_BOUNDARY boundary;
	boundary = readDumpBoundary (file_dump, boundary);

	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	printf("Number of atoms in the trajectory file: %d\n", nAtoms);

	atoms = initializeAtoms (atoms, nAtoms);

	rewind (file_dump);
	file_status = fgetc (file_dump);

	float rdf_maxdist = atof (argv[4]), rdf_binWidth = atof (argv[5]);
	int rdf_nBins = ceil (rdf_maxdist / rdf_binWidth);
	RDF *rdf_atomType2;

	int currentTimestep = 0, nAtoms_atomType2;
	float density_atomType2;

	while (file_status != EOF)
	{
		fprintf(stdout, "Scanning timestep: %d                           \r", currentTimestep);
		fflush (stdout);

		atoms = readTimestep (file_dump, atoms, nAtomEntries, &boundary);

		if (currentTimestep == 0)
		{
			nAtoms_atomType2 = countNAtoms_byType (atoms, nAtoms, nAtoms_atomType2, atomType2);
			rdf_atomType2 = (RDF *) malloc (rdf_nBins * sizeof (RDF));
			density_atomType2 = nAtoms_atomType2 / (boundary.xhi - boundary.xlo) * (boundary.yhi - boundary.ylo) * (boundary.zhi - boundary.zlo);
			rdf_atomType2 = initializeRDF (rdf_atomType2, rdf_nBins, rdf_binWidth);
		}

		currentTimestep++;
		file_status = fgetc (file_dump);
	}


	fclose (file_dump);
	return 0;
}