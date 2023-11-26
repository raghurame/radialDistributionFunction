#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

int getNLines (FILE *inputRDF, int nLines)
{
	nLines = 0;
	char lineString[3000];

	while (fgets (lineString, 2000, inputRDF) != NULL)
	{
		nLines++;
	}

	rewind (inputRDF);

	return nLines;
}

void getValues (float **rdfValues, float **distance, int nLines, FILE *inputRDF)
{
	char lineString[3000];
	int currentIndex = 0;

	while (fgets (lineString, 3000, inputRDF) != NULL)
	{
		sscanf (lineString, "%f %*f %f\n", &(*distance)[currentIndex], &(*rdfValues)[currentIndex]);
		currentIndex++;
	}
}

void getFirstMinima (float *rdfValues, float *distance, int nLines, int *firstPeak_start, int *firstPeak_end)
{
	float previousRDF = 0, currentRDF = 0;
	(*firstPeak_start) = 0; (*firstPeak_end) = 0;

	for (int i = 0; i < nLines; ++i)
	{

		currentRDF = rdfValues[i];

		if (currentRDF < previousRDF)
		{
			(*firstPeak_start) = i;
			break;
		}

		previousRDF = rdfValues[i];
	}

	currentRDF = 0; previousRDF = 99999;

	for (int i = (*firstPeak_start); i < nLines; ++i)
	{
		currentRDF = rdfValues[i];

		if (currentRDF > previousRDF)
		{
			(*firstPeak_end) = i;
			break;
		}

		previousRDF = rdfValues[i];
	}
}

int main(int argc, char const *argv[])
{
	FILE *inputRDF;
	inputRDF = fopen (argv[1], "r");

	int nLines = getNLines (inputRDF, nLines), firstPeak_start, firstPeak_end;
	float *rdfValues, *distance;
	rdfValues = (float *) malloc (nLines * sizeof (float));
	distance = (float *) malloc (nLines * sizeof (float));

	getValues (&rdfValues, &distance, nLines, inputRDF);
	getFirstMinima (rdfValues, distance, nLines, &firstPeak_start, &firstPeak_end);

	printf("First RDF peak:\n  start: %d\n  end: %d\n\n", firstPeak_start, firstPeak_end);

	fclose (inputRDF);
	return 0;
}