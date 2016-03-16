#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "lib.h"
#include "fileio.h"

using namespace std;

DataFile::DataFile(int num_of_files, string* input_file_names, int cols, int SkipRows)
{
    data_storage = new float**[num_of_files];
    Rows = new int[num_of_files];

    for(int i=0; i<num_of_files; i++)
    {
        data_storage[i] = Read(input_file_names[i], cols, SkipRows, i);
    }
}

float** DataFile::Read(string input_file_names, int cols, int SkipRows, int index)
{
    fstream input_file;
    input_file.open(input_file_names.c_str(), ios::in);
    if(!input_file.is_open())  {error(1);}

    Rows[index] = NumOfRows(input_file) - SkipRows;
    Cols = cols;

	float** data = (float**)matrix(Rows[index], Cols, sizeof(float));

	string comments;
	//read comments lines
	for(int i=0; i<SkipRows; i++)
		getline(input_file, comments);

	for(int i=0; i < Rows[index]; i++)
	{
		for (int j=0; j< Cols; j++)
			input_file >> data[i][j];
	}
	input_file.close();

	return data;
}

int DataFile::NumOfRows(fstream& file)
{
  int length;
  char * buffer;

  // get length of file:
  file.seekg (0, ios::end);
  length = file.tellg();
  file.seekg (0, ios::beg);

  // allocate memory:
  buffer = new char [length];

  // read data as a block:
  file.read (buffer,length);

  //count instances of endlines
  int end_lines = 0;
  for(int i=0; i<length; i++)
  {
      if(buffer[i] == '\n')
        end_lines++;
  }

  delete[] buffer;
  file.seekg(0, ios::beg);
  return end_lines;
}

void DataFile::error(int TypeOfError)
{
	exit(EXIT_FAILURE);
}
