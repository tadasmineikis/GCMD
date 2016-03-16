#include <fstream>

using namespace std;

class DataFile
{
	public:
		DataFile(int num_of_files, string* input_file_names, int cols, int SkipRows);
		~DataFile();

		float** Read(string input_file_names, int cols, int SkipRows, int index);
		int NumOfRows(fstream& file);
		void error(int);

		float ***data_storage;

		int *Rows, Cols;
};
