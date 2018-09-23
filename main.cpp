#include <iostream>
#include <cstdlib>
#include <exception>
#include <fstream>

class Matrix
{
	public:

    explicit Matrix(int n)
	{
		this->_n = n;
		this->_a = new double[n];
	}

    explicit Matrix(std::istream &stream)
    {
        int n;
        stream >> n;

        this->_n = n;
        this->_a = new double[n];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                stream >> a(i, j);
            }
        }

    }

    friend std::ostream& operator<< (std::ostream & os, Matrix & A)
    {
        for (int i = 0; i < A._n; ++i) {
            for (int j = 0; j < A._n; ++j) {
                os << A.a(i, j);
            }
        }
        os << std::endl;
    }

    ~Matrix()
    {
        delete(_a);
    }

	int n()
    {
        return _n;
    }

	double & a(int i, int j)
    {
        return _a[_n * i + j];
    }

	private:
	int _n;
	double * _a;
};

int main(int argc, char * argv[])
{
	if (argc != 2) {
        std::cerr << "Usage: ./gauss <matrix_filename>\n";
        return 1;
    }

	char* filename = argv[1];

	std::ifstream file;
    file.open(filename);
	
	Matrix mat = Matrix(file);

	std::cout << mat;

    file.close();
	return 0;
}


// ------------------------------------------------------------
