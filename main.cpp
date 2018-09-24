#include <iostream>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>

class Matrix
{
	public:

    explicit Matrix(int n, int m)
	{
		this->_n = n;
        this->_m = m;
		this->_a = new double[n * m];
	}

    explicit Matrix(std::istream &stream)
    {
        int n;
        stream >> n;

        this->_n = n;
        this->_m = n;
        this->_a = new double[n * n];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                stream >> a(i, j);
            }
        }
    }

    friend std::ostream& operator<< (std::ostream & os, Matrix * A)
    {
        for (int i = 0; i < A->_n; ++i) {
            for (int j = 0; j < A->_m; ++j) {
                os << std::setw(3) << A->a(i, j) << ' ';
            }
            os << std::endl;
        }

        return os;
    }

    Matrix * operator | (Matrix & b)
    {
        int n = std::max(_n, b._n);
        int m = _m + b._m;

        Matrix * compound = new Matrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (i < _n && j < _m)
                    compound->a(i, j) = a(i, j);
                else if (i < b._n && j >= _m)
                    compound->a(i, j) = b.a(i, j - _m);
                else
                    compound->a(i, j) = 0;
            }
        }

        return compound;
    }

    void addLine(int i0, int i1, double c = 1)
    {
        for (int j = 0; j < _m; ++j) {
            a(i0, j) += a(i1, j) * c;
        }
    }

    void subtractLine(int i0, int i1, double c = 1)
    {
        for (int j = 0; j < _m; ++j) {
            a(i0, j) -= a(i1, j) * c;
        }
    }

    void mulLine(int i, double c)
    {
        for (int j = 0; j < _m; ++j) {
            a(i, j) *= c;
        }
    }

    void divLine(int i, double c)
    {
        for (int j = 0; j < _m; ++j) {
            a(i, j) /= c;
        }
    }

    void swapLines(int i1, int i2)
    {
        for (int j = 0; j < _m; ++j) {
            double t = a(i1, j);
            a(i1, j) = a(i2, j);
            a(i2, j) = t;
        }
    }

    void gauss()
    {
        for (int j = 0; j < _m; ++j)
        {
            int in0;
            for (in0 = j; in0 < _n; in0++)
                if (a(in0, j) != 0)
                    break;
            if (in0 == _n)
                continue;

            if (in0 != j)
                swapLines(in0, j);

            divLine(j, a(j, j));

            for (int i = 0; i < _n; ++i) {
                if (i == j || a(i, j) == 0)
                    continue;

                subtractLine(i, j, a(i, j));
            }
        }
    }

    ~Matrix()
    {
        delete(_a);
    }

	int n()
    {
        return _n;
    }

    int m() const {
        return _m;
    }

    inline double & a(int i, int j)
    {
        return _a[_m * i + j];
    }

	private:
	int _n;
    int _m;
	double * _a;
};

int main(int argc, char * argv[])
{
	if (argc != 2) {
        std::cerr << "Usage: ./gauss <matrix_filename>" << std::endl;
        return 1;
    }

	char* filename = argv[1];

	std::ifstream file(filename);
    if (!file.good())
    {
        std::cerr << "Wrong file: " << filename << std::endl;
        return 1;
    }
	
	Matrix * mat = new Matrix(file);
    Matrix * vec = new Matrix(mat->n(), 1);
    for (int i = 0; i < mat->n(); ++i) {
        vec->a(i, 0) = i + 1;
    }

    Matrix * expanded = *mat | *vec;

	std::cout << expanded << "\n";

    expanded->gauss();

    std::cout << expanded;

    //file.close();
	return 0;
}


// ------------------------------------------------------------
