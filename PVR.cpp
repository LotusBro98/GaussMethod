#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string.h>

class Matrix
{
	public:

    explicit Matrix(int n, int m)
	{
		this->_n = n;
        this->_m = m;
		this->_a = new double[n * m];
	}

	//Read matrix from stream
    explicit Matrix(std::istream &stream)
    {
        int n;
		int m;
		stream >> n;
		m = n;

        this->_n = n;
        this->_m = m;
        this->_a = new double[n * m];

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                stream >> a(i, j);
            }
        }
    }

	//Generate matrix with aij = f(i, j)
	Matrix(double (*f)(int i, int j), int n, int m) : Matrix(n, m)
	{
		for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
				a(i, j) = f(i, j);
	}

    friend std::ostream& operator<< (std::ostream & os, Matrix * A)
    {
        for (int i = 0; i < A->_n; ++i) {
            for (int j = 0; j < A->_m; ++j) {
                os << std::setw(3) << A->a(i, j) << ' ';
            }
            os << std::endl;
        }
        os << std::endl;
        return os;
    }

    void divLine(int i, double c)
    {
        for (int j = 0; j < _m; ++j) {
            a(i, j) /= c;
        }
    }

    friend Matrix * operator * (Matrix & A, Matrix & B)
    {
        if (A._m != B._n)
            return nullptr;

        Matrix * C = new Matrix(A._n, B._m);
        for (int i = 0; i < C->_n; ++i) {
            for (int j = 0; j < C->_m; ++j) {
                C->a(i, j) = 0;
                for (int k = 0; k < A._m; ++k) {
                    C->a(i, j) += A.a(i, k) * B.a(k, j);
                }
            }
        }

        return C;
    }

	Matrix * operator -= (Matrix & A)
	{
		if (A._n != _n || A._m != _m)
            return nullptr;

        for (int i = 0; i < _n; ++i)
            for (int j = 0; j < _m; ++j)
                a(i, j) -= A.a(i, j);

		return this;
	}

    double vectorNormSqrt()
    {
        double len2 = 0;
        if (_n == 1)
            for (int j = 0; j < _m; ++j)
                len2 += a(0, j) * a(0, j);
        else if (_m == 1)
            for (int i = 0; i < _n; ++i)
                len2 += a(i, 0) * a(i, 0);
        else
            return -1;

        return std::sqrt(len2);
    }

    Matrix * copy()
    {
        Matrix * A = new Matrix(_n, _m);

        std::copy(_a, _a + _n * _m,  A->_a);

        return A;
    }

	//Multiply line i of this matrix and column j of v
	double mulLines(int i, Matrix * v, int j = 0)
	{
		if (v->_n != _m)
			throw;

		double sum = 0;
		for (int k = 0; k < _m; k++)
			sum += a(i, k) * v->a(k, j);
		
		return sum;
	}

	//Multiply and store result in this matrix
	void setMul(Matrix * A, Matrix * B)
	{
	   	if (A->_m != B->_n || _n != A->_n || _m != B->_m)
			throw;

        for (int i = 0; i < _n; ++i) {
            for (int j = 0; j < _m; ++j) {
                a(i, j) = 0;
                for (int k = 0; k < A->_m; ++k) {
                    a(i, j) += A->a(i, k) * B->a(k, j);
                }
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

Matrix * solveIter(Matrix * A, Matrix * b, double w = 1, double eps = 1e-4)
{
	if (A->n() != b->n())
		throw;
	
	Matrix * B = A->copy();
	Matrix * c = b->copy();
	for (int i = 0; i < A->n(); i++)
	{
		B->divLine(i, A->a(i, i));
		c->a(i, 0) /= A->a(i, i);
	}
	//Bij = Aij/Aii
	//ci = bi/Aii

	Matrix * x = new Matrix(A->m(), 1);	
	Matrix * r = new Matrix(A->m(), 1);
	do
	{
		//B = L + E + U
		//x(k+1) = x(k) - w*(L * x(k+1) + (E + U) * x(k) - c)
		//...
		//xi := xi(k+1)
		//...
		for (int i = 0; i < x->n(); i++)
			x->a(i, 0) = x->a(i, 0) - w * ( B->mulLines(i, x) - c->a(i, 0) );

		r->setMul(A, x);
		*r -= *b;
		std::cerr << "|r| = " << r->vectorNormSqrt() << "\n";
	}
	while (r->vectorNormSqrt() > eps);

	return x; 
}

// ------------------------------------------------------------

double stolbGen(int i, int j)
{
	return i + 1;
}

int main()
{
	std::ifstream file;

	file.open("matrix_3.txt");
	Matrix * A = new Matrix(file);
	file.close();

	Matrix * xPrimal = new Matrix(stolbGen, A->m(), 1);
	Matrix * b = (*A) * (*xPrimal);

	Matrix * xIter = solveIter(A, b, 1.5, 1e-4);
	std::cout << "\n" << xIter;
}
