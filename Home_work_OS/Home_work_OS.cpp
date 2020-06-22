#include <iostream>
#include <omp.h>
#include <random>
#include <Windows.h>
#include <time.h>

void Task1()
{
#pragma omp parallel num_threads(8)
	{
		std::cout << "Hello World! " << " Thread id: " << omp_get_thread_num();
		std::cout << " Number of threads: " << omp_get_num_threads() << std::endl;
	}
}

void Task2()
{
	omp_set_num_threads(3);
#pragma omp parallel if(omp_get_max_threads() > 1)
	{
		std::cout << "Id: " << omp_get_thread_num();
		std::cout << " Threads number: " << omp_get_num_threads() << std::endl;
	}


	omp_set_num_threads(3);
#pragma omp parallel if(omp_get_max_threads() > 1)
	{
		std::cout << "Id: " << omp_get_thread_num();
		std::cout << " Threads number: " << omp_get_num_threads() << std::endl;
	}
}


void Task3()
{
	int a = 100, b = 100;
	std::cout << "First Test!" << std::endl;
	std::cout << "Before a: " << a << " , b: " << b << std::endl;
	omp_set_num_threads(2);
#pragma omp parallel firstprivate (b) private (a) 
	{
		a = 0;
		a += omp_get_thread_num();
		b += omp_get_thread_num();
		std::cout << "In a: " << a << " , b: " << b << std::endl;
	}
	std::cout << "After a: " << a << " , b: " << b << std::endl;

	std::cout << "Second Test!" << std::endl;
	std::cout << "Before a: " << a << " , b: " << b << std::endl;
	omp_set_num_threads(4);
#pragma omp parallel shared (a) private (b) 
	{
		b = 0;
		a -= omp_get_thread_num();
		b -= omp_get_thread_num();
		std::cout << "In a: " << a << " , b: " << b << std::endl;
	}
	std::cout << "After a: " << a << " , b: " << b << std::endl;
}


void fill(int* mas, int size)
{
	for (int i = 0; i < size; ++i)
		mas[i] = rand() % size;
}

void show(int* mas, int size)
{
	for (int i = 0; i < size; ++i)
		std::cout << mas[i] << " ";
	std::cout << std::endl;
}


void Task4()
{
	const int size = 10;
	int a[size], b[size];
	fill(a, size);
	fill(b, size);
	int max = -1, min = size+1;
#pragma omp parallel shared(a,b) num_threads(2)
	{
		
		if (omp_get_thread_num())
		{
			for (int i = 0; i < size; ++i)
			{
				if (b[i] > max) max = b[i];
			}
		}
		else
		{
			for (int i = 0; i < size; ++i)
			{
				if (a[i] < min) min = a[i];
			}
		}
	
	}
	show(a, size);
	std::cout << "min: " << min << std::endl;
	show(b, size);
	std::cout << "max: " << max << std::endl;
}

void Task5()
{
	const int n = 6, m = 8;
	int d[n][m];
	for (int i = 0; i < n; ++i) fill(d[i], m);
	double a=0, b=0, c= 0;
#pragma omp parallel shared(a,b,c)
	{
#pragma omp sections firstprivate(a,b,c)
		{
#pragma omp section
			{
				std::cout << "First section, id: " << omp_get_thread_num() << std::endl;
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < m; ++j)
						a += d[i][j];
				std::cout << "Result: " << a / (double)(n * m) << std::endl;
			}
#pragma omp section
			{
				a = d[0][0]; // min
				b  = d[0][0]; // max
				std::cout << "Second section, id: " << omp_get_thread_num() << std::endl;
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < m; ++j)
					{
						if (d[i][j] > b) b = d[i][j];
						if (d[i][j] < a) a = d[i][j];
					}
				std::cout << "Result, min " << a <<" max " << b << std::endl;
			}
#pragma omp section
			{
				std::cout << "Third section, id: " << omp_get_thread_num() << std::endl;
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < m; ++j)
						if (d[i][j] % 3 == 0) a += 1;
				std::cout << "Result: " << a  << std::endl;
			}
		}
	}

	for (int i = 0; i < n; ++i) show(d[i], m);
	std::cout << "End" << std::endl;
}

void Task6()
{
	const int n = 10;
	int a[n], b[n];
	fill(a, n);
	fill(b, n);
	double mean_a = 0, mean_b = 0;
	double temp = 0;

	omp_set_dynamic(1);
#pragma omp parallel for reduction(+ : temp)
	for (int i = 0; i < n; ++i)
		temp += a[i];
	mean_a = temp / n;

	temp = 0;
#pragma omp parallel for reduction(+ : temp)
	for (int i = 0; i < n; ++i)
		temp += b[i];
	mean_b = temp / n;


	std::cout << "Result, mean A:  " << mean_a << " mean B: " << mean_b << std::endl;
	show(a,n);
	show(b,n);
}


void Task7()
{
	const int n = 12;
	int a[n], b[n], c[n];
	int i;
	std::cout << "First parallel zone!\n\n";
	omp_set_num_threads(3);
#pragma omp parallel for schedule(static, 2)
	for (i = 0; i < n; ++i)
	{
		std::cout << "Num_threads: " << omp_get_num_threads() << " Thread: " << omp_get_thread_num();
		std::cout << " Iteration: " << i << std::endl;
		a[i] = i;
		b[i] = i;
	}
	std::cout  << "\nSecond parallel zone!\n\n";
	omp_set_num_threads(4);
#pragma omp parallel for schedule(dynamic, 2)
	for (i = 0; i < n; ++i)
	{
		std::cout << "Num_threads: " << omp_get_num_threads() << " Thread: " << omp_get_thread_num();
		std::cout << " Iteration: " << i << std::endl;
		c[i] = a[i] + b[i];
	}

	std::cout  << "\nResult\n\n";
	show(a, n);
	show(b, n);
	show(c, n);
}

bool check(int* v1, int* v2, int n)
{
	for (int i = 0; i < n; ++i)
		if (v1[i] != v2[i]) return false;
	return true;
}

bool check(int** m1, int** m2, int n, int m)
{
	for (int i = 0; i < n; ++i)
		for(int j = 0; j < m; ++j)
			if (m1[i][j] != m2[i][j]) return false;
	return true;
}


int* parallel_mult(int** matr, int* v, int n, int m)
{
	int* result = new int[n];
	int temp = 0;
	std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
	double start = omp_get_wtime();
#pragma omp parallel for private(temp)
	for (int i = 0; i < n; ++i)
	{
		temp = 0;
		for (int j = 0; j < m; ++j)
			temp += v[j] * matr[i][j];
#pragma omp critical
		{
			result[i] = temp;
		}
	}
	double end = omp_get_wtime();
	std::cout << "Parallel time: " << (end - start) << std::endl;
	return result;
}

int* standart_mult(int** matr, int* v, int n, int m)
{
	int* result = new int[n];
	int temp;
	double start = omp_get_wtime();
	for (int i = 0; i < n; ++i)
	{
		temp = 0;
		for (int j = 0; j < m; ++j)
			temp += v[j] * matr[i][j];

		result[i] = temp;
	}
	double end = omp_get_wtime();
	std::cout << "Oneway time: " << end - start << std::endl;
	return result;
}

void Task8()
{
	int n = 10000, m = 10000;
	int* v = new int[m];
	int** matr = new int*[n];
	for (int i = 0; i < n; ++i)
		matr[i] = new int[m];
	fill(v, m);
	for (int i = 0; i < n; ++i) fill(matr[i], m);

	omp_set_num_threads(3);
	int* result1 = parallel_mult(matr, v, n, m);
	int* result2 = standart_mult(matr, v, n, m);
	std::cout << "Check: " << check(result1, result2, n) << std::endl;
}

void Task9()
{
	int n = 6, m = 8;
	int** d = new int* [n];
	for (int i = 0; i < n; ++i) d[i] = new int[m];

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			d[i][j] = rand() % (n * m);

	int max = d[0][0], min = d[0][0];

#pragma omp parallel for
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
		{
#pragma omp critical
			{
				if (d[i][j] > max) max = d[i][j];
				else if (d[i][j] < min) min = d[i][j];
			}
		}


	std::cout << "Max: " << max << " Min: " << min << std::endl;
	for (int i = 0; i < n; ++i) show(d[i], m);

	for (int i = 0; i < n; ++i) delete[] d[i];
	delete[] d;
}

void Task10()
{
	int n = 30, counter = 0;
	int* a = new int[n];
	fill(a, n);

#pragma parallel for
	for (int i = 0; i < n; ++i)
	{
		if (a[i] % 9 == 0)
		{
#pragma omp atomic
			++counter;
		}
	}

	std::cout << "Counter: " << counter << std::endl;
	delete[] a;
}

void Task11()
{
	int n = 100000, max = -1;
	bool flag = false;
	int* a = new int[n];
	fill(a, n);

#pragma parallel for
	for (int i = 0; i < n; ++i)
	{
		if (a[i] % 7 == 0)
		{
#pragma omp critical
			{
				if (!flag || a[i] > max) max = a[i];
			}
		}
	}

	std::cout << "Max: " << max << std::endl;
	delete[] a;
}

int main()
{
	return 0;
}

