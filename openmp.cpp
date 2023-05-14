#include<iostream>
#include<omp.h>
#include <windows.h>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX¡¢AVX2
using namespace std;
const int N=2500;
const int p=1;
float** m;
//float** n;
LARGE_INTEGER freq, t1, t2;
const int NUM_THREADS = 7;

void m_reset()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
            m[i][j] = 0;
        m[i][i] = 1.0;
        for (int j = i + 1; j < N; j++)
            m[i][j] = rand();
    }
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
            for (int j = 0; j < N; j++)
                m[i][j] = m[k][j] + m[i][j];
}

void normal()//
{
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
    for (int k = 0; k < N; k++)
        {
                for (int j = k + 1; j < N; j++)
                {
                        m[k][j] = m[k][j] / m[k][k];

                }
                m[k][k] = 1.0;
                for (int i = k + 1; i < N; i++)
                {
                        for (int j = k + 1; j < N; j++)
                        {
                                m[i][j] -= m[i][k] * m[k][j];
                        }
                        m[i][k] = 0;
                }

        }
    QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "normal_time: " << sumtime/p << "ms" << endl;
}

void openmp()
{
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
    #pragma omp parallel //if(parallel),num_threads(NUM_THREADS),private(i,j,k,tmp)
    for (int k = 0; k < N; k++)
        {


            float tmp = m[k][k];
            for (int j = k + 1; j < N; j++)
                {
                        m[k][j] = m[k][j] / tmp;

                }
                m[k][k] = 1.0;
                #pragma omp for
                for (int i = k + 1; i < N; i++)
                {
                        float tmp= m[i][k];
                        for (int j = k + 1; j < N; j++)
                        {
                                m[i][j] -= tmp * m[k][j];
                        }
                        m[i][k] = 0;
                }

        }

    QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     cout << "openmp_time: " << sumtime/p<< "ms" << endl;
     }

}

void openmp_dynamic()
{
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
     #pragma omp parallel //if(parallel),num_threads(NUM_THREADS),private(i,j,k,tmp)
    for (int k = 0; k < N; k++)
        {



            float tmp=m[k][k];
            for (int j = k + 1; j < N; j++)
            {
                    m[k][j] = m[k][j] / tmp;

            }
            m[k][k] = 1.0;
            #pragma omp for schedule(dynamic,1)
            for (int i = k + 1; i < N; i++)
            {


                float tmp= m[i][k];
                for (int j = k + 1; j < N; j++)
                {
                    m[i][j] -= tmp * m[k][j];
                }
                m[i][k] = 0;
            }

        }
    QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "openmp_dynamic_time: " << sumtime/p << "ms" << endl;
}


void openmp_simd()
{
    __m128 t5,t6,t3,t4;
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x=0;x<p;x++)
    {
        m_reset();
    QueryPerformanceCounter(&t1);
    #pragma omp parallel num_threads(NUM_THREADS)
	for (int k = 0; k < N; k++)
        {
        t5 = _mm_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 4 < N; j += 4)
		{
			t6 = _mm_loadu_ps(&m[k][j]);
			t6 = _mm_div_ps(t6, t5);
			_mm_storeu_ps(&m[k][j], t6);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;
        #pragma omp for
		for (int i = k + 1; i < N; i++)
		{
			t5 = _mm_set1_ps(m[i][k]);
            int j = k + 1;
			for (; j + 4 < N; j += 4)
			{
				t6 = _mm_loadu_ps(&m[k][j]);
				t3 = _mm_loadu_ps(&m[i][j]);
				t4 = _mm_mul_ps(t5,t6);
				t3 = _mm_sub_ps(t3, t4);
				_mm_storeu_ps(&m[i][j], t3);
			}
			for (; j < N; j++)
				m[i][j] -= m[i][k] * m[k][j];
			m[i][k] = 0;
		}
	}
	QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
    }
    cout << "openmp_simd_time: " << sumtime/p << "ms" << endl;
}

void openmp_simd_dy()
{
    __m128 t5,t6,t3,t4;
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x=0;x<p;x++)
    {
        m_reset();
    QueryPerformanceCounter(&t1);
    #pragma omp parallel num_threads(NUM_THREADS)
	for (int k = 0; k < N; k++)
        {
        t5 = _mm_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 4 < N; j += 4)
		{
			t6 = _mm_loadu_ps(&m[k][j]);
			t6 = _mm_div_ps(t6, t5);
			_mm_storeu_ps(&m[k][j], t6);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;
        #pragma omp for schedule(dynamic,5)
		for (int i = k + 1; i < N; i++)
		{
			t5 = _mm_set1_ps(m[i][k]);
            int j = k + 1;
			for (; j + 4 < N; j += 4)
			{
				t6 = _mm_loadu_ps(&m[k][j]);
				t3 = _mm_loadu_ps(&m[i][j]);
				t4 = _mm_mul_ps(t5,t6);
				t3 = _mm_sub_ps(t3, t4);
				_mm_storeu_ps(&m[i][j], t3);
			}
			for (; j < N; j++)
				m[i][j] -= m[i][k] * m[k][j];
			m[i][k] = 0;
		}
	}
	QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
    }
    cout << "openmp_simd_dy_time: " << sumtime/p << "ms" << endl;
}

void openmp_avx()
{
    __m256 t5,t6,t3,t4;
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x=0;x<p;x++)
    {
        m_reset();
    QueryPerformanceCounter(&t1);
    #pragma omp parallel num_threads(NUM_THREADS)
	for (int k = 0; k < N; k++)
        {
        t5 = _mm256_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 8 < N; j += 8)
		{
			t6 = _mm256_loadu_ps(&m[k][j]);
			t6 = _mm256_div_ps(t6, t5);
			_mm256_storeu_ps(&m[k][j], t6);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;
        #pragma omp for
		for (int i = k + 1; i < N; i++)
		{
			t5 = _mm256_set1_ps(m[i][k]);
            int j = k + 1;
			for (; j + 8 < N; j += 8)
			{
				t6 = _mm256_loadu_ps(&m[k][j]);
				t3 = _mm256_loadu_ps(&m[i][j]);
				t4 = _mm256_mul_ps(t5,t6);
				t3 = _mm256_sub_ps(t3, t4);
				_mm256_storeu_ps(&m[i][j], t3);
			}
			for (; j < N; j++)
				m[i][j] -= m[i][k] * m[k][j];
			m[i][k] = 0;
		}
	}
	QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
    }
    cout << "openmp_avx_time: " << sumtime/p << "ms" << endl;
}

void openmp_avx_dy()
{
    __m256 t5,t6,t3,t4;
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x=0;x<p;x++)
    {
        m_reset();
    QueryPerformanceCounter(&t1);
    #pragma omp parallel num_threads(NUM_THREADS)
	for (int k = 0; k < N; k++)
        {
        t5 = _mm256_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 8 < N; j += 8)
		{
			t6 = _mm256_loadu_ps(&m[k][j]);
			t6 = _mm256_div_ps(t6, t5);
			_mm256_storeu_ps(&m[k][j], t6);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;
        #pragma omp for schedule(dynamic)
		for (int i = k + 1; i < N; i++)
		{
			t5 = _mm256_set1_ps(m[i][k]);
            int j = k + 1;
			for (; j + 8 < N; j += 8)
			{
				t6 = _mm256_loadu_ps(&m[k][j]);
				t3 = _mm256_loadu_ps(&m[i][j]);
				t4 = _mm256_mul_ps(t5,t6);
				t3 = _mm256_sub_ps(t3, t4);
				_mm256_storeu_ps(&m[i][j], t3);
			}
			for (; j < N; j++)
				m[i][j] -= m[i][k] * m[k][j];
			m[i][k] = 0;
		}
	}
	QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
    }
    cout << "openmp_avx_dy_time: " << sumtime/p << "ms" << endl;
}



void print()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
            cout<<m[i][j]<<" ";
        cout<<endl;
    }
}


int main()
{

    m = new float* [N];
	for (int i = 0; i < N; i++)
    {
		m[i] = new float[N];
	}

    //normal();
    //openmp();
    openmp_dynamic();
    //openmp_simd();
    //openmp_simd_dy();
    //openmp_avx();
    //openmp_avx_dy();
	return 0;
}
