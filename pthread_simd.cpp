#include<iostream>
#include<pthread.h>
#include <windows.h>
#include<semaphore.h>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX、AVX2
using namespace std;
const int N=1000;
const int p=1;
float** m;
//float** n;
LARGE_INTEGER freq, t1, t2;

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

void normal()//普通
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
const int worker_count = 7;
sem_t sem_main;
sem_t sem_workerstart[worker_count];
sem_t sem_workerend[worker_count];

void unalign_sse()//不对齐SSE
{
    __m128 t5,t6,t3,t4;
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x=0;x<p;x++)
    {
        m_reset();
    QueryPerformanceCounter(&t1);
	for (int k = 0; k < N; k++) //除法优化
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
    cout << "unalign_sse_time: " << sumtime/p << "ms" << endl;
}





typedef struct
{
	int k;
	int t_id;
	int worker_count;
} threadParam_t;






//pthread+SSE算法:
//线程函数:
void* threadFuncSSE(void* param)
{
	threadParam_t* p = (threadParam_t*)param;
	int k = p->k;
	int t_id = p->t_id;
	int worker_count = p->worker_count;
for (int k = 0; k < N; k++)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞
        //循环划分任务
	for (int i = k + 1 + t_id; i < N; i += worker_count)
    {
		__m128 t3 = _mm_loadu_ps(&m[i][k]);
		int j = k + 1;
		for (; j + 4 <= N; j += 4)
        {
			__m128 t4 = _mm_loadu_ps(&m[k][j]);
			__m128 t5 = _mm_loadu_ps(&m[i][j]);
			__m128 t6 = _mm_mul_ps(t4, t3);
			t5 = _mm_sub_ps(t5, t6);
			_mm_storeu_ps(&m[i][j], t5);
		}
		for (; j < N; j++)//结尾处理
				m[i][j] -= m[i][k] * m[k][j];
		m[i][k] = 0;
	}
	sem_post(&sem_main);
    sem_wait(&sem_workerend[t_id]);
    }
	pthread_exit(NULL);
}

void Pthread_SSE()
{
     double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);



    sem_init(&sem_main, 0, 0);

    for (int i = 0; i < worker_count; i++)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }



    pthread_t handles[worker_count]; // 创建对应的 Handle
    threadParam_t param[worker_count]; // 创建对应的线程数据结构

	for (int t_id = 0; t_id < worker_count; t_id++)
    {

        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFuncSSE, (void *)&param[t_id]);
    }

	for (int k = 0; k < N; ++k)
        {
		//主线程做除法操作,使用SIMD优化;
		__m128 vt = _mm_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 4 <= N; j += 4)
        {
			__m128 va = _mm_loadu_ps(&m[k][j]);
			va = _mm_div_ps(va, vt);
			_mm_storeu_ps(&m[k][j], va);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;

        //开始唤醒工作线程
        for (int t_id = 0; t_id < worker_count; t_id++)
        {
            sem_post(&sem_workerstart[t_id]);
        }
        //主线程进入睡眠
        for(int t_id = 0;t_id < worker_count; t_id++)
        {
            sem_wait(&sem_main);
        }

        // 再次唤醒工作线程
        for (int t_id = 0; t_id < worker_count; t_id++)
        {
            sem_post(&sem_workerend[t_id]);
        }
    }
		//主线程挂起等待所有的工作线程完成此轮消去工作
		for (int t_id = 0; t_id < worker_count; t_id++)
        {
			pthread_join(handles[t_id], NULL);
		}

	sem_destroy(&sem_main);
    sem_destroy(sem_workerend);
    sem_destroy(sem_workerstart);

	 QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "simd_time: " << sumtime/p << "ms" << endl;
}


void* threadFuncAVX256(void* param)
{
	threadParam_t* p = (threadParam_t*)param;
	int k = p->k;
	int t_id = p->t_id;
	int worker_count = p->worker_count;
for (int k = 0; k < N; k++)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞
        //循环划分任务
	for (int i = k + 1 + t_id; i < N; i += worker_count)
    {
		__m256 t3 = _mm256_loadu_ps(&m[i][k]);
		int j = k + 1;
		for (; j + 8 <= N; j += 8)
        {
            __m256 t4 = _mm256_loadu_ps(&m[k][j]);
			__m256 t5 = _mm256_loadu_ps(&m[i][j]);
			__m256 t6 = _mm256_mul_ps(t4, t3);
			t5 = _mm256_sub_ps(t5, t6);
			_mm256_storeu_ps(&m[i][j], t5);
		}
		for (; j < N; j++)//结尾处理
				m[i][j] -= m[i][k] * m[k][j];
		m[i][k] = 0;
	}
	sem_post(&sem_main);
    sem_wait(&sem_workerend[t_id]);
    }
	pthread_exit(NULL);
}

void Pthread_AVX256()
{
     double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
	for (int k = 0; k < N; ++k)
        {
		//主线程做除法操作,使用SIMD优化;
		__m256 vt = _mm256_set1_ps(m[k][k]);
		int j = k + 1;
		for (; j + 8 <= N; j += 8)
        {
			__m256 va = _mm256_loadu_ps(&m[k][j]);
			va = _mm256_div_ps(va, vt);
			_mm256_storeu_ps(&m[k][j], va);
		}
		for (; j < N; j++)
			m[k][j] /= m[k][k];
		m[k][k] = 1.0;

		 //开始唤醒工作线程
        for (int t_id = 0; t_id < worker_count; t_id++)
        {
            sem_post(&sem_workerstart[t_id]);
        }
        //主线程进入睡眠
        for(int t_id = 0;t_id < worker_count; t_id++)
        {
            sem_wait(&sem_main);
        }

        // 再次唤醒工作线程
        for (int t_id = 0; t_id < worker_count; t_id++)
        {
            sem_post(&sem_workerend[t_id]);
        }
        }
		//主线程挂起等待所有的工作线程完成此轮消去工作
		for (int t_id = 0; t_id < worker_count; t_id++)
        {
			pthread_join(handles[t_id], NULL);
		}

	sem_destroy(&sem_main);
    sem_destroy(sem_workerend);
    sem_destroy(sem_workerstart);

	 QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "simd_time: " << sumtime/p << "ms" << endl;
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

    normal();

    Pthread_SSE();
   // Pthread_AVX256();
	return 0;
}
