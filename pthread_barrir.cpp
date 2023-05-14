#include<iostream>
#include<pthread.h>
#include<semaphore.h>
#include <windows.h>
using namespace std;
const int N=2500;
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





typedef struct
{
	int k;
	int t_id;
} threadParam_t;

const int worker_count = 7;

//barrier算法
pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;
void* _threadFunc(void* param)
{

	threadParam_t* p = (threadParam_t*)param;

	int k = p->k; //消去的轮次
	int t_id = p->t_id; //线程编号

	for (int k = 0; k < N; k++)
    {
		if (t_id == 0)
		{
            for (int j = k + 1; j < N; j++)
            {
				m[k][j] = m[k][j] / m[k][k];
			}
            m[k][k] = 1.0;
        }
        pthread_barrier_wait(&barrier_Divsion);

        for (int i = k + 1 + t_id; i < N; i += worker_count)
        {
            //消去
            for (int j = k + 1; j < N; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k] = 0.0;
        }
        pthread_barrier_wait(&barrier_Elimination);
    }

	pthread_exit(NULL);
}


void barrier()
{
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);

    pthread_barrier_init(&barrier_Divsion, NULL, worker_count);
    pthread_barrier_init(&barrier_Elimination, NULL, worker_count);

    //创建线程

    pthread_t handles[worker_count];
    threadParam_t param[worker_count];

	for (int t_id = 0; t_id < worker_count; t_id++)
    {

        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, _threadFunc, (void *)&param[t_id]);
    }

    for (int t_id = 0; t_id < worker_count; t_id++)
    {
        pthread_join(handles[t_id],NULL);
    }

    //销毁barrier
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);

    QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "barrier_time: " << sumtime/p << "ms" << endl;
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

    barrier();
    //print();
	return 0;
}
