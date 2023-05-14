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




//并行
const int worker_count = 7;
sem_t sem_main;
sem_t sem_workerstart[worker_count];
sem_t sem_workerend[worker_count];

typedef struct
{
	int k;
	int t_id;
} threadParam_t;

//线程函数:
void* threadFunc(void* param)
{

	threadParam_t* p = (threadParam_t*)param;
	int k = p->k;
	int t_id = p->t_id;

	for (int k = 0; k < N; k++)
    {
        sem_wait(&sem_workerstart[t_id]); // 阻塞
        //循环划分任务
        for (int i = k + 1 + t_id; i < N; i += worker_count)
        {

            for (int j = k + 1; j < N; ++j)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];

            m[i][k] = 0.0;
        }
        sem_post(&sem_main);
        sem_wait(&sem_workerend[t_id]);
    }

	pthread_exit(NULL);
}


void pthread_static()
{

    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
    //初始化信号量
    sem_init(&sem_main, 0, 0);

    for (int i = 0; i < worker_count; i++)
    {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }

    //创建线程

    pthread_t handles[worker_count]; // 创建对应的 Handle
    threadParam_t param[worker_count]; // 创建对应的线程数据结构

	for (int t_id = 0; t_id < worker_count; t_id++)
    {

        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (void *)&param[t_id]);
    }

    for (int k = 0; k < N; k++)
    {
        //主线程做除法操作
        for (int j = k + 1; j < N; j++)
        {
            m[k][j] = m[k][j] / m[k][k];
        }
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

    for (int t_id = 0; t_id < worker_count; t_id++)
    {
        pthread_join(handles[t_id],NULL);
    }

    //销毁所有信号量
    sem_destroy(&sem_main);
    sem_destroy(sem_workerend);
    sem_destroy(sem_workerstart);

     QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "static_time: " << sumtime/p << "ms" << endl;
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

    pthread_static();
    //print();
	return 0;
}
