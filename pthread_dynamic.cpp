
#include<iostream>
#include<iomanip>
#include<pthread.h>
#include <windows.h>
using namespace std;
const int N=1500;
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
	int k;//消去的轮次
	int t_id;//线程id
	int worker_count;

} threadParam_t;

//线程函数:
void* threadFunc(void* param)
{
	//获取参数:
	threadParam_t* p = (threadParam_t*)param;
	int k = p->k; //消去的轮次
	int t_id = p->t_id; //线程编号
	int worker_count = p->worker_count;

	for (int i = k + 1 + t_id; i < N; i += worker_count)
    {
		for (int j = k + 1; j < N; ++j)
		{
			m[i][j] = m[i][j] - m[i][k] * m[k][j];
		}
		m[i][k] = 0;
	}
	pthread_exit(NULL);
}

void dynamic()
{
    double sumtime=0;
    QueryPerformanceFrequency(&freq);
    for(int x = 0 ;x < p ; x++ )
    {
    m_reset();
    QueryPerformanceCounter(&t1);
	for (int k = 0; k < N; ++k)
    {
		//主线程做除法操作
		for (int j = k + 1; j < N; j++)
		{
			m[k][j] = m[k][j] / m[k][k];
		}
		m[k][k] = 1.0;

		//工作线程进行消去操作
		int worker_count = 7;
		pthread_t* handles = new pthread_t[worker_count];
		threadParam_t* param = new threadParam_t[worker_count];

		//分配任务
		for (int t_id = 0; t_id < worker_count; t_id++)
        {
			param[t_id].k = k;
			param[t_id].t_id = t_id;
			param[t_id].worker_count = worker_count;
		}

		//创建线程
		for (int t_id = 0; t_id < worker_count; t_id++)
        {
			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
		}

		//主线程挂起等待所有的工作线程完成此轮消去工作
		for (int t_id = 0; t_id < worker_count; t_id++)
        {
			pthread_join(handles[t_id], NULL);
		}
	}
	QueryPerformanceCounter(&t2);
    sumtime += (t2.QuadPart - t1.QuadPart) * 1000.0 / freq.QuadPart;
     }
    cout << "dynamic_time: " << sumtime/p << "ms" << endl;
}

int main()
{

    m = new float* [N];
	for (int i = 0; i < N; i++)
    {
		m[i] = new float[N];
	}


	normal();
    dynamic();


	return 0;
}

