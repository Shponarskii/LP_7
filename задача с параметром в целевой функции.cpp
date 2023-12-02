#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
using namespace std;

void readfile(int**& A, int &m, int &n, int *&b, int *&c,int *&c1)  // Считывание из файла 
{
    string l;
    ifstream f("input.txt");
    if (f.is_open())
    {
        getline(f, l);
        for (int i = 0; i < l.size(); i++)
            if (l[i] == ' ')
                n++;
        while (!f.eof())
        {
            getline(f, l);
            m++;
        }
        m -= 3;
        f.seekg(0);
        A = new int* [m];
        b = new int[m];
        c = new int[n];
        c1 = new int[n];
        for (int i = 0; i < m; i++)
            A[i] = new int[n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                f >> A[i][j];
        for (int i = 0; i < m; i++)
            f >> b[i];
        for (int i = 0; i < n; i++)
            f >> c[i];
        for (int i = 0; i < n; i++)
            f >> c1[i];
    }
    f.close();
}

void createM(int**& p, int m,int n,int a)  // Создание матрицы
{
    p = new int* [m];
    for (int i = 0; i < m; i++)
        p[i] = new int[n];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            p[i][j] = a;
}

void createA(int*& p, int n, int a)  // Создание массива
{
    p = new int[n];
    for (int i = 0; i < n; i++)
        p[i] = a;
}

void printM(int** A, int m, int n)  // Вывод матрицы
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}

void printA(int* a, int n)  // Вывод массива
{
    for (int i = 0; i < n; i++)
        cout << a[i] << " ";
    cout << endl;
}

void copy_M(int** B, int**& B1, int m, int n)  // Копирование матрицы
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            B1[i][j] = B[i][j];
}

void kol(int ** A, int**& p, int m, int n)  // Подсчет кол-ва единичных столбцов матрицы
{
    for (int j = 0; j < n; j++)
    {
        int k0 = 0, k1 = 0, q1 = 0;
        for (int i = 0; i < m; i++)
        {
            if (A[i][j] == 1)
            {
                k1++;
                q1 = i;
            }
            if (A[i][j] == 0)
                k0++;
        }
        if (k1 == 1 && k0 == (m - 1))
        {
            p[q1][0] = 1;
            p[q1][1] = j;
        }
    }
}

int Evklid(int a, int b)  // Наибольший общий делитель двух чисел
{
    while (a != 0 && b != 0)
    {
        if (a >= b)
            a %= b;
        else
            b %= a;
    }
    return a + b;
}

int Nok(int a, int b)  // Наименьшее общее кратное
{
    return a * b / Evklid(a, b);
}

void reduce(int &a, int &b)  // Сокращение дроби
{
    int t = abs(a);
    int r = Evklid(t, b);
    t = t / r;
    b = b / r;
    a = a < 0 ? -t : t;
}

void out(int a, int b)  // Вывод дроби
{
    if (b == 1 || a == 0) cout << a;
    else cout << a << "/" << b;
}

void output(int m, int n, int** p, int* c, int* b, int** B, int** A, int* d, int* d1, int* c1, int* d2, int* d3)  // Вывод симплекс-таблицы
{
    for (int i = 0; i < m; i++)
    {
        cout << "A" << p[i][1] + 1 << setw(5) << c[p[i][1]] << setw(5) << c1[p[i][1]] << setw(5);
        out(b[i], B[i][n]);
        for (int j = 0; j < n; j++)
        {
            cout << setw(5);
            out(A[i][j], B[i][j]);
        }
        cout << endl;
    }
    cout << "           d" << setw(5);
    out(d[n], d1[n]);
    for (int i = 0; i < n; i++)
    {
        cout << setw(5);
        out(d[i], d1[i]);
    }
    cout << endl;
    cout << "          d1" << setw(5);
    out(d2[n], d3[n]);
    for (int i = 0; i < n; i++)
    {
        cout << setw(5);
        out(d2[i], d3[i]);
    }
    cout << endl << endl;
}

void Gauss(int m, int n, int**& B, int**& A, int*& b, int i, int w)  // метод Гаусса
{
    if (A[i][w] < 0)
    {
        for (int j = 0; j < n; j++)
            A[i][j] = -A[i][j];
        b[i] = -b[i];
    }
    int g0 = B[i][w];
    for (int j = 0; j < n + 1; j++)
        B[i][j] *= A[i][w];
    for (int j = 0; j < n; j++)
        A[i][j] *= g0;
    b[i] *= g0;
    for (int s = 0; s < m; s++)
    {
        int r1 = A[s][w];
        int r2 = B[s][w];
        for (int j = 0; j < n; j++)
            if (s != i)
            {
                int t0 = Nok(B[s][j], r2 * B[i][j]);
                A[s][j] = A[s][j] * t0 / B[s][j] - (r1 * A[i][j] * t0) / (r2 * B[i][j]);
                B[s][j] = t0;
            }
        if (s != i)
        {
            int t1 = Nok(B[s][n], r2 * B[i][n]);
            b[s] = b[s] * t1 / B[s][n] - (r1 * b[i] * t1) / (r2 * B[i][n]);
            B[s][n] = t1;
        }
    }
    for (int s = 0; s < m; s++)
        for (int j = 0; j < n; j++)
            reduce(A[s][j], B[s][j]);
    for (int s = 0; s < m; s++)
        reduce(b[s], B[s][n]);
}

int beg(int**& A, int*& b, int**& p, int**& B, int &m, int n)  // Формирование начального базиса
{
    for (int i = 0; i < m; i++)
    {
        while (p[i][0] == 0)
        {
            int w = 0;
            int flag = 0;
            for (int j = 0; j < n; j++)
                if (A[i][j] != 0)
                {
                    w = j;
                    flag = 1;
                    break;
                }
            if (flag == 0)
            {
                if (b[i] != 0)
                {
                    cout << "Задача неразрешима при любом значении q";
                    return -1;
                }
                else
                {
                    if (i == (m - 1))
                        break;
                    else
                    {
                        swap(A[i], A[m - 1]);
                        swap(B[i], B[m - 1]);
                        swap(p[i], p[m - 1]);
                        swap(b[i], b[m - 1]);
                        m--;
                    }
                }
            }
            else
            {
                Gauss(m, n, B, A, b, i, w);
                p[i][0] = 1;
                p[i][1] = w;
            }
        }
    }
    return 1;
}

int escape(int m, int*& b, int n, int**& A, int**& B, int**& p)  // Избавление от отрицательных коэффициентов в столбце b
{
    int Min = 1, q1 = 0;
    for (int i = 0; i < m; i++)
        if (b[i] < 0)
            if (b[i] < Min)
            {
                Min = b[i];
                q1 = i;
            }
    while (Min != 1)
    {
        int Min1 = 1, q2 = 0;
        for (int j = 0; j < n; j++)
            if (A[q1][j] < 0)
                if (A[q1][j] < Min1)
                {
                    Min1 = A[q1][j];
                    q2 = j;
                }
        if (Min1 == 1)
        {
            cout << "Задача неразрешима при любом значении q";
            return -1;
        }

        for (int j = 0; j < n; j++)
            A[q1][j] = -A[q1][j];
        b[q1] = -b[q1];
        Gauss(m, n, B, A, b, q1, q2);

        p[q1][1] = q2;
        Min = 1, q1 = 0;
        for (int i = 0; i < m; i++)
            if (b[i] < 0)
                if (b[i] < Min)
                {
                    Min = b[i];
                    q1 = i;
                }
    }
    return 1;
}

void delta(int m, int n, int** A, int** B, int* c, int** p, int* b, int*& d, int*& d1, int* c1, int r1, int r2)  // Вычисление дельт
{
    for (int j = 0; j < n; j++)
    {
        int Max = -1;
        for (int i = 0; i < m; i++)
            if (B[i][j] > Max)
                Max = B[i][j];
        for (int i = 0; i < m; i++)
            d[j] += (c[p[i][1]] * r2 + c1[p[i][1]] * r1) * A[i][j] * Max / B[i][j];
        d[j] -= (c[j] * r2 + c1[j] * r1) * Max;
        d1[j] = r2 * Max;
    }

    int Max = -1;
    for (int i = 0; i < m; i++)
        if (B[i][n] > Max)
            Max = B[i][n];
    for (int i = 0; i < m; i++)
        d[n] += (c[p[i][1]] * r2 + c1[p[i][1]] * r1) * b[i] * Max / B[i][n];
    d1[n] = r2 * Max;

    for (int i = 0; i < n + 1; i++)
        reduce(d[i], d1[i]);
}

int rules(int n, int*& d, int*& d1, int m, int**& A, int u, int*& b, int**& p, int**& B, int* c, int* c1,int r1,int r2)  // Правила симплекс-метода
{
    for (int i = 0; i < n + 1; i++)
    {
        d[i] = 0;
        d1[i] = 1;
    }
    int flag = 0;
    double Min1 = DBL_MAX;
    int r = 0;
    for (int i = 0; i < m; i++)
    {
        if (A[i][u] > 0)
        {
            flag++;
            if ((b[i] * B[i][u]) / float((A[i][u] * B[i][n])) < Min1)
            {
                Min1 = (b[i] * B[i][u]) / float(A[i][u] * B[i][n]);
                r = i;
            }
        }
    }
    if (flag == 0)
        return -1;
    else
    {
        p[r][1] = u;
        Gauss(m, n, B, A, b, r, u);
        delta(m, n, A, B, c, p, b, d, d1, c1, r1, r2);
    }
    return 3;
}

int Simplex(int n, int*& d, int*& d1, int m, int**& A, int*& b, int**& p, int**& B, int* c, int* c1,int r1,int r2,int *d2,int *d3)  // Симплекс-метод
{
    double Min = 1;
    int u = 0;
    for (int j = 0; j < n; j++)
        if (d[j] / float(d1[j]) < Min)
        {
            Min = d[j];
            u = j;
        }
    while (Min < 0)
    {
        if (rules(n, d, d1, m, A, u, b, p, B, c, c1, r1, r2) == -1)
            return -1;
        Min = 1;
        u = 0;
        for (int j = 0; j < n; j++)
            if (d[j] / float(d1[j]) < Min)
            {
                Min = d[j];
                u = j;
            }
        for (int i = 0; i < n + 1; i++)
        {
            d2[i] = 0;
            d3[i] = 1;
        }
        delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
        output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);
    }
    return 1;
}

void otvet(int m, int n, int** p, int* b, int** B, int* d, int* d1, int* d2, int* d3)  // Вывод оптимального решения
{
    cout << endl;
    cout << "Оптимальный план : x = (";
    for (int j = 0; j < n; j++)
    {
        int fl = 0;
        for (int i = 0; i < m; i++)
            if (j == p[i][1])
            {
                if (j == n - 1)
                    out(b[i], B[i][n]);
                else
                {
                    out(b[i], B[i][n]);
                    cout << " ; ";
                }
                fl++;
                break;
            }
        if (fl == 0)
        {
            if (j == n - 1) cout << 0;
            else cout << 0 << " ; ";
        }
    }
    cout << "), значение целевой функции : l = ";
    out(d[n], d1[n]);
    if (d2[n] > 0)
    {
        cout << "+";
        if (d3[n] != 1)
        {
            cout << "(";
            out(d2[n], d3[n]);
            cout << ")";
        }
        else out(d2[n], d3[n]);
        cout << "q";
    }
    else
    {
        if (d2[n] < 0)
        {
            if (d3[n] != 1)
            {
                cout << "-(";
                cout << abs(d2[n]) << "/" << d3[n] << ")";
            }
            else out(d2[n], d3[n]);
            cout << "q";
        }
    }
    cout << endl << endl;
}

void parametr1(int m, int n, int**& A, int**& B, int* c1, int**& p, int*& b, int*& d2, int*& d3, int*& d, int*& d1, int* c, int& s1, int& t1, int& t2)  // Рассматриваем q, меньшее λ_1/λ_2
{
    int λ_1, λ_2 = 1;
    for (int i = 0; i < n + 1; i++)
    {
        d2[i] = 0;
        d3[i] = 1;
    }
    delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
    int flag2 = 0;
    int u1 = -INT_MAX, u2 = 1;
    int q1 = -1;
    for (int i = 0; i < n; i++)
        if (d2[i] > 0)
        {
            flag2++;
            if (-d[i] * d3[i] / float(d1[i] * d2[i]) > u1 / float(u2))
            {
                int r1 = -d[i] * d3[i], r2 = d1[i] * d2[i];
                if (r2 < 0)
                {
                    r1 = -r1;
                    r2 = -r2;
                }
                reduce(r1, r2);
                u1 = r1;
                u2 = r2;
                q1 = i;
            }
        }
    if (flag2 == 0)
        λ_1 = -INT_MAX;
    else
    {
        λ_1 = u1;
        λ_2 = u2;
    }
    if (λ_1 == -INT_MAX)
    {
        int λ1, λ2 = 1;
        int flag1 = 0;
        int u3 = INT_MAX, u4 = 1;
        int q2 = -1;
        for (int i = 0; i < n; i++)
            if (d2[i] < 0)
            {
                flag1++;
                if ((-d[i] * d3[i]) / float(d1[i] * d2[i]) < u3 / float(u4))
                {
                    int r1 = -d[i] * d3[i], r2 = d1[i] * d2[i];
                    if (r2 < 0)
                    {
                        r1 = -r1;
                        r2 = -r2;
                    }
                    reduce(r1, r2);
                    u3 = r1;
                    u4 = r2;
                    q2 = i;
                }
            }
        if (flag1 == 0)
            λ1 = INT_MAX;
        else
        {
            λ1 = u3;
            λ2 = u4;
        }
        output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);
        cout << "При q в (-inf ; ";
        if (λ1 != INT_MAX)
        {
            out(λ1, λ2);
            cout << "]:";
        }
        else cout << "+inf)";
        otvet(m, n, p, b, B, d, d1, d2, d3);
    }
    else
    {
        t1 = λ_1;
        t2 = λ_2;
        s1 = 1;
    }
    while (λ_1 != -INT_MAX)
    {
        if (rules(n, d, d1, m, A, q1, b, p, B, c, c1, 0, 1) == -1)
        {
            cout << "При q в (-inf ; ";
            out(λ_1, λ_2);
            cout << ") целевая функция не ограничена сверху на допустимом множестве" << endl << endl;
            λ_1 = -INT_MAX;
        }
        else
        {
            int y1 = λ_1, y2 = λ_2;
            for (int i = 0; i < n + 1; i++)
            {
                d2[i] = 0;
                d3[i] = 1;
            }
            delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
            output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);

            flag2 = 0, u1 = -INT_MAX, u2 = 1, q1 = -1;
            for (int i = 0; i < n; i++)
                if (d2[i] > 0)
                {
                    flag2++;
                    if (-d[i] * d3[i] / float(d1[i] * d2[i]) > u1 / float(u2))
                    {
                        int r1 = -d[i] * d3[i], r2 = d1[i] * d2[i];
                        if (r2 < 0)
                        {
                            r1 = -r1;
                            r2 = -r2;
                        }
                        reduce(r1, r2);
                        u1 = r1;
                        u2 = r2;
                        q1 = i;
                    }
                }
            if (flag2 == 0)
                λ_1 = -INT_MAX;
            else
            {
                λ_1 = u1;
                λ_2 = u2;
            }

            cout << "При q в ";
            if (λ_1 == -INT_MAX)
                cout << "(-inf";
            else
            {
                cout << "[";
                out(λ_1, λ_2);
            }
            cout << " ; ";
            out(y1, y2);
            cout << "]:";
            otvet(m, n, p, b, B, d, d1, d2, d3);
        }
    }
}

void parametr2(int m, int n, int**& A, int**& B, int* c1, int**& p, int*& b, int*& d2, int*& d3, int*& d, int*& d1, int* c, int& s1, int& t1, int& t2)  // Рассматриваем q, большее λ1/λ2
{
    int λ1, λ2 = 1;
    for (int i = 0; i < n + 1; i++)
    {
        d[i] = 0;
        d1[i] = 1;
    }
    delta(m, n, A, B, c, p, b, d, d1, c1, 0, 1);
    for (int i = 0; i < n + 1; i++)
    {
        d2[i] = 0;
        d3[i] = 1;
    }
    delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
    int flag1 = 0;
    int u3 = INT_MAX, u4 = 1;
    int q2 = -1;
    for (int i = 0; i < n; i++)
        if (d2[i] < 0)
        {
            flag1++;
            if ((-d[i] * d3[i]) / float(d1[i] * d2[i]) < u3 / float(u4))
            {
                int r1 = -d[i] * d3[i], r2 = d1[i] * d2[i];
                if (r2 < 0)
                {
                    r1 = -r1;
                    r2 = -r2;
                }
                reduce(r1, r2);
                u3 = r1;
                u4 = r2;
                q2 = i;
            }
        }
    if (flag1 == 0)
        λ1 = INT_MAX;
    else
    {
        λ1 = u3;
        λ2 = u4;
    }

    if (s1 == 1)
    {
        output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);
        cout << "При q в [";
        out(t1, t2);
        cout << " ; ";
        if (λ1 == INT_MAX)
            cout << "inf):";
        else
        {
            out(λ1, λ2);
            cout << "]:";
        }
        otvet(m, n, p, b, B, d, d1, d2, d3);
    }

    while (λ1 != INT_MAX)
    {
        if (rules(n, d, d1, m, A, q2, b, p, B, c,c1,0,1) == -1)
        {
            cout << "При q в (";
            out(λ1, λ2);
            cout << " ; +inf) целевая функция не ограничена сверху на допустимом множестве" << endl << endl;
            λ1 = INT_MAX;
        }
        else
        {
            int y1 = λ1, y2 = λ2;
            for (int i = 0; i < n + 1; i++)
            {
                d2[i] = 0;
                d3[i] = 1;
            }
            delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
            output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);

            flag1 = 0, u3 = INT_MAX, u4 = 1, q2 = -1;
            for (int i = 0; i < n; i++)
                if (d2[i] < 0)
                {
                    flag1++;
                    if ((-d[i] * d3[i]) / float(d1[i] * d2[i]) < u3 / float(u4))
                    {
                        int r1 = -d[i] * d3[i], r2 = d1[i] * d2[i];
                        if (r2 < 0)
                        {
                            r1 = -r1;
                            r2 = -r2;
                        }
                        reduce(r1, r2);
                        u3 = r1;
                        u4 = r2;
                        q2 = i;
                    }
                }
            if (flag1 == 0)
                λ1 = INT_MAX;
            else
            {
                λ1 = u3;
                λ2 = u4;
            }
            cout << "При q в [";
            out(y1, y2);
            cout << " ; ";
            if (λ1 == INT_MAX)
                cout << "+inf):";
            else
            {
                out(λ1, λ2);
                cout << "]:";
            }
            otvet(m, n, p, b, B, d, d1, d2, d3);
        }
    }
}

int main()
{
    setlocale(LC_ALL, "rus");
    int** A, * b, * c, * c1, ** p, ** B, * d, * d1, * d2, * d3;
    int m = 1, n = 1;
    int s1 = 0, t1 = 0, t2 = 0;
    int r1 = -2, r2 = 1;
    int flag1 = 0, flag2 = 0;

    readfile(A, m, n, b, c, c1);
    createM(p, m, 2, 0);
    createM(B, m, n + 1, 1);
    createA(d, n + 1, 0);
    createA(d1, n + 1, 1);
    createA(d2, n + 1, 0);
    createA(d3, n + 1, 1);
    kol(A, p, m, n);

    if (beg(A, b, p, B, m, n) == -1)
        return 2;

    if (escape(m, b, n, A, B, p) == -1)
        return 2;
 
    cout << "Возьмем q = 0 и проведем решение симплекс-методом: " << endl << endl;
    delta(m, n, A, B, c, p, b, d, d1, c1, r1, r2);
    delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
    output(m, n, p, c, b, B, A, d, d1, c1, d2, d3);

    int** A1, * b1, ** p1, ** B1;

    createM(A1, m, n, 0);
    createM(B1, m, n + 1, 1);
    createM(p1, m, 2, 0);
    createA(b1, m, 0);

    copy_n(b, m, b1);
    copy_M(B, B1, m, n + 1);
    copy_M(p, p1, m, 2);
    copy_M(A, A1, m, n);

    if (Simplex(n, d, d1, m, A1, b1, p1, B1, c, c1, r1, r2, d2, d3) == -1)
    {
        cout << "При q = ";
        out(r1, r2);
        cout << " целевая функция задачи не ограничена на допустимом множестве" << endl << endl;
        for (int i = 0; i < n + 1; i++)
        {
            d[i] = 0;
            d1[i] = 1;
        }
        for (int i = 0; i < n + 1; i++)
        {
            d2[i] = 0;
            d3[i] = 1;
        }
        delta(m, n, A1, B1, c, p1, b1, d, d1, c1, 0, 1);
        delta(m, n, A1, B1, c1, p1, b1, d2, d3, c1, 0, 1);
        output(m, n, p1, c, b1, B1, A1, d, d1, c1, d2, d3);

        int k = -1;
        for (int i = 0; i < n; i++)
            if ((d[i] / float(d1[i]) + (r1 / float(r2)) * (d2[i] / float(d3[i]))) < 0)
                k = i;
        if (d2[k] == 0)
        {
            cout << "При q в (-inf ; +inf) целевая функция не ограничена на допустимом множестве" << endl;
            return -1;
        }
        else
        {
            r1 = -d[k] * d3[k], r2 = d1[k] * d2[k];
            if (r2 < 0)
            {
                r1 = -r1;
                r2 = -r2;
            }
            reduce(r1, r2);
            cout << "Возьмем q = ";
            out(r1, r2);
            cout << " и проведем решение симплекс-методом" << endl << endl;
            copy_n(b, m, b1);
            copy_M(B, B1, m, n + 1);
            copy_M(p, p1, m, 2);
            copy_M(A, A1, m, n);
            for (int i = 0; i < n + 1; i++)
            {
                d[i] = 0;
                d1[i] = 1;
            }
            delta(m, n, A1, B1, c, p1, b1, d, d1, c1, r1, r2);
            if (d2[k] > 0) flag1++;
            else flag2++;
            int gr = d2[k];
            while (Simplex(n, d, d1, m, A1, b1, p1, B1, c, c1, r1, r2, d2, d3) == -1)
            {
                for (int i = 0; i < n + 1; i++)
                {
                    d2[i] = 0;
                    d3[i] = 1;
                }
                delta(m, n, A, B, c1, p, b, d2, d3, c1, 0, 1);
                int s = -1;
                for (int i = 0; i < n; i++)
                    if ((d[i]/d1[i]+(r1/r2)*(d2[i]/d3[i])) < 0)
                        s = i;
                copy_n(b, m, b1);
                copy_M(B, B1, m, n + 1);
                copy_M(p, p1, m, 2);
                copy_M(A, A1, m, n);
                if (gr > 0)
                {
                    flag1++;
                    if (d2[s] <= 0)
                    {
                        cout << "При q в (-inf ; +inf) целевая функция не ограничена на допустимом множестве" << endl;
                        return -1;
                    }
                }
                if (gr < 0)
                {
                    flag2++;
                    if (d2[s] >= 0)
                    {
                        cout << "При q в (-inf ; +inf) целевая функция не ограничена на допустимом множестве" << endl;
                        return -1;
                    }
                }
                r1 = -d[s] * d3[s], r2 = d1[s] * d2[s];
                if (r2 < 0)
                {
                    r1 = -r1;
                    r2 = -r2;
                }
                reduce(r1, r2);
            }

            for (int i = 0; i < n + 1; i++)
            {
                d2[i] = 0;
                d3[i] = 1;
            }
            delta(m, n, A1, B1, c1, p1, b1, d2, d3, c1, 0, 1);
            for (int i = 0; i < n + 1; i++)
            {
                d[i] = 0;
                d1[i] = 1;
            }
            delta(m, n, A1, B1, c, p1, b1, d, d1, c1, 0, 1);
            output(m, n, p1, c, b1, B1, A1, d, d1, c1, d2, d3);

            if (flag2 > 0)
            {
                cout << "При q в (";
                out(r1, r2);
                cout << " ; +inf) целевая функция не ограничена сверху на допустимом множестве" << endl << endl;
            }
            if (flag1 > 0)
            {
                cout << "При q в (-inf ; ";
                out(r1, r2);
                cout << ") целевая функция не ограничена сверху на допустимом множестве" << endl << endl;
            }

            int** A2, * b2, ** p2, ** B2;

            createM(A2, m, n, 0);
            createM(B2, m, n + 1, 1);
            createM(p2, m, 2, 0);
            createA(b2, m, 0);

            copy_n(b1, m, b2);
            copy_M(B1, B2, m, n + 1);
            copy_M(p1, p2, m, 2);
            copy_M(A1, A2, m, n);

            cout << "-----------------------------------------------------------" << endl << endl;
            parametr1(m, n, A1, B1, c1, p1, b1, d2, d3, d, d1, c, s1, t1, t2);
            parametr2(m, n, A2, B2, c1, p2, b2, d2, d3, d, d1, c, s1, t1, t2);
        }
    }

    else
    {
        int** A2, * b2, ** p2, ** B2;
        cout << "При q = ";
        out(r1, r2);
        cout << " задача имеет решение" << endl << endl;

        createM(A2, m, n, 0);
        createM(B2, m, n + 1, 1);
        createM(p2, m, 2, 0);
        createA(b2, m, 0);

        copy_n(b1, m, b2);
        copy_M(B1, B2, m, n + 1);
        copy_M(p1, p2, m, 2);
        copy_M(A1, A2, m, n);

        cout << "-----------------------------------------------------------" << endl << endl;
        parametr1(m, n, A1, B1, c1, p1, b1, d2, d3, d, d1, c, s1, t1, t2);
        parametr2(m, n, A2, B2, c1, p2, b2, d2, d3, d, d1, c, s1, t1, t2);
    }
}
