#include <iostream>
#include <vector>
#include <cmath> 
using namespace std;


//функции для N и для M
double func_N(double x) {
    return x * x * x * x ;
}

double func_M(double x) {
    return x * x * x ;
}

//трапеция
double trapezoid(double a, double b, double n, double (*func)(double)) {
    double h = (b - a) / (n);
    double primit = 0.0;
    std::vector<double> y;

    double x = a;
    while (x <= b)
    {
        y.push_back((*func)(x));
        x += h;
    }

    for (int i = 0; i < y.size(); ++i)
    {
        if (i == 0 || i == y.size() - 1)
        {
            primit += y[i] / 2;
        }
        else
        {
            primit += y[i];
        }
    }

    primit *= h;
    return primit;
}

//симпосн
double simpson(double a, double b, int n, double (*func)(double)) {
    double h = (b - a) / (n);
    double primit = 0.0;
    std::vector<double> y;

    double x = a;
    while (x <= b) {
        y.push_back((*func)(x));
        x += h;
    }

    for (int i = 0; i < y.size(); ++i) {
        if (i == 0 || i == y.size() - 1) {
            primit += y[i];
        }
        else if (i % 2 == 0) {
            primit += 2 * y[i];
        }
        else {
            primit += 4 * y[i];
        }
    }

    primit *= h / 3.0;
    return primit;
}

//Ньютон 3/8
double newton(double a, double b, int n, double (*func)(double)) {
    double h = (b - a) / (n);
    double primit = 0.0;
    std::vector<double> y;

    double x = a;
    while (x <= b) {
        y.push_back((*func)(x));
        x += h;
    }

    for (int i = 0; i < y.size(); i++) {
        if (i == 0 || i == y.size() - 1) {
            primit += y[i];
        }
        //else if (i % 3 == 0) {
        //    primit += 2 * y[i];
        //}
        else {
            primit += 3 * y[i];
        }
    }

    primit *= (3 * h / 8);
    //primit -= 0.0625;
    return primit;
}


//квадратура ньютона котеса
double kotes(double a, double b, int n, double (*func)(double)) {
    std::vector<double> H;
    int N;

    if (n == 1) {
        H = { 1, 1 };
        N = 2;
    }
    else if (n == 2) {
        H = { 1, 4, 1 };
        N = 6;
    }
    else if (n == 3) {
        H = { 1, 3, 3, 1 };
        N = 8;
    }
    else if (n == 4) {
        H = { 7, 32, 12, 32, 7 };
        N = 90;
    }
    else if (n == 5) {
        H = { 19, 75, 50, 50, 75, 19 };
        N = 288;
    }
    else if (n == 6) {
        H = { 41, 216, 27, 272, 27, 216, 41 };
        N = 840;
    }
    else if (n == 7) {
        H = { 751, 3577, 1323, 2989, 2989, 1323, 3577, 751 };
        N = 17280;
    }
    else if (n == 8) {
        H = { 989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989 };
        N = 28350;
    }
    else {
        std::cerr << "Invalid value of n!" << std::endl;
        return 0.0;
    }

    double h = (b - a) / (n);
    double primit = 0.0;
    std::vector<double> y;

    double x = a;
    while (x <= b) {
        y.push_back((*func)(x));
        x += h;
    }


    for (int i = 0; i < y.size(); ++i) {
        primit += H[i] * y[i] / N;
    }


    primit *= (b - a); // *sum_y;
    
    return primit;
}


//функция полиномов лежандра
double polinom(int n, double tk) {

    if (n == 0) {
        return 1;
    }

    if (n == 1) {
        return tk;
    }

    if (n == 2) {
        double a = (1 / 2.0) * (3 * tk * tk - 1);
        return a;
    }

    if (n == 3) {
        return (1 / 2.0) * (5 * tk * tk * tk - 3 * tk);
    }

    if (n == 4) {
        return (1 / 8.0) * (35 * tk * tk * tk * tk - 30 * tk * tk + 3);
    }

}




//квадратура Гауса в узлах Лежандра
double gaus(double a, double b, int n, double (*func)(double)) {
    //double h = (b - a) / (n-1);
    //double x = a, t;

    double primit = 0.0;
    double c, f;
    vector<double> t;

    if (n == 1) {
        t = { 0 };
    }

    if (n == 2) {
        t = { -sqrt(3.0) / 3.0, sqrt(3.0) / 3.0 };
    }

    if (n == 3) {
        t = { -sqrt(15.0) / 5.0, 0, sqrt(15.0) / 5.0 };
    }

    if (n == 4) {
        t = { -sqrt(525 + 70 * sqrt(30)) / 35.0, -sqrt(525 - 70 * sqrt(30)) / 35.0, sqrt(525 + 70 * sqrt(30)) / 35.0, sqrt(525 - 70 * sqrt(30)) / 35.0 };
    }

    if (n == 5) {
        t = { -sqrt(245 + 14 * sqrt(70)) / 21.0, -sqrt(245 - 14 * sqrt(70)) / 21.0, sqrt(245 + 14 * sqrt(70)) / 21.0, sqrt(245 - 14 * sqrt(70)) / 21.0 , 0 };
    }

    for (int k = 0; k < n; k++) {
        primit += (b - a) / 2 * 2 * (1 - t[k] * t[k]) / (n * n * polinom(n - 1, t[k]) * polinom(n - 1, t[k])) * func((b + a) / 2 + (b - a) / 2 * t[k]);
    }

    //it = polinom(n - 1, t[0]);

    return primit;
}





int main() {
    setlocale(LC_ALL, "Russian");

    double h = 2;
    //double a = -1 * h / 2, b = h / 2, n = 4;
    double a = 0, b = 2, n = 4;

    cout << "Квадратура Гауса для N: " << gaus(a, b, 5, func_N) << endl;
    cout << "Квадратура Гауса для M: " << gaus(a, b, 5, func_M) << endl;
    cout << endl;

    //cout << "Составная формула трапеции при n = 3: " << trapezoid(a, b, 3, func_M) << endl;
    //cout << "Погрешность при n = 3: " << trapezoid(a, b, 3, func_M) - 2 / 3.0 << endl << endl;
    //
    //cout << "Составная формула трапеции при n = 6: " << trapezoid(a, b, 6, func_M) << endl;
    //cout << "Погрешность при n = 6: " << trapezoid(a, b, 6, func_M) - 2 / 3.0 << endl << endl;
    //
    //
    //cout << "Составная формула Симпсона при n = 3: " << simpson(a, b, 3, func_N) << endl;
    //cout << "Погрешность при n = 3: " << simpson(a, b, 3, func_N) - 2 / 5.0 << endl << endl;
    //
    //cout << "Составная формула Симпсона при n = 6: " << simpson(a, b, 6, func_N) << endl;
    //cout << "Погрешность при n = 6: " <<  simpson(a, b, 6, func_N) - 2 / 5.0 << endl << endl;

    cout << "Формула Ньютона 3/8 при n = 3: " << newton(a, b, 3, func_N) << endl;
    cout << "Погрешность при n = 3: " << newton(a, b, 3, func_N) - 2 / 5.0 << endl << endl;

    cout << "Формула Ньютона 3/8 при n = 6: " << newton(a, b, 6, func_N) << endl;
    cout << "Погрешность при n = 6: " << newton(a, b, 6, func_N) - 2 / 5.0 << endl << endl;

    cout << "Составная формула трапеции для N: " << trapezoid(a, b, n, func_N) << endl;
    cout << "Составная формула трапеции для M: " << trapezoid(a, b, n, func_M) << endl;
    cout << endl;

    cout << "Составная формула Симпсона для N: " << simpson(a, b, n, func_N) << endl;
    cout << "Составная формула Симпсона для M: " << simpson(a, b, n, func_M) << endl;
    cout << endl;

    cout << "Формула Ньютона 3/8 для N: " << newton(a, b, n, func_N) << endl;
    cout << "Формула Ньютона 3/8 для M: " << newton(a, b, n, func_M) << endl;
    cout << endl;

    if (n >= 1.0 && n <= 8.0) {
        cout << "Квадратура Ньютона-Котеса для N: " << kotes(a, b, n, func_N) << endl;
        cout << "Квадратура Ньютона-Котеса для M: " << kotes(a, b, n, func_M) << endl;
    }




    return 0;
}




















    //cout << "Составная формула трапеции при n = 9: " << trapezoid(a, b, 9, func_N) << endl;
    //cout << "Погрешность при n = 9: " << trapezoid(a, b, 9, func_M) - 2 / 3.0 << endl << endl;





//cout << "Составная формула трапеции при n = 3: " << trapezoid(a, b, 3, func_N) << endl;
//cout << "Погрешность при n = 3: " << trapezoid(a, b, 3, func_N) - 2 / 5.0 << endl << endl;
//
//cout << "Составная формула трапеции при n = 5: " << trapezoid(a, b, 5, func_N) << endl;
//cout << "Погрешность при n = 5: " << trapezoid(a, b, 5, func_N) - 2 / 5.0 << endl << endl;
//
//cout << "Составная формула трапеции при n = 9: " << trapezoid(a, b, 9, func_N) << endl;
//cout << "Погрешность при n = 9: " << trapezoid(a, b, 9, func_N) - 2 / 5.0 << endl << endl;


/*

double ck(int n, double tk) {
    if (n == 1) {
        return (2*(1 - pow(tk, 2))) / pow(n, 2);
    }

    if (n == 2) {
        return (2 * (1 - pow(tk, 2))) / (pow(n, 2) * pow(tk, 2));
    }

    if (n == 3) {
        return (2 * (1 - pow(tk, 2))) / (pow(n, 2)* pow((1 / 2.0) * (3 * pow(tk, 2) - 1), 2));
    }

    if (n == 4) {
        return (2 * (1 - pow(tk, 2))) / (pow(n, 2) * pow((1 / 2.0) * (5 * pow(tk, 3) - 3 * tk), 2));
    }

    if (n == 5) {
        //return (2.0 * (1.0 - pow(tk, 2.0))) / (pow(n, 2.0) * pow((1.0 / 8.0) * (35.0 * pow(tk, 5) - 70.0 * pow(tk, 3.0) + 15.0 * tk), 2.0));
        double c = 35 * pow(tk, 5) - 70;
        double d = pow(tk, 3);

        double a = (1.0 / 8.0) * (c * d + 15 * tk);
        double b = (pow(n, 2) * pow((1.0 / 8.0) * (35 * pow(tk, 5) - 70 * pow(tk, 3) + 15 * tk), 2));
        return (2 * (1 - pow(tk, 2))) / b;
    }

}

//квадратура Гауса в узлах Лежандра
double gaus(double a, double b, int n, double (*func)(double)) {
    double h = (b - a) / (n-1);
    double x = a, t;
    double it = 0.0;
    double c, f;


    for (int k = 1; k <= n; k++) {
        x = x + h;
        t = ((2 * x) - b - a) / (b - a);
        c = ck(n, t);
        f = (*func)((b + a) / 2 + t * (b - a) / 2);
        it += c + f;
    }

    it *= (b-a)/2;

    return it;
}



double newton(double a, double b, int n, double (*func)(double)) {
    double h = (b - a) / n;
    double it = 0.0;
    std::vector<double> y;

    double x = a;
    while (x <= b) {
        y.push_back((*func)(x));
        x += h;
    }

    for (size_t i = 0; i < y.size(); ++i) {
        if (i % 4 == 0 || i % 4 == 3) {
            it += y[i];
        }
        else {
            it += 3 * y[i];
        }
    }

    it *= (3 * h / 8);
    return it;
}
*/