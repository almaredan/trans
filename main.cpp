#include <cstdlib>
#include <iostream>

#include <math.h>
#include <vector>

using namespace std;

double integrate(double(*f)(double), double x0, double h) {
    return f(x0)*h  + ( f(x0+h/2) + f(x0-h/2) - 2*f(x0) )*h/6;
}

double f_x(double x) {
    return x*x*x;
}


double u_big(vector<double> & u, int i, double h) { 					//Метод, считающий большие величины с применением МНК
    double u_tmp;
    double k;
    double b;
    int n = u.size();
    if ((i == 0) || (i == n)) {
        if (i == 0) {
            k = (-3*u[0] + u[1] + 2*u[2])/(5*h);
            b = u[0];
        } else {
            k = (3*u[n] - u[n-1] - 2*u[n-2])/(5*h);
            b = u[n-2];
        }
    } else {
        k = (u[i+2] - u[i])/(2*h);
        b = u[i];
    }
    u_tmp = k*(i*h - (i-1/2)*h) + b;
    return u_tmp;
}

vector<double> scheme(vector<double> & u, int N, double a,
        double x1, double x2) {
    double n = 10.;
    double Cu = 0.5;
    double h = (x2 - x1) / N;
    double tau = 0.01;
    double u_left, u_right;

    while (n > 1e-15) {
        vector<double> u_buff;
        u_buff.resize( N );
        u_buff[0] = 1.;
        u_left = u_big(u, 0, h);
        for (int x = 1; x < N; ++x) {
            u_right = u_big(u, x, h);
            u_buff[x] = u[x] - a*tau/h * (u_right - u_left)
                    + tau/h*integrate( [](double y){ return f_x(y); }, x*h + x1, h);
                    u_left = u_right;
        }
        double e = .0;
        for (int i = 1; i < N; i++) {
            e = (abs((u_buff[i] - u[i])/u[i]) > e) ? abs((u_buff[i] - u[i])/u[i]) : e;
        }
        n = (n > e) ? e : n;
        u = u_buff;
    }

    return u;
}

vector<double>  discreting(int N, double x1, double x2) {
    double h = (x2 - x1) / N;
    vector<double> x_i;
    for (double x = x1; x <= x2; x += h) {
        x_i.push_back(x);
    }

    return x_i;
}

int main (int argc, char** argv) {
    double a = 1.;
    double x1 = 0;
    double x2 = 2;

    int N_r = 100;
    int N_m = 200;

    vector<double> u;
    vector<double> v;

    u.push_back(1.);
    for (int i = 1; i < N_r; i++) {
        u.push_back(2.);
    }

    v.push_back(1.);
    for (int i = 1; i < N_m; i++) {
        v.push_back(2.);
    }

    scheme(u, N_r, a, x1, x2);
    scheme(v, N_m, a, x1, x2);

    cout << "Грубое " << endl;

    double h = (x2 - x1) / N_r;
    for (int x = 0; x < N_r; x++) {
        cout << x*h << " " << u[x] << endl;
    }

    double h1 = (x2 - x1) / N_m;
    cout << "\nТочное" << endl;

    for (int x = 0; x < N_m; x++) {
        cout << x*h1 << " " << v[x] << endl;
    }

    double u_r = integrate([](double y){ return f_x(y); }, x1, 0.01) + 1;
    vector<double> x_ir = discreting(N_r, x1, x2);
    vector<double> x_im = discreting(N_m, x1, x2);
    
    double i1 = 0, i2 = 0;
    
    //Не забудь заменить здесь и ниже константу рядом с f_x(), если она изменится
    //Получить её из интегрирования
    double err_r = abs((u[0] - u_r)/u[0]);
    int n = u.size();
    for (int i = 1; i <  n; i++) {
        double err_tmp = abs(u[i] - f_x(x_ir[i])*x_ir[i]/4 - 1)/abs(u[i]);
        err_r = (err_r > err_tmp) ? err_r : err_tmp;
    }
    
    double err_m = abs((v[0] - u_r)/v[0]);
    n = v.size();
    for (int i = 1; i <  n; i++) {
        double err_tmp = abs(v[i] - f_x(x_im[i])*x_im[i]/4 - 1)/abs(v[i]);
        err_m = (err_m > err_tmp) ? err_m : err_tmp;
    }

    double u_ex = 1.25;
    double p = logf((u[49] - u_ex)/(v[99]-u_ex)) / logf(2);

    cout << "\n Порядок аппроксимации" << " " << p << endl;
    cout << "Error: " << err_r <<" " << err_m << endl;

    return 0;
}
