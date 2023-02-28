#include <iostream>
#include <vector>
#include <fstream>
#include<math.h>
using namespace std;

//PENDOLO TRIPLO RUNGE-KUTTA


//Definisco costanti

double g=9.806;

double l1=1.;
double l2=1.;
double l3=1.;

double m1=1.;
double m2=1.;
double m3=1.;

double At;

//Definisco vettore y 

vector<double> yi(6,0.);


//Definisco funzioni per equazioni del moto

double funz1 ( double t1, double t2, double t3, double dt1, double dt2, double dt3){
    return -(g*pow(m2,2)*sin(t1) + g*pow(m2,2)*sin(t1 - 2.*t2) + 2.*pow(dt2,2)*l2*pow(m2,2)*sin(t1 - t2) + 2.*g*m1*m2*sin(t1) + 
           g*m1*m3*sin(t1) + g*m2*m3*sin(t1) - (g*m1*m3*sin(t1 - 2.*t2 + 2*t3))/2. - (g*m1*m3*sin(t1 + 2*t2 - 2.*t3))/2.
           + g*m2*m3*sin(t1 - 2.*t2) + pow(dt1,2)*l1*pow(m2,2)*sin(2*t1 - 2.*t2) + 2.*pow(dt2,2)*l2*m2*m3*sin(t1 - t2) + 
           pow(dt3,2)*l3*m2*m3*sin(t1 - t3) + pow(dt1,2)*l1*m2*m3*sin(2*t1 - 2.*t2) + pow(dt3,2)*l3*m2*m3*sin(t1 - 2*t2 + t3))/(l1*(2.*m1*m2 + m1*m3 + m2*m3 - pow(m2,2)*cos(2*t1 - 2*t2) + pow(m2,2) - m2*m3*cos(2*t1 - 2*t2) - m1*m3*cos(2.*t2 - 2.*t3)));
}

double funz2 (double t1, double t2, double t3, double dt1, double dt2, double dt3){
    return (g*pow(m2,2)*sin(2.*t1 - t2) - g*pow(m2,2)*sin(t2) + g*m1*m2*sin(2.*t1 - t2) + (g*m1*m3*sin(2.*t1 - t2))/2. + 
           g*m2*m3*sin(2.*t1 - t2) + 2.*pow(dt1,2)*l1*pow(m2,2)*sin(t1 - t2) - g*m1*m2*sin(t2) - (g*m1*m3*sin(t2))/2. - g*m2*m3*sin(t2) - 
           (g*m1*m3*sin(2.*t1 + t2 - 2.*t3))/2. - (g*m1*m3*sin(t2 - 2*t3))/2. + pow(dt2,2)*l2*pow(m2,2)*sin(2.*t1 - 2.*t2) + 
           2.*pow(dt1,2)*l1*m1*m2*sin(t1 - t2) + pow(dt1,2)*l1*m1*m3*sin(t1 - t2) + 2.*pow(dt1,2)*l1*m2*m3*sin(t1 - t2) - 2.*pow(dt3,2)*l3*m1*m3*sin(t2 - t3) - 
           pow(dt3,2)*l3*m2*m3*sin(t2 - t3) + pow(dt3,2)*l3*m2*m3*sin(2*t1 - t2 - t3) + pow(dt2,2)*l2*m2*m3*sin(2.*t1 - 2.*t2) -
           pow(dt2,2)*l2*m1*m3*sin(2.*t2 - 2.*t3) - pow(dt1,2)*l1*m1*m3*sin(t1 + t2 - 2*t3))/(l2*(2*m1*m2 + m1*m3 + m2*m3 - pow(m2,2)*cos(2.*t1 - 2.*t2) +
           pow(m2,2) - m2*m3*cos(2.*t1 - 2.*t2) - m1*m3*cos(2.*t2 - 2.*t3)));
}

double funz3 (double t1, double t2, double t3, double dt1, double dt2, double dt3){
    return (m1*(g*m2*sin(2.*t1 - t3) - g*m3*sin(t3) - g*m2*sin(2.*t1 - 2.*t2 + t3) - g*m3*sin(2.*t1 - 2.*t2 + t3) - g*m2*sin(t3) + g*m2*sin(2.*t2 - t3) + g*m3*sin(2.*t1 - t3) + 
    g*m3*sin(2.*t2 - t3) + 2.*pow(dt1,2)*l1*m2*sin(t1 - t3) + 2.*pow(dt1,2)*l1*m3*sin(t1 - t3) + 4.*pow(dt2,2)*l2*m2*sin(t2 - t3) + 4.*pow(dt2,2)*l2*m3*sin(t2 - t3) + 2.*pow(dt3,2)*l3*m3*sin(2.*t2 - 2.*t3) - 
    2.*pow(dt1,2)*l1*m2*sin(t1 - 2.*t2 + t3) - 2.*pow(dt1,2)*l1*m3*sin(t1 - 2.*t2 + t3)))/(2.*l3*(2.*m1*m2 + m1*m3 + m2*m3 - pow(m2,2)*cos(2.*t1 - 2.*t2) + pow(m2,2) - m2*m3*cos(2.*t1 - 2.*t2) - m1*m3*cos(2.*t2 - 2.*t3)));
}



//Definisco funzione f e funzioni per calcolo vettoriale

vector<double> f(vector<double> y) {
	vector<double> funz{y[3],y[4],y[5],funz1(y[0],y[1],y[2],y[3],y[4],y[5]),funz2(y[0],y[1],y[2],y[3],y[4],y[5]),funz3(y[0],y[1],y[2],y[3],y[4],y[5])};
    return funz;
}

vector<double> sum(vector<double> a,vector<double> b) {
	vector<double> somma(6,0.);
    for (int i=0; i<6; i++){
        somma[i]=a[i]+b[i];
    }
    return somma;
}

vector<double> pps(vector<double> a, double b) {
	vector<double> prodotto(6,0.);
    for (int j=0; j<6; j++){
        prodotto[j]=a[j]*b;
    }
    return prodotto;
}


// Definisco funzioni per algoritmo Runge-kutta

vector<double> Y1(vector<double> y) {
	vector<double> y1(6,0.);
    for (int q=0; q<6; q++){
        y1[q]=y[q];
    }
    return y1;
}


vector<double> Y2(vector<double> y) {
	vector<double> y2=sum(y,pps(f(Y1(y)),At/2.));
    return y2;

}

vector<double> Y3(vector<double> y) {
	vector<double> y3=sum(y,pps(f(Y2(y)),At/2.));
    return y3;

}

vector<double> Y4(vector<double> y) {
	vector<double> y4=sum(y,pps(f(Y3(y)),At));
    return y4;

}

//Definisco funzioni energia

double Ek ( double t1, double t2, double t3, double dt1, double dt2, double dt3){
    double Ec=0.5*pow(l1,2)*pow(dt1,2)*(m1+m2+m3)+0.5*pow(l2,2)*pow(dt2,2)*(m2+m3)+0.5*pow(l3,2)*pow(dt3,2)+
    pow(l2,2)*dt1*dt2*cos(t1-t2)*(m2+m3)+l2*l3*dt2*dt3*cos(t2-t3)*m3+l1*l3*dt1*dt3*cos(t1-t3)*m3;
    return Ec;
}

double Ep ( double t1, double t2, double t3){
    double Epot= g*(l1+l2+l3)*(m1+m2+m3)-g*(l1*(m1+m2+m3)*cos(t1)+l2*(m2+m3)*cos(t2)+l3*m3*cos(t3));
    return Epot;   
}
 

//Funzione main

int main() {

    double t10=0.;
    double t20=0.;
    double t30=0.;
    double dt10=0.;
    double dt20=0.;
    double dt30=0.;
    double N;
 

    
    cout << "Il programma computa le posizioni di un pendolo triplo in un intervallo di 10 s" << endl;
    cout << "Numero di punti:";
    cin >> N;

    At=10./N;
    
    cout << "Valore Theta1 iniziale [rad]:";
    cin >> t10;
    cout << "Valore Theta2 iniziale [rad]:";
    cin >> t20;
    cout << "Valore Theta3 iniziale [rad]:";
    cin >> t30;
    cout << "Valore velocità angolare iniziale 1 [rad/s]:";
    cin >> dt10;
    cout << "Valore velocità angolare iniziale 2 [rad/s]:";
    cin >> dt20;
    cout << "Valore velocità angolare iniziale 3 [rad/s]:";
    cin >> dt30;

    yi[0]=t10;
    yi[1]=t20;
    yi[2]=t30;
    yi[3]=dt10;
    yi[4]=dt20;
    yi[5]=dt30;

    ofstream testo("pendolo_triplo_angoli.txt");
    ofstream testo1("pendolo_triplo.txt");
    ofstream testo2("Energia.txt");





    testo << 0 << ' ' << yi[0] << ' ' << yi[1] << ' ' << yi[2] << endl;
    testo1 << 0 << ' ' << 0 << ' ' <<l1*sin(yi[0]) << ' ' << -l1*cos(yi[0]) << ' ' << l1*sin(yi[0])+l2*sin(yi[1]) << ' ' << -(l1*cos(yi[0])+l2*cos(yi[1])) << ' ' << l1*sin(yi[0])+l2*sin(yi[1])+l3*sin(yi[2]) << ' ' << -(l1*cos(yi[0])+l2*cos(yi[1])+l3*cos(yi[2])) << endl << endl;
    testo2 << 0  << ' ' << Ek(yi[0],yi[1],yi[2],yi[3],yi[4],yi[5])+Ep(yi[0],yi[1],yi[2]) << endl;
    
    //Calcolo posizioni e velocità e energia n+1 salvando risultati su file di testo

    for (int w=1; w<=N; w++){

        yi=sum(yi,pps(sum(sum(f(Y1(yi)),pps(f(Y2(yi)),2)),sum(pps(f(Y3(yi)),2),f(Y4(yi)))),At/6.));
        testo << At*w << ' ' << yi[0] << ' ' << yi[1] << ' ' << yi[2] << endl;
        testo1 << 0 << ' ' << 0 << ' ' << l1*sin(yi[0]) << ' ' << -l1*cos(yi[0]) << ' ' << l1*sin(yi[0])+l2*sin(yi[1]) << ' ' << -(l1*cos(yi[0])+l2*cos(yi[1])) << ' ' << l1*sin(yi[0])+l2*sin(yi[1])+l3*sin(yi[2]) << ' ' << -(l1*cos(yi[0])+l2*cos(yi[1])+l3*cos(yi[2])) << endl << endl;
        testo2 << At*w << ' ' << Ek(yi[0],yi[1],yi[2],yi[3],yi[4],yi[5])+Ep(yi[0],yi[1],yi[2]) << endl;
    }

}