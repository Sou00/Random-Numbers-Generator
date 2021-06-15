//Igor Łonak
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

int seria1[19][19]={
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 3, 3, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 2, 3, 3, 4, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 2, 3, 4, 4, 5, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 2, 3, 3, 4, 5, 5, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0},
        {2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 9, 0, 0, 0, 0, 0, 0, 0},
        {2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 9, 10, 0, 0, 0, 0, 0, 0},
        {2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 9, 10, 10, 11, 0, 0, 0, 0, 0},
        {2, 3, 4, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 11, 11, 0, 0, 0, 0},
        {2, 3, 4, 5, 6, 7, 7, 8, 9, 9, 10, 10, 11, 11, 12, 12, 0, 0, 0},
        {2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 10, 11, 11, 12, 12, 13, 13, 0, 0},
        {2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13, 13, 14, 14, 0},
        {2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15}
};

int seria2[19][19]={
        {4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 8, 9, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 10, 11, 12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 10, 11, 12, 13, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 10, 11, 12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 12, 13, 14, 14, 15, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 12, 13, 14, 15, 16, 16, 17, 0, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 12, 13, 14, 15, 16, 17, 17, 18, 0, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 12, 13, 15, 16, 16, 17, 18, 19, 19, 0, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 13, 14, 15, 16, 17, 18, 18, 19, 20, 20, 0, 0, 0, 0, 0},
        {5, 7, 9, 11, 13, 14, 15, 16, 17, 18, 19, 20, 20, 21, 22, 0, 0, 0, 0},
        {5, 7, 9, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21, 22, 23, 0, 0, 0},
        {5, 7, 9, 11, 13, 14, 15, 17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 0, 0},
        {5, 7, 9, 11, 13, 14, 15, 17, 18, 19, 20, 21, 22, 22, 23, 24, 24, 25, 0},
        {5, 7, 9, 11, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 24, 25, 26, 26}
};



int factorial(int x){

    int y=1;
    for(int i=1;i<=x;i++){
        y*=i;
    }
    return y;
}

void printI(int* tab,int amount){
    for (int i = 0; i < amount; ++i) {
        cout<<tab[i]<<", ";
    }
    cout<<endl;
}
void printI(unsigned int* tab,int amount){
    for (int i = 0; i < amount; ++i) {
        cout<<tab[i]<<", ";
    }
    cout<<endl;
}

void printD(double* tab,int amount){
    cout<<fixed<<setprecision(6);
    for (int i = 0; i < amount; ++i) {
        cout<<tab[i]<<", ";
    }
    cout<<endl;
}

unsigned int* G(int a, int mod, int amount, unsigned int x){
    unsigned int *g=new unsigned int[amount];
    g[0]= x;
    for (int i = 1; i < amount; ++i) {
        g[i]=a*g[i-1]%mod;
    }

    return g;
}

double* J(int a, int mod, int amount, unsigned int x){
    unsigned int * g= G(a,mod,amount,x);
    double *j= new double[amount];

    for (int i = 0; i < amount; ++i) {
        j[i]=(double)g[i]/mod;
    }

    delete g;

    return j;
}

int* B(int a,int mod,int amount, double p, unsigned int x){
    double * j= J(a,mod,amount,x);
    int* b = new int[amount];

    for (int i = 0; i < amount; ++i) {
        if(j[i]<=p)
            b[i]=1;
        else
            b[i]=0;
    }
    delete j;

    return b;
}

int* D(int a, int mod, int amount, double p, unsigned int x){
    int* b;
    int* d= new int[amount];

    for (int i = 0; i < amount; ++i) {
        b= B(a,mod,11,p,x* 7*(i+1));
        d[i]=0;
        for (int j=0 ; j < 10; ++j) {
            d[i]+=b[j+1];
        }
        delete b;
    }

    return d;
}

int* P(int a,int mod,int amount, unsigned int x,double lambda){
    double L= exp(-lambda);
    int k=0;
    double p=1.0;
    int temp=0;
    int* pois=new int[amount];
    double * j= J(a,mod,amount, x);

    for (int i = 0; i < amount; ++i) {
        p=1.0;
        k=0;
        while(p>L){
            k++;
            p*=j[temp++];
            if(temp>=amount) {
                delete j;
                j = J(a, mod, amount, x*7*(i+1));
                temp=0;
            }
        }
        pois[i]=k-1;
    }

    delete j;

    return pois;

}

double* W (int a, int mod, int amount, unsigned int x){

    double * j= J(a,mod,amount,x);
    double* w = new double [amount];

    for (int i = 0; i < amount; ++i) {
        w[i]= -log(1-j[i]);
    }

    delete j;

    return w;
}

double* N (int a,int mod, int amount, unsigned int x){
    double * j= J(a,mod,amount,x);
    double* z=new double [amount];

    for (int i = 0; i < amount; i++) {
        double r= -2*log(j[i]);
        double omega=2*M_PI*j[i+1];
        z[i]= sqrt(r)* cos(omega);
        i++;
        if(i<amount)
            z[i]= sqrt(r)* sin(omega);
    }

    delete j;

    return z;


}

void insertionSort(int arr[], int n)
{
    int i, key, j;
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void insertionSort(double arr[], int n)
{
    int i, j;
    double key;
    for (i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

double mediana(unsigned int* tab,int amount){
    double m=0;
    int tab1[amount];
    for (int i = 0; i < amount; ++i) {
        tab1[i]=tab[i];
    }
    insertionSort(tab1,amount);
    if(amount%2==0)
        m=tab1[(amount-1)/2];
    else
        m=(tab1[(amount-2)/2]+tab1[(amount)/2])/2;

    return m;
}

double mediana(double* tab, int amount){
    double m=0;
    double tab1[amount];
    for (int i = 0; i < amount; ++i) {
        tab1[i]=tab[i];
    }
    insertionSort(tab1,amount);
    if(amount%2==0)
        m=tab1[(amount-1)/2];
    else
        m=(tab1[(amount-2)/2]+tab1[(amount)/2])/2;

    return m;
}

void seria(unsigned int* tab,int amount){
    double m=mediana(tab,amount);
    int u=1;
    char tab1[amount];
    int n1=0;
    int n2=0;

    for (int i = 0; i < amount; ++i) {
        if(tab[i]>m) {
            tab1[i] = 'a';
            n1++;
        }
        else if(tab[i]<m) {
            tab1[i] = 'b';
            n2++;
        }
        else
            tab1[i]='-';
    }

    for (int i = 0; i < amount; ++i) {
        cout<<tab1[i]<<" ";
    }
    cout <<endl;

    for (int i = 1; i < amount; ++i) {
        if(tab1[i-1]!=tab1[i] && tab1[i-1]!='-')
            u++;
    }
    cout<<"ilość podciągów: "<< u<<" mediana: "<<m <<endl;
    cout<<"liczba a: "<<n1 <<" liczba b:"<<n2<<endl;


    int k1=seria1[n1-2][n2-2];
    int k2=seria2[n1-2][n2-2];

    cout<<"Akceptowalny przedział: ("<<k1<<";"<<k2<<")"<<endl;

    if(u>k1 && u<k2)
        cout<<"Jest ok"<<endl;
    else
        cout<<"Jest nieok"<<endl;

    cout<<"---------------------------------------------------------------------------------"<<endl;

}

void seria(double* tab,int amount){
    double m=mediana(tab,amount);
    int u=1;
    char tab1[amount];
    int n1=0;
    int n2=0;

    for (int i = 0; i < amount; ++i) {
        if(tab[i]>m) {
            tab1[i] = 'a';
            n1++;
        }
        else if(tab[i]<m) {
            tab1[i] = 'b';
            n2++;
        }
        else
            tab1[i]='-';
    }

    for (int i = 0; i < amount; ++i) {
        cout<<tab1[i]<<" ";
    }
    cout <<endl;

    for (int i = 1; i < amount; ++i) {
        if(tab1[i-1]!=tab1[i] && tab1[i-1]!='-')
            u++;
    }
    cout<<"ilość podciągów: "<< u<<" mediana: "<<m <<endl;
    cout<<"liczba a: "<<n1 <<" liczba b:"<<n2<<endl;


    int k1=seria1[n1-2][n2-2];
    int k2=seria2[n1-2][n2-2];

    cout<<"Akceptowalny przedział: ("<<k1<<";"<<k2<<")"<<endl;


    if(u>k1 && u<k2)
        cout<<"Jest ok"<<endl;
    else
        cout<<"Jest nieok"<<endl;

    cout<<"---------------------------------------------------------------------------------"<<endl;


}

double critical[10]={3.841, 5.991, 7.815, 9.488, 11.071, 12.592, 14.067, 15.507, 16.919, 18.307};

void chiB(int* b,int amount, double p){
    double obs[2];// 0  1
    obs[0]=0.0;
    obs[1]=0.0;
    for (int i = 0; i < amount; ++i) {
        if(b[i]==0){
            obs[0]++;
        }
        else
            obs[1]++;
    }
    double exp[2];
    exp[0]=amount*(1-p);
    exp[1]=amount*p;

    double chi=pow(obs[0]-exp[0],2)/exp[0]+pow(obs[1]-exp[1],2)/exp[1];

    cout<<"Observed:"<<endl;
    cout<<"0 1"<<endl;

    for (int i = 0; i < 2; ++i) {
        cout<<obs[i]<<" ";
    }

    cout<<endl;

    cout<<"Expected:"<<endl;
    cout<<"0 1"<<endl;

    for (int i = 0; i < 2; ++i) {
        cout<<exp[i]<<" ";
    }
    cout<<endl;

    cout<<"Chi-kwadrat: "<<chi<<endl;
    if(chi<critical[0])
        cout<<"Jest OK"<<endl;
    else
        cout<<"Nie jest OK"<<endl;


    cout<<"---------------------------------------------------------------------------------"<<endl;

}

void chiD (int* d,int amount, double p){
    double obs[4]; //0-2 3-5 6-8 reszta
    for (int i = 0; i < 4; ++i) {
        obs[i]=0.0;
    }

    for (int i = 0; i < amount; ++i) {
        if(d[i]<=2)
            obs[0]++;
        else if(d[i]<=5)
            obs[1]++;
        else if(d[i]<=8)
            obs[2]++;
        else obs[3]++;
    }
    double exp[4];

    double n= factorial(10);



    for (int i = 6; i < 9; ++i) {
        exp[2]+=amount*(n/(factorial(i)*factorial(10-i))* pow(p, i)*pow(1-p,10-i));
    }
    for (int i = 3; i < 6; ++i) {
        exp[1]+=amount*(n/(factorial(i)*factorial(10-i))* pow(p, i)*pow(1-p,10-i));
    }
    for (int i = 0; i < 3; ++i) {
        exp[0]+=amount*(n/(factorial(i)*factorial(10-i))* pow(p, i)*pow(1-p,10-i));
    }

    exp[3]=amount-(exp[0]+exp[1]+exp[2]);

    cout<<"Observed:"<<endl;
    cout<<"0-2           3-5        6-8        reszta"<<endl;
    for (int i = 0; i < 4; ++i) {
        cout<<obs[i]<<" ";
    }

    cout<<endl;

    cout<<"Expected:"<<endl;
    cout<<"0-2           3-5        6-8        reszta"<<endl;

    for (int i = 0; i < 4; ++i) {
        cout<<exp[i]<<" ";
    }
    cout<<endl;

    double chi=0.0;

    for (int i = 0; i < 4; ++i) {
        chi+=pow(obs[i]-exp[i],2)/exp[i];
    }

    cout<<"Chi-kwadrat: "<<chi<<endl;
    if(chi<critical[2])
        cout<<"Jest OK"<<endl;
    else
        cout<<"Nie jest OK"<<endl;

    cout<<"---------------------------------------------------------------------------------"<<endl;

}

void chiP (int* p, int amount, double lambda){
    double obs[4];//0-2 3-5 6-8 reszta
    for (int i = 0; i < 4; ++i) {
        obs[i]=0.0;
    }

    for (int i = 0; i < amount; ++i) {
        if(p[i]<=2)
            obs[0]++;
        else if(p[i]<=5)
            obs[1]++;
        else if(p[i]<=8)
            obs[2]++;
        else obs[3]++;
    }

    double exp[4];

    for (int i = 0; i < 3; ++i) {
        exp[0]+=amount* pow(lambda,i)*std::exp(-lambda)/ factorial(i);
    }
    for (int i = 3; i < 6; ++i) {
        exp[1]+=amount* pow(lambda,i)*std::exp(-lambda)/ factorial(i);
    }
    for (int i = 6; i < 9; ++i) {
        exp[2]+=amount* pow(lambda,i)*std::exp(-lambda)/ factorial(i);
    }

    exp[3]=amount-(exp[0]+exp[1]+exp[2]);

    cout<<"Observed:"<<endl;
    cout<<"0-2           3-5        6-8        reszta"<<endl;

    for (int i = 0; i < 4; ++i) {
        cout<<obs[i]<<" ";
    }

    cout<<endl;

    cout<<"Expected:"<<endl;
    cout<<"0-2           3-5        6-8        reszta"<<endl;

    for (int i = 0; i < 4; ++i) {
        cout<<exp[i]<<" ";
    }
    cout<<endl;

    double chi=0.0;

    for (int i = 0; i < 4; ++i) {
        chi+=pow(obs[i]-exp[i],2)/exp[i];
    }

    cout<<"Chi-kwadrat: "<<chi<<endl;
    if(chi<critical[2])
        cout<<"Jest OK"<<endl;
    else
        cout<<"Nie jest OK"<<endl;

    cout<<"---------------------------------------------------------------------------------"<<endl;

}

int main() {

    int mod =pow(2,31)-1;
    int a=pow(7,5);
    int amount=10000;
    int x= pow(2,25)-1;

    //Generator G
    unsigned int* g= G(a,mod,amount,x);
    cout<<"Generator G:"<<endl;
    printI(g,amount);

    //Generator J
    double* j = J(a,mod,amount,x);
    cout<<"Generator J:"<<endl;
    printD(j,amount);

    //Bernoulli
    double p=0.3;
    int* b = B(a,mod,amount,p,x);
    cout<<"BERNOULLI: "<<endl;
    printI(b,amount);

    //Dwumianowy
    int* d= D(a,mod,amount,p,x);
    cout<<"DWUMIANOWY: "<<endl;
    printI(d,amount);

    //Poisson
    double lambda=3;
    int*pois= P(a,mod,amount,x,lambda);
    cout<<"POISSON: "<<endl;
    printI(pois,amount);

    //Wykładniczy
    double *w= W(a,mod,amount,x);
    cout<<"Wykładniczy:"<<endl;
    printD(w,amount);

    //Normalny
    double* n= N(a,mod,amount,x);
    cout<<"Normalny: "<<endl;
    printD(n,amount);

    cout<<"---------------------------------------------------------------------------------"<<endl;
    //Testy serii
    cout<<"Test serii dla G(n=5):"<<endl;
    seria(g,5);
    cout<<"Test serii dla G(n=10):"<<endl;
    seria(g,10);
    cout<<"Test serii dla G(n=15):"<<endl;
    seria(g,15);
    cout<<"Test serii dla G(n=20):"<<endl;
    seria(g,20);

    cout<<"Test serii dla J(n=5):"<<endl;
    seria(j,5);
    cout<<"Test serii dla J(n=10):"<<endl;
    seria(j,10);
    cout<<"Test serii dla J(n=15):"<<endl;
    seria(j,15);
    cout<<"Test serii dla J(n=20):"<<endl;
    seria(j,20);

    cout<<"---------------------------------------------------------------------------------"<<endl;
    //Testy chi-kwadrat
    cout<<"Test chi-kwadrat dla Bernoulliego(n=500, p=0.3):"<<endl;
    chiB(b,500,p);
    cout<<"Test chi-kwadrat dla Bernoulliego(n=1000, p=0.3):"<<endl;
    chiB(b,1000,p);
    cout<<"Test chi-kwadrat dla Bernoulliego(n=5000, p=0.3):"<<endl;
    chiB(b,5000,p);
    cout<<"Test chi-kwadrat dla Bernoulliego(n=10000, p=0.3):"<<endl;
    chiB(b,10000,p);

    cout<<"Test chi-kwadrat dla Dwumianowego(n=500,k=10,p=0.3):"<<endl;
    chiD(d,500,p);
    cout<<"Test chi-kwadrat dla Dwumianowego(n=1000,k=10,p=0.3):"<<endl;
    chiD(d,1000,p);
    cout<<"Test chi-kwadrat dla Dwumianowego(n=5000,k=10,p=0.3):"<<endl;
    chiD(d,5000,p);
    cout<<"Test chi-kwadrat dla Dwumianowego(n=10000,k=10,p=0.3):"<<endl;
    chiD(d,10000,p);

    cout<<"Test chi-kwadrat dla Poisson(n=500,lambda=3):"<<endl;
    chiP(pois,500,lambda);
    cout<<"Test chi-kwadrat dla Poisson(n=1000,lambda=3):"<<endl;
    chiP(pois,1000,lambda);
    cout<<"Test chi-kwadrat dla Poisson(n=5000,lambda=3):"<<endl;
    chiP(pois,5000,lambda);
    cout<<"Test chi-kwadrat dla Poisson(n=10000,lambda=3):"<<endl;
    chiP(pois,10000,lambda);


    return 0;
}
