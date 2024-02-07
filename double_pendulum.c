#include <stdio.h>
#include <math.h>
#define g 9.81
#define PI 3.141592
#define N 4
#define n 2

double h1 (double phi, double psi, double Pphi, double Ppsi)
{
    return (Pphi*Ppsi*sin(phi-psi))/(1+pow(sin(phi-psi), 2));
}

double h2 (double phi, double psi, double Pphi, double Ppsi)
{
    return (Pphi*Pphi + 2*Ppsi*Ppsi - 2*Pphi*Ppsi*cos(phi-psi))/(2*pow( 1 + pow(sin(phi - psi) ,2),2));
}

double func_phi (double phi, double psi, double Pphi, double Ppsi)
{
    return (Pphi - Ppsi*cos(phi-psi))/(1+pow(sin(phi-psi), 2));
}

double func_psi (double phi, double psi, double Pphi, double Ppsi)
{
    return (-1.*Pphi*cos(phi-psi) + 2*Ppsi)/(1+pow(sin(phi-psi), 2));
}

double func_Pphi(double phi, double psi, double Pphi, double Ppsi, double h1, double h2)
{
    return -1.*2*g*sin(phi)-h1+ h2*sin(2*(phi-psi));
}

double func_Ppsi(double phi, double psi, double Pphi, double Ppsi, double h1, double h2)
{
    return -1.*g*sin(psi)+h1- h2*sin(2*(phi-psi));
}

double lyapunov (double psi[n], double psiPunto_mapa[n], double delta0, double h)
{
    double sol;

    sol = sqrt(pow(psi[1]-psi[0], 2)+ pow(h*(psiPunto_mapa[1]-psiPunto_mapa[0]), 2));

    return log(sol/delta0);
}

int main()
{
    double h = 0.001;
    double E=1.;
    int k;

    double phi[n], psi[n];
    for(k=0; k<n; k++)
    {
        phi[k] = k*0.01;
        psi[k] = k*0.01;
    }

    double phiPunto[n], psiPunto[n], Pphi[n], Ppsi[n];
    for(k=0; k<n; k++)
    {
        phiPunto[k] = sqrt(E- 2*g*(1-cos(phi[k])) - g*(1-cos(psi[k])));;
        psiPunto[k] = 0.;

        Pphi[k] = 2*phiPunto[k];
        Ppsi[k] = phiPunto[k]*cos(phi[k]-psi[k]);
    }
    
    double exponente = 0.;
    double t, k1[N],k2[N], k3[N], k4[N], x[2], y[2], rx, ry, dist, Energia[n], phiPunto_mapa[n], psiPunto_mapa[n], delta0;

    delta0 = sqrt(pow(psi[1]-psi[0], 2)+ pow(psiPunto[1]-psiPunto[0], 2));

    int i, j;

    FILE *f1, *f2, *f3, *f4, *f5, *f6, *f7;

    f1 = fopen("Pendulo_E_1_Varios.txt", "w");
    //f2 = fopen("Distancia_E_3.txt", "w");
    f3 = fopen("E_1_Varios.txt", "w");
    f4 = fopen("Poincare_E_1_psi_psiPunto.txt", "w");
    f5 = fopen("Exponente_Lyapunov_E_1.txt", "w");
    f6 = fopen("Poincare_E_1_psi_phiPunto.txt", "w");
    f7 = fopen("Poincare_E_1_psi_phi.txt", "w");

    t=0.;
    for(i=0; i<1e4; i++)
    {
        for(k=0; k<n; k++)
        {
            k1[0] = h*func_phi(phi[k], psi[k], Pphi[k], Ppsi[k]);
            k1[1] = h*func_psi(phi[k], psi[k], Pphi[k], Ppsi[k]);
            k1[2] = h*func_Pphi(phi[k], psi[k], Pphi[k], Ppsi[k], h1(phi[k], psi[k], Pphi[k], Ppsi[k]), h2(phi[k], psi[k], Pphi[k], Ppsi[k]));
            k1[3] = h*func_Ppsi(phi[k], psi[k], Pphi[k], Ppsi[k], h1(phi[k], psi[k], Pphi[k], Ppsi[k]), h2(phi[k], psi[k], Pphi[k], Ppsi[k]));

            k2[0] = h*func_phi(phi[k] + k1[0]/2, psi[k] +k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2);
            k2[1] = h*func_psi(phi[k] + k1[0]/2, psi[k] +k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2);
            k2[2] = h*func_Pphi(phi[k] + k1[0]/2, psi[k] +k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2, h1(phi[k] + k1[0]/2, psi[k] + k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2), h2(phi[k] + k1[0]/2, psi[k] + k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2));
            k2[3] = h*func_Ppsi(phi[k] + k1[0]/2, psi[k] +k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2, h1(phi[k] + k1[0]/2, psi[k] + k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2), h2(phi[k] + k1[0]/2, psi[k] + k1[1]/2, Pphi[k] + k1[2]/2, Ppsi[k] + k1[3]/2));

            k3[0] = h*func_phi(phi[k] + k2[0]/2, psi[k] +k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2);
            k3[1] = h*func_psi(phi[k] + k2[0]/2, psi[k] +k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2);
            k3[2] = h*func_Pphi(phi[k] + k2[0]/2, psi[k] +k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2, h1(phi[k] + k2[0]/2, psi[k] + k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2), h2(phi[k] + k2[0]/2, psi[k] + k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2));
            k3[3] = h*func_Ppsi(phi[k] + k2[0]/2, psi[k] +k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2, h1(phi[k] + k2[0]/2, psi[k] + k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2), h2(phi[k] + k2[0]/2, psi[k] + k2[1]/2, Pphi[k] + k2[2]/2, Ppsi[k] + k2[3]/2));

            k4[0] = h*func_phi(phi[k] + k3[0], psi[k] +k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]);
            k4[1] = h*func_psi(phi[k] + k3[0], psi[k] +k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]);
            k4[2] = h*func_Pphi(phi[k] + k3[0], psi[k] +k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3], h1(phi[k] + k3[0], psi[k] + k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]), h2(phi[k] + k3[0], psi[k] + k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]));
            k4[3] = h*func_Ppsi(phi[k] + k3[0], psi[k] +k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3], h1(phi[k] + k3[0], psi[k] + k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]), h2(phi[k] + k3[0], psi[k] + k3[1], Pphi[k] + k3[2], Ppsi[k] + k3[3]));

            phi[k] = phi[k]+ (1./6.)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); 
            psi[k] = psi[k] + (1./6.)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
            Pphi[k] = Pphi[k] + (1./6.)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
            Ppsi[k] = Ppsi[k] + (1./6.)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
            
            //Obtengo las posiciones en cartesianas de cada bola de cada péndulo
            x[0] = sin(phi[k]);
            y[0] = -1.*cos(phi[k]);

            x[1] = sin(phi[k]) + sin(psi[k]);
            y[1] = -1.*cos(phi[k]) -1.*cos(psi[k]);

            //Muestro en un fichero la posición cada 10 iteraciones
            if(i%10 ==0)
            {
                for(j=0; j<=1; j++)
                {
                    fprintf(f1, "%lf,\t%lf\n", x[j], y[j]);
                }
            }

            //Calculo la energía del sistema, debe conservarse
            Energia[k] = (Pphi[k]*Pphi[k] + 2*Ppsi[k]*Ppsi[k] - 2*Pphi[k]*Ppsi[k]*cos(phi[k]-psi[k]))/(2*(1+pow((sin(phi[k]-psi[k])), 2))) + 2*g*(1-cos(phi[k])) + g*(1-cos(psi[k]));
            fprintf(f3, "%lf\n", Energia[k]);
            fprintf(f3, "\n");

            //Obtengo las coordenadas necesarias para obtener los distintos mapas de Pincaré y muestro cada 10 iteraciones
            phiPunto_mapa[k] = (Pphi[k] - Ppsi[k]*cos(phi[k]-psi[k]))/(1+pow(sin(phi[k]-psi[k]), 2));
            psiPunto_mapa[k] = (-Pphi[k]*cos(phi[k]-psi[k]) + 2*Ppsi[k])/(1+pow(sin(phi[k]-psi[k]), 2));
            if(i%10 ==0)
            {
                fprintf(f4, "%lf,\t%lf\n", psi[k], psiPunto_mapa[k]);
                fprintf(f6, "%lf,\t%lf\n", psi[k], phiPunto_mapa[k]);
                fprintf(f7, "%lf,\t%lf\n", psi[k], phi[k]);
            }

        }

        //Hago un salto de línea en el fichero de las posicines y mapas de Poincaré
        if(i%10 ==0)
        {
            fprintf(f1, "\n");
            fprintf(f4, "\n");
        }
        
        //Obtengo el exponente de Lyapunov cada 10 iteraciones
        if(i%10 ==0 && i != 0)
        {
            exponente += lyapunov(psi, psiPunto_mapa, delta0, h)/(i*h);
            fprintf(f5, "%i \t%lf\n", i, exponente);
        }
        
        //Compruebo que la distancia entre ambas bolas del péndulo es siempre la misma e igual a 1 (solo se ha comprobado para el caso de tener un único péndulo doble)
        rx = x[1] - x[0];
        ry = y[1] - y[0];

        dist =sqrt(rx*rx + ry*ry);
        //fprintf(f2, "eje x: %lf, eje y: %lf, distancia: %lf\n", rx, ry, dist);

        t=t+h;
    }    

    fclose(f1);
    //fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);

    return 0;
}