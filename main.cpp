#include <iostream>
using namespace std;

double** create_matrix(int dim, int n);
void delete_matrix(int dim, int n, double** state);
double ising_first(double** state, int dim, int n);
void ising_diff(double** state, double* energy);

int main()
{
	// MÅ TA INN
	// - Dimensjoner
	// - Antall partikler i hver dimensjon
	// - Hvor mange ganger Monte Carlo skal kjøres
	// - En grense for hvor høy sannsynligheten skal være før noe gjennomføres?
	// - Temperaturintervallet funksjonen skal kjøres over.


	// Skal gjøres for mange temperaturer. Må ha en for-løkke for dette ytterst.

	// Lagre verdier i en peke-matrise innenfor temperaturløkka? Og så skrive til til for hver temperatur?
	// Ha en egen klasse for systemet? Virker som en god ide. Kan lagre energier osv.

	// MONTE CARLO
	// Må lage en random number generator.
	// Gir ut #dim forskjellige tall (i første omgang 2). Disse skaleres så til en verdi mellom 1 og N (antall spin i hver dimensjon).
	// Dette beskriver da en spesifikk partikkel som evt. skal endre spin.
	// Monte Carlo kjøres i en for-løkke ca 10^6 ganger.  


	// ISING
	// Regner ut energien til hele systemet slik det er nå (dette kan sannsynligvis lagres på), og energien hvis partikkelen funnet
	// i Monte Carlo skifter spinn. Dette tas med inn i Metropolis. Skal bruke periodiske grensebetingelser. 
	// Bør være en egen funksjon.


	// METROPOLIS
	// Regner ut sannsynligheten for at endringen blir gjennomført, gitt av W = AT. Må ha en grense for hvor høy sannsynligheten skal
	// være for at endringen skal gjøres?



	// RESULTATER
	// Gi ut gjennomsnittsenergi, gjennomsnittsmagnetisering, varmekapasitet og suseptebilitet.


	return 0;
}

double** create_matrix(int dim, int n){

    double **state;
        state = new double*[dim];

        for(int i=0; i<total_planets; i++){
            state[i] = new double[n];
        }

    return state;
}

void delete_matrix(int dim, int n, double** state){

    for(int i=0; i<n; i++){
        delete [] state[i];
    }
    delete [] state;
}

double ising_first(double **state, int dim, int n){ // MÅ TILPASSES ANDRE DIMENSJONER ENN 2

    double J; // Coupling constant
    double energy = 0.0;
    double temp = 0.0;

    for(int x=0; x<n; x++){

        for(int y=0; y<n; y++){

            // KANT VENSTRE
            if(x == 0)
                temp += state[n-1][y];
            else
                temp += state[x-1][y];

            // KANT OPPE
            if(y == 0)
                temp += state[x][n-1];
            else
                temp += state[x][y-1];

            // KANT HØYRE
            if(x == (n-1))
                temp += state[0][y];
            else
                temp += state[x+1][y];

            // KANT NEDE
            if( y == (n-1))
                temp += state[x][0];
            else
                temp += state[x][y+1];

            energy += state[x][y]*temp;
            temp = 0.0;
        }
    }
    return (-J*energy);
}

void ising_diff(double **state, double *energy){

    double J; // Coupling constant
    int temp = 0.0;

    // KANT VENSTRE
    if(x == 0)
        temp += state[n-1][y];
    else
        temp += state[x-1][y];

    // KANT OPPE
    if(y == 0)
        temp += state[x][n-1];
    else
        temp += state[x][y-1];

    // KANT HØYRE
    if(x == (n-1))
        temp += state[0][y];
    else
        temp += state[x+1][y];

    // KANT NEDE
    if( y == (n-1))
        temp += state[x][0];
    else
        temp += state[x][y+1];

    *energy = (2*J*state[x][y]*temp);
}
