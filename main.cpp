#include <iostream>
#include <cmath>
using namespace std;

int** create_matrix_2(int n);
void delete_matrix_2(int n, int** state);
void print_matrix_to_screen(int n, int** state);
void initialize(int** state, double temperature, int& E, int& M, int n);
int periodic(int position, int N, int jump);
void Metropolis(int N, double random_number, int **spins, int& E, int& M, double *w);
double r2();


int main()
{
	// MÅ TA INN
	// - Dimensjoner
	// - Antall partikler i hver dimensjon
	// - Hvor mange ganger Monte Carlo skal kjøres
	// - En grense for hvor høy sannsynligheten skal være før noe gjennomføres?
	// - Temperaturintervallet funksjonen skal kjøres over.

    //int dimension = 2;
    int N = 2;
    int** spins;
    double final_temp = 3, initial_temp = 1;
    double temp_step = (final_temp-initial_temp)/pow(10,1);
    int energy = 0, magnetization = 0;


    // Skal gjøres for mange temperaturer. Må ha en for-løkke for dette ytterst.

    for(int temperature = initial_temp; temperature < final_temp; temperature += temp_step){

        // Precalculate probabilities according to the boltzmann distribution for the different energy changes
        double w[17];
        for(int de = -8; de < 8; de ++)
            w[de] = 0.0;
        for(int de = -8; de < 8; de +=4)
            w[de] = exp(-de/temperature);

    // Lagre verdier i en peke-matrise innenfor temperaturløkka? Og så skrive til til for hver temperatur?

        spins = create_matrix_2(N);

        initialize(spins, temperature, energy, magnetization, N);

        cout << "The energy initial energy is " << energy << ", and the magnetization is " << magnetization << ".\n";


	// MONTE CARLO
	// Må lage en random number generator.
	// Gir ut #dim forskjellige tall (i første omgang 2). Disse skaleres så til en verdi mellom 1 og N (antall spin i hver dimensjon).
	// Dette beskriver da en spesifikk partikkel som evt. skal endre spin.
	// Monte Carlo kjøres i en for-løkke ca 10^6 ganger.  

        for(int i=0; i < pow(10, 1); i++){

            Metropolis(N, r2(), spins, energy, magnetization, w);
            print_matrix_to_screen(N, spins);

        }
	// RESULTATER
	// Gi ut gjennomsnittsenergi, gjennomsnittsmagnetisering, varmekapasitet og suseptebilitet.
    }

    delete_matrix_2(N, spins);

	return 0;
}


int** create_matrix_2(int n){

    int **state;
    state = new int*[n];

    for(int i=0; i < n; i++){
        state[i] = new int[n];
    }

    for(int y=0; y<n; y++){
        for(int x=0; x<n; x++){
            state[x][y] = 0;
        }
    }

    return state;
}


void delete_matrix_2(int n, int** state){

    for(int i=0; i<n; i++){
        delete [] state[i];
    }
    delete [] state;
}


void print_matrix_to_screen(int n, int **state){

    cout << "\n";
    for(int y = 0; y < n; y++){
        for(int x = 0; x < n; x++){

            cout << state[x][y] << "  ";
        }
        cout << "\n";
    }
}


void initialize(int **state, double temperature, int& E, int& M, int n){ // Uses periodic boundary conditions.

    int J = 1; // Coupling constant
    int temp_energy = 0;

    // Initialiserer til at alle spin peker opp hvis temperaturen er lav. Skal ellers fortsette som det var i temperaturen før.
    for(int y = 0; y < n; y++){
        for(int x = 0; x < n; x++){
            if(temperature < 1.5)
                state[x][y] = 1;
        }
    }

    E = M = 0;

    for(int y=0; y<n; y++){
        for(int x=0; x<n; x++){

            // KANT VENSTRE
            if(x == 0)
                temp_energy += state[n-1][y];
            else
                temp_energy += state[x-1][y];

            // KANT OPPE
            if(y == 0)
                temp_energy += state[x][n-1];
            else
                temp_energy += state[x][y-1];

            // KANT HØYRE
            if(x == (n-1))
                temp_energy += state[0][y];
            else
                temp_energy += state[x+1][y];

            // KANT NEDE
            if( y == (n-1))
                temp_energy += state[x][0];
            else
                temp_energy += state[x][y+1];

            E += -J*state[x][y]*temp_energy;
            temp_energy = 0;

            M += state[x][y];
        }
    }
}


periodic(int position, int N, int jump){ // Mulig å gjøre dette uten if? Det er tydeligvis dårlig

    if((position+jump) < 0){
        return (N-1);
    }
    else if((position+jump) > (N-1)){
        return 0;
    }
    else{
        return (position + jump);
    }
}


void Metropolis(int N, double random_number, int **spins, int& E, int& M, double *w){


    // METROPOLIS
    // Regner ut sannsynligheten for at endringen blir gjennomført, gitt av W = AT. Må ha en grense for hvor høy sannsynligheten skal
    // være for at endringen skal gjøres?

    int deltaE = 0;
    int random_spin_x = rand()%2;
    int random_spin_y = rand()%2;

    deltaE = 2*spins[random_spin_x][random_spin_y]*
         (spins[random_spin_x][periodic(random_spin_y, N, 1)] + spins[random_spin_x][periodic(random_spin_y, N, -1)] +
         spins[periodic(random_spin_x, N, 1)][random_spin_y] + spins[periodic(random_spin_x, N, -1)][random_spin_y]);

    cout << "Spin to be flipped: " << random_spin_x << ", " << random_spin_y << endl;
    cout << "DeltaE: " << deltaE << endl;
    cout << "Compare " << random_number << " and " << w[deltaE+8] << endl;

    if(random_number > w[deltaE+8]){
         spins[random_spin_x][random_spin_y] *= -1;
         E += deltaE;
         M += 2*spins[random_spin_x][random_spin_y];
    }
}


double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}
