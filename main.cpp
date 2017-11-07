#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
ofstream ofile;
ofstream ofile2;


int** create_matrix_2(int n);
void delete_matrix_2(int n, int** state);
void print_matrix_to_screen(int n, int** state);
void initialize(int** state, double temperature, int& E, int& M, int n);
int periodic(int position, int N, int jump);
void Metropolis(int N, double random_number, int **spins, int& E, int& M, double *w);
double r2();
void output(int N, int mc_cycles, double temperature, double* expectation_values);


int main()
{
    int N = 20;
    int** spins;
    spins = create_matrix_2(N);

    double final_temp = 4;
    double initial_temp = 1;
    int mc_cycles = pow(10,5);
    double temp_step = 0.50;
    int energy = 0;
    int magnetization = 0;
    double expectation_values[5]; // E, E*E, M, M*M, |M|
    double w[17];

    srand(time(NULL));

    ofile.open("Resultater.txt");
    //ofile2.open("MC-utvikling.txt");

    for(double temperature = initial_temp; temperature <= final_temp; temperature += temp_step){

        for(int i=0; i<5; i++){
            expectation_values[i] = 0.0;}
        for(int de = -8; de <= 8; de ++){
            w[de+8] = 0.0;}
        for(int de = -8; de <= 8; de +=4){
            w[de+8] = exp(-de/temperature);}

        initialize(spins, temperature, energy, magnetization, N);

	// MONTE CARLO
	// Må lage en random number generator.
	// Gir ut #dim forskjellige tall (i første omgang 2). Disse skaleres så til en verdi mellom 1 og N (antall spin i hver dimensjon).
	// Dette beskriver da en spesifikk partikkel som evt. skal endre spin.
	// Monte Carlo kjøres i en for-løkke ca 10^6 ganger.  

        for(int cycles = 0; cycles < mc_cycles; cycles++){
            Metropolis(N, r2(), spins, energy, magnetization, w);
            expectation_values[0] += energy;
            expectation_values[1] += energy*energy;
            expectation_values[2] += magnetization;
            expectation_values[3] += magnetization*magnetization;
            expectation_values[4] += fabs(magnetization);

            //ofile2 << setiosflags(ios::showpoint | ios::uppercase);
            //ofile2 << setw(15) << setprecision(8) << expectation_values[0]/((double)(mc_cycles));
            //ofile2 << setw(15) << setprecision(8) << expectation_values[3]/((double)(mc_cycles));
            //ofile2 << endl;

        }

        cout << "Temperature " << temperature << endl;

        print_matrix_to_screen(N, spins);
        output(N, mc_cycles, temperature, expectation_values);
    }

    delete_matrix_2(N, spins);
    ofile.close();
    //ofile2.close();
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
    int sum = 0;

    // Initialiserer til at alle spin peker opp hvis temperaturen er lav. Skal ellers fortsette som det var i temperaturen før.
    for(int y = 0; y < n; y++){
        for(int x = 0; x < n; x++){
            if(temperature < 1.3)
                state[x][y] = 1;
        }
    }

    E = M = 0;

    for(int y=0; y<n; y++){
        for(int x=0; x<n; x++){

            sum = state[x][periodic(y, n, -1)] + state[x][periodic(y, n, +1)]
                    + state[periodic(x, n, -1)][y] + state[periodic(x, n, +1)][y];

            E += -J*state[x][y]*sum;
            sum = 0;

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

    for(int x = 0; x<N; x++){
        for(int y = 0; y<N; y++){

            int deltaE = 0;
            int random_spin_x = rand()%N;
            int random_spin_y = rand()%N;

            deltaE = 2*spins[random_spin_x][random_spin_y]*
                (spins[random_spin_x][periodic(random_spin_y, N, 1)] + spins[random_spin_x][periodic(random_spin_y, N, -1)] +
                spins[periodic(random_spin_x, N, 1)][random_spin_y] + spins[periodic(random_spin_x, N, -1)][random_spin_y]);

            if(random_number >= w[deltaE+8]){
                spins[random_spin_x][random_spin_y] *= -1;
                E += deltaE;
                M += 2*spins[random_spin_x][random_spin_y];
            }
        }
    }
}


double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}


void output(int N, int mc_cycles, double temperature, double* expectation_values)
{
double norm = 1.0/((double) (mc_cycles)); // divided by number of cycles
double E_ExpectationValues = expectation_values[0]*norm;
double E2_ExpectationValues = expectation_values[1]*norm;
double M_ExpectationValues = expectation_values[2]*norm;
double M2_ExpectationValues = expectation_values[3]*norm;
double Mabs_ExpectationValues = expectation_values[4]*norm;

// all expectation values are per spin, divide by 1/NSpins/NSpins

double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/(N*N);
double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/(N*N);

ofile << setiosflags(ios::showpoint | ios::uppercase);
ofile << setw(15) << setprecision(8) << temperature;
ofile << setw(15) << setprecision(8) << E_ExpectationValues/(N*N);
ofile << setw(15) << setprecision(8) << Evariance/(temperature*temperature);
ofile << setw(15) << setprecision(8) << M_ExpectationValues/(N*N);
ofile << setw(15) << setprecision(8) << Mvariance/temperature;
ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/(N*N);
ofile << endl;
}
