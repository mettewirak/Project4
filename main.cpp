#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
using namespace std;

std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

ofstream ofile;
ofstream ofile2;
ofstream ofile3;

int** create_matrix(int N);
void delete_matrix(int N, int** state);
void print_matrix_to_screen(int N, int** state);
void initialize(int** state, double temperature, int& E, int& M, int N);
int periodic(int position, int N, int jump);
void Metropolis(int N, int **spins, int& E, int& M, double *w, int &accepted_configs);
void update_probabilities(int energy, int* counter, int& count_total, int possible_energies);
void output_probabilities(int *counter, int& count_total, int possible_energies, int N);
void output(int N, int mc_cycles, double temperature, double* expectation_values);


int main()
{
    int N = 20;
    int** spins;
    spins = create_matrix(N);

    double initial_temp = 1;
    double final_temp = 1;
    double temp_step = 0.50;

    int mc_cycles = pow(10,5);
    int accepted_configs = 0;

    int energy = 0;
    int magnetization = 0;
    double expectation_values[5]; // E, E*E, M, M*M, |M|
    double w[17];

    // Variables needed for finding the probabilities of the energies at a given temperature. Can be commented out if that is not to be done.
    int possible_energies = N*N*2*2 + 1;
    int* counter;
    int count_total=0;
    counter = new int[possible_energies];
    for(int i = 0; i<possible_energies; i++)
        counter[i] = 0;

    ofile.open("Properties T=1 dt=0.5 MC=5.txt");
    ofile2.open("Probabilities T=2.4 dt=.5 MC=5.txt");
    //ofile3.open("Accepted configs T=1 dt=0.5 MC=5.txt");

    for(double temperature = initial_temp; temperature <= final_temp; temperature += temp_step){

        // Precalculate values for probabilities(deltaE)
        for(int i=0; i<5; i++){
            expectation_values[i] = 0.0;}
        for(int de = -8; de <= 8; de ++){
            w[de+8] = 0.0;}
        for(int de = -8; de <= 8; de +=4){
            w[de+8] = exp(-de/temperature);}

        initialize(spins, temperature, energy, magnetization, N);

        // Starting Monte Carlo
        for(int cycles = 0; cycles < mc_cycles; cycles++){

            Metropolis(N, spins, energy, magnetization, w, accepted_configs);

            if(cycles>(0.1*mc_cycles)){
                update_probabilities(energy, counter, count_total, possible_energies);
            }

            expectation_values[0] += energy;
            expectation_values[1] += energy*energy;
            expectation_values[2] += magnetization;
            expectation_values[3] += magnetization*magnetization;
            expectation_values[4] += fabs(magnetization);

            /*ofile3 << setiosflags(ios::showpoint | ios::uppercase);
            ofile3 << setw(15) << setprecision(8) << cycles;
            ofile3 << setw(15) << setprecision(8) << accepted_configs;
            ofile3 << endl;
            */
        }

        output_probabilities(counter, count_total, possible_energies, N);
        output(N, mc_cycles, temperature, expectation_values);
    }

    delete_matrix(N, spins);
    ofile.close();
    ofile2.close();
    //ofile3.close();

    return 0;
}


int** create_matrix(int N){

    int **state;
    state = new int*[N];

    for(int i=0; i < N; i++){
        state[i] = new int[N];
    }

    for(int y=0; y<N; y++){
        for(int x=0; x<N; x++){
            state[x][y] = 0;
        }
    }

    return state;
}


void delete_matrix(int N, int** state){

    for(int i=0; i<N; i++){
        delete [] state[i];
    }

    delete [] state;
}


void print_matrix_to_screen(int N, int **state){

    cout << "\n";
    for(int y = 0; y < N; y++){
        for(int x = 0; x < N; x++){

            cout << setiosflags(ios::showpoint | ios::uppercase);
            cout << setw(3) << setprecision(1) << state[x][y];
        }
        cout << endl;
    }
}


void initialize(int **state, double temperature, int& E, int& M, int N){

    int J = 1;                              // Coupling constant

    for(int y = 0; y < N; y++){             // If the simulation is only run on 1 temperature, the if-loop is not necessary.
        for(int x = 0; x < N; x++){
            if(state[x][y]==0){             // Only applied if this is the first time initialized is called
                state[x][y] = 1;
            }
            else if(temperature < 1.5){
                state[x][y] = 1;
            }
        }
    }

    E = 0;
    M = 0;

    for(int y=0; y < N; y++){
        for(int x=0; x < N; x++){
            E += -J*state[x][y]*(state[x][periodic(y, N, -1)] + state[periodic(x, N, -1)][y]);
            M += state[x][y];
        }
    }
}


periodic(int position, int N, int jump){

    if((position+jump) < 0){
        return (N-1);}
    else if((position+jump) > (N-1)){
        return 0;}
    else{
        return (position + jump);}
}


void Metropolis(int N, int **spins, int& E, int& M, double *w, int& accepted_configs){

    for(int x = 0; x<N; x++){
        for(int y = 0; y<N; y++){

            int deltaE = 0;
            double random_number = RandomNumberGenerator(gen);
            int random_spin_x = int(RandomNumberGenerator(gen)*N);
            int random_spin_y = int(RandomNumberGenerator(gen)*N);

            deltaE = 2*spins[random_spin_x][random_spin_y]*
                    (spins[random_spin_x][periodic(random_spin_y, N, -1)] + spins[random_spin_x][periodic(random_spin_y, N, 1)] +
                    spins[periodic(random_spin_x, N, 1)][random_spin_y] + spins[periodic(random_spin_x, N, -1)][random_spin_y]);

            if(random_number <= w[deltaE+8]){
                spins[random_spin_x][random_spin_y] *= -1;
                E += deltaE;
                M += 2*spins[random_spin_x][random_spin_y];
                accepted_configs +=1;
            }
        }
    }
}


void update_probabilities(int energy, int *counter, int &count_total, int possible_energies){

    counter[energy + ((possible_energies-1)/2)] += 1;
    count_total += 1;
}


void output_probabilities(int *counter, int &count_total, int possible_energies, int N){

    for(int i=0; i<possible_energies; i++){
        ofile2 << setiosflags(ios::showpoint | ios::uppercase);
        ofile2 << setw(15) << setprecision(8) << ((double)(i-((possible_energies-1)/2))/(double)(N*N));
        ofile2 << setw(15) << setprecision(8) << counter[i];
        ofile2 << setw(15) << setprecision(8) << (double)counter[i] / (double)(count_total);
        ofile2 << endl;
    }
}


void output(int N, int mc_cycles, double temperature, double* expectation_values)
{
    double norm = 1.0/((double) (mc_cycles));
    double E_ExpectationValues = expectation_values[0]*norm;
    double E2_ExpectationValues = expectation_values[1]*norm;
    double M_ExpectationValues = expectation_values[2]*norm;
    double M2_ExpectationValues = expectation_values[3]*norm;
    double Mabs_ExpectationValues = expectation_values[4]*norm;

    double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/(N*N);
    double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/(N*N);

    // All output expectation values are per spin.
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << E_ExpectationValues/(N*N);
    ofile << setw(15) << setprecision(8) << Evariance/(temperature*temperature);
    ofile << setw(15) << setprecision(8) << M_ExpectationValues/(N*N);
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/(N*N);
    ofile << endl;
}
