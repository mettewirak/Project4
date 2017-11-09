
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <mpi.h>
using namespace std;
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
ofstream ofile;
ofstream ofile2;


int** create_matrix_2(int n);
void delete_matrix_2(int n, int** state);
void print_matrix_to_screen(int n, int** state);
void initialize(int** state, double temperature, int& E, int& M, int n);
int periodic(int position, int N, int jump);
void Metropolis(int N, int **spins, int& E, int& M, double *w);
double r2();
void output(int N, int mc_cycles, double temperature, double* expectation_values);


int main(int argc, char* argv[])
{
    int numprocs, my_rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int N = 20;
    int** spins;
    spins = create_matrix_2(N);

    double final_temp = 5;
    double initial_temp = 1;
    int mc_cycles = pow(10,5);
    double temp_step = 0.50;
    int energy = 0;
    int magnetization = 0;
    double expectation_values[5]; // E, E*E, M, M*M, |M|
    double expectation_values_global[5];
    double w[17];

    srand(time(NULL)+6886456*my_rank);
    if(my_rank==0){
        ofile.open("Resultater.txt");}
    //ofile.open("Resultater_T_1.txt");
    //ofile.open("Resultater_T_2_4.txt");
    //ofile.open("Resultater_T_1_opp.txt");
    //ofile.open("Resultater_T_2_4_opp.txt");
    //ofile2.open("MC-utvikling.txt");

    for(double temperature = initial_temp; temperature <= final_temp; temperature += temp_step){

        for(int i=0; i<5; i++){
            expectation_values[i] = 0.0;}
        for(int de = -8; de <= 8; de ++){
            w[de+8] = 0.0;}
        for(int de = -8; de <= 8; de +=4){
            w[de+8] = exp(-de/temperature);}

        /*for(int i=0; i<17; i++){
            cout << "w[" << i << "]: " << w[i] << endl;
        }*/

        initialize(spins, temperature, energy, magnetization, N);
        //print_matrix_to_screen(N, spins);

        for(int cycles = 0; cycles < mc_cycles; cycles++){
            Metropolis(N, spins, energy, magnetization, w);
            expectation_values[0] += energy;
            expectation_values[1] += energy*energy;
            expectation_values[2] += magnetization;
            expectation_values[3] += magnetization*magnetization;
            expectation_values[4] += fabs(magnetization);
            //output(N, cycles, temperature, expectation_values);

            //ofile2 << setiosflags(ios::showpoint | ios::uppercase);
            //ofile2 << setw(15) << setprecision(8) << expectation_values[0]/((double)(mc_cycles));
            //ofile2 << setw(15) << setprecision(8) << expectation_values[3]/((double)(mc_cycles));
            //ofile2 << endl;
            //print_matrix_to_screen(N, spins);
        }
        MPI_Allreduce( expectation_values, expectation_values_global, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int i=0;i<5;i++){expectation_values_global[i]/=2;}
        //print_matrix_to_screen(N, spins);

        if(my_rank==0){output(N, mc_cycles, temperature, expectation_values_global);}
    }

    delete_matrix_2(N, spins);
    ofile.close();
    //ofile2.close();
    MPI_Finalize ();

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

    // Initialiserer til at alle spin peker opp hvis temperaturen er lav. Skal ellers fortsette som det var i temperaturen før.
    for(int y = 0; y < n; y++){
        for(int x = 0; x < n; x++){
            //if(temperature < 1.3)
            state[x][y] = 1; //Alle spinn peker opp
            //int i =rand()%2; //Random
            //if(i==1){state[x][y]=1;}
            //else{state[x][y]=-1;}

        }
    }

    E =0;
    M = 0;

    for(int y=0; y<n; y++){
        for(int x=0; x<n; x++){
            E += -J*state[x][y]*(state[x][periodic(y, n, -1)] + state[periodic(x, n, -1)][y]);
            M += state[x][y];
        }
    }
}


int periodic(int position, int N, int jump){ // Mulig å gjøre dette uten if? Det er tydeligvis dårlig

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


void Metropolis(int N, int **spins, int& E, int& M, double *w){

    for(int x = 0; x<N; x++){
        for(int y = 0; y<N; y++){

            double random_number = RandomNumberGenerator(gen);
            int deltaE = 0;
            int random_spin_x = int(RandomNumberGenerator(gen)*N);
            int random_spin_y = int(RandomNumberGenerator(gen)*N);

            deltaE = 2*spins[random_spin_x][random_spin_y]*
                    (spins[random_spin_x][periodic(random_spin_y, N, 1)] + spins[random_spin_x][periodic(random_spin_y, N, -1)] +
                    spins[periodic(random_spin_x, N, 1)][random_spin_y] + spins[periodic(random_spin_x, N, -1)][random_spin_y]);
            //double a=2*spins[random_spin_x][random_spin_y]*(spins[random_spin_x][periodic(random_spin_y, N, 1)]);
            //double b=spins[random_spin_x][periodic(random_spin_y, N, -1)];
            //cout<< a<< "b=" <<b<<endl;
            if(random_number <= w[deltaE+8]){
                spins[random_spin_x][random_spin_y] *= -1;
                E += deltaE;
                M += 2*spins[random_spin_x][random_spin_y];
                //cout<< E<<" de "<<deltaE<<endl;
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


