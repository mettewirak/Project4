#include <iostream>
#include "lib.h"
using namespace std;

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
    int p= ran0(0)
cout<< p
	return 0;
}

//double MonteCarlo(double *x){



//}
