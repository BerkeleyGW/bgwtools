//#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iomanip>

using namespace std;

int main (int argc, char** argv) {
	double kx,ky,kz;
	int Nb,b;
	FILE *fp;
	int buf1,buf2;
	int kpts_cnt;

	if (argc!=2){
		cout << "Expecting at least one argument." <<endl;
		exit(1);
	}

	//setiosflags(ios::fixed);
	cout.setf(ios::fixed);
	cout.precision(5);

	typedef double vect[3];
	vect* Kpts=NULL;
	vect v;

	double *E_LDA=NULL, *E_GW=NULL;
	int offset;

	fp = fopen(argv[1],"r");
	kpts_cnt = 0;
	while(fscanf(fp,"%lf %lf %lf %d",v,v+1,v+2,&Nb)==4){
		Kpts  = (vect*)   realloc(Kpts , sizeof(vect)     *(kpts_cnt+1));
		E_LDA = (double*) realloc(E_LDA, sizeof(double)*Nb*(kpts_cnt+1));
		E_GW  = (double*) realloc(E_GW , sizeof(double)*Nb*(kpts_cnt+1));
			
		memcpy(Kpts+kpts_cnt, v, sizeof(vect));
		for (b=0; b<Nb; b++){
			offset = kpts_cnt*Nb;
			if (fscanf(fp,"%d %d %lf %lf", &buf1, &buf2,\
				E_LDA+offset+b, E_GW+offset+b) != 4){
				cout<<"Could not read band information!"<<endl;
				exit(1);
			}
		}
		kpts_cnt++;
	}

	for (int k = 0; k < kpts_cnt; ++k){
		for (int i=0; i<3; i++){
			cout << Kpts[k][i] << " ";
		}
		for (int i=0; i<Nb; i++){
			offset = Nb*k;
			//cout << E_LDA[offset+i] << " ";
			cout << E_GW [offset+i] << " ";
		}
		cout << endl;
	}
	free(Kpts);
	free(E_LDA);
	free(E_GW);
}
