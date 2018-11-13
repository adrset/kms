#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>

#include <TCanvas.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TGraph.h>
#include "Vec3.hpp"
using namespace std;

const double a = 0.38;
const double T_0 = 100;
const double k = 0.00831;
const double m = 39.948;
const double R = 0.38;
const double L = 2.3;
const double f = 10000.0;
const double e = 1.0;
const double PI = 3.14159265359;

const unsigned int STEPS_O = 1000;
const unsigned int STEPS_XYZ = 25;
const unsigned int STEPS_OUT = 200;
const unsigned int n = 3;
const unsigned int N = n * n * n;

double pown(double a, int power){
	double ret = a;
	
	for(int ii = 0; ii < power-1; ii++ ){
		ret *= a;
	}

	return ret;
}

Double_t getMax(Double_t tab[], int size){
	Double_t max = tab[0];
	for(int i = 1;i< size ; i++){
		if(tab[i]>max)	
			max = tab[i];
		
	}
	
	return max;
	
}

Double_t getMin(Double_t tab[], int size){
	Double_t min = tab[0];
	for(int i = 1;i< size ; i++){
		if(tab[i]<min)	
			min = tab[i];
		
	}
	
	return min;
	
}

void calculateForces(std::vector<Vec3>& F_i, std::vector<Vec3>& rVectors, double& V, double& P){
	V = 0;
	P = 0;
	
	for(unsigned int i=0;i<N;i++){
		F_i[i].zero();
	}
	// GLOWNA PETLA
	for (unsigned int ii=0;ii<N;ii++){
		
		// V - (10) i dodac do V_calk
		double r_i = rVectors[ii].getLength(); // dlugosc wektora r_i
		//std::cout<< r_i << std::endl;
		double temp = r_i < L ? 0.0 : f*0.5 * (r_i - L) * (r_i - L);
		V += temp; // 													akumulacja V
		
		// F - (14) i dodac do F_i
		Vec3 tempF = r_i < L ? Vec3() : f*(L-r_i) * rVectors[ii] / r_i;
		F_i[ii] = F_i[ii] + tempF;
		
		// Cisnienie chwilowe (15)
		
		P += tempF.getLength();
		//std::cout<<r_i<<"<"<<L<<std::endl;
			// WEWNETRZNA PETLA
			// <!!!!!!!!>
			for(unsigned int jj=0;jj<ii;jj++){
				Vec3 len;
				len.setR(rVectors[ii]);
				len = len - rVectors[jj];
				double length = len.getLength();
				// pot. par atomowych (9)
				double V_p = e * (pown((R / length),12) - 2.0 *  pown((R / length),6));
				V += V_p;// 											akumulacja V
				
				// Sily miedzyatomowe Fi Fj (13)
				Vec3 F_ij;
				F_ij = 12.0 * e *((pown(R/length,12) - pown(R/length,6))) * ((rVectors[ii] - rVectors[jj]) / (length * length));
				F_i[ii] =  F_i[ii] + F_ij;
				F_i[jj] =  F_i[jj] + F_ij * -1.0;
				
			}
			// <!!!!!!!!>
		
	}
	P *= 1.0 / (4.0 * PI * L * L);
	
}

int main(int argc, char** argv){
	if(argc<4)
		return 1;

	bool paint = false;
	if(argc>4)
		paint = true;
	std::vector<std::string> argList;
	for(int i=0;i<argc;i++){
		
		argList.push_back(argv[i]);
	    //now you can access argList[n]
	}
	 
	double dt = std::stod(argList[2]);
	double span = std::stod(argList[3]);
	int STEPS = span;
	//std::cout<<"dt: " << dt << std::endl;
	std::ofstream xyz("output/xyz" + argList[1] + ".dat");
	std::ofstream out("output/out" + argList[1] + ".dat");
	srand48(time(NULL));
	Vec3 b0 (a, 0.0, 0.0);
	Vec3 b1 (a / 2.0, a / 2.0 * sqrt(3), 0.0);
	Vec3 b2 (a / 2.0, a / 6.0 * sqrt(3), a * sqrt(2.0/ 3.0));
	
	ofstream plik("dane");
	unsigned int writes = 0;
	std::vector<Vec3> rVectors;
	std::vector<Vec3> pVectors;
	Vec3 pTotal(0,0,0);
	Vec3 tmp1, tmp2, tmp3;
	for(unsigned int ii=0;ii<n;ii++){
		for(unsigned int jj=0;jj<n;jj++){
			for(unsigned int kk=0;kk<n;kk++){
				
				double xx=drand48(), yy=drand48(), zz=drand48();
				
				xx = xx > 0.5 ? 1.0 : -1.0;
				yy = yy > 0.5 ? 1.0 : -1.0;
				zz = zz > 0.5 ? 1.0 : -1.0;
				
				//double kineticEnergy = -0.5 * k * T_0 * log(drand48());	
				Vec3 p(xx * sqrt(2*m*(-0.5 * k * T_0 * log(drand48()))), yy * sqrt(2*m*(-0.5 * k * T_0 * log(drand48()))), zz * sqrt(2*m*(-0.5 * k * T_0 * log(drand48()))));
				pVectors.push_back(p);	
				pTotal = p + pTotal;	
				
				tmp1.setR(b0);
				tmp1= tmp1 * (ii - (( n - 1 ) / 2.0 ));
				tmp2.setR(b1);
				tmp2 = tmp2 * (jj - (( n - 1 ) / 2.0 ));
				tmp3.setR(b2);
				tmp3 = tmp3 * (kk - (( n - 1 ) / 2.0 ));
				
				Vec3 r(tmp1);
				
				r = r + tmp2;
				r = r + tmp3;
				
				rVectors.push_back(r);

				
			}
		}
	}
	pTotal = pTotal / ((float)N);
	
	for(unsigned int ii=0;ii<N;ii++){
		pVectors[ii] = pVectors[ii] - pTotal;
	}
	
	//std::cout<< pTotal << std::endl;		
	for(unsigned int ii=0;ii<N;ii++){
		//std::cout<< rVectors[ii] << std::endl;		
	}
	
	// Lab 2
	// Null Forces P and V
	
	double V = 0.0;
	double P = 0.0;
	double E_kin = 0;
	double T = 0;
	double potentialTmp = 0;
	double energyTmp = 0;
	double pTmp = 0;
	double tTmp = 0;
	std::vector<Vec3> F_i;
	for(unsigned int i=0;i<N;i++){
		F_i.push_back(Vec3(0,0,0));
	}
	
	calculateForces(F_i, rVectors, V, P);
	std::cout<<a<<"\t"<<V<<std::endl;
	for(unsigned int i=0;i<N;i++){
				
	}
	
	for( int s = 1;s< STEPS +STEPS_O; s++){
		
		for(unsigned int i=0;i<N;i++){
			pVectors[i] = pVectors[i] + 0.5 * F_i[i] * dt;
			rVectors[i] = rVectors[i] + 1.0 / m * pVectors[i] * dt;	
		}
		
		calculateForces(F_i, rVectors, V, P);
		
		for(unsigned int i=0;i<N;i++){
			pVectors[i] = pVectors[i] + 0.5 * F_i[i] * dt;	
		}
		
		E_kin = 0;
		T = 0;
		for(unsigned int i=0;i<N;i++){
			E_kin += pVectors[i].getLengthSquared() / (2.0 * m);
		}
		T = 2.0 / (3.0 * N * k) * E_kin;
		if(s > STEPS_O){
			potentialTmp += V;
			energyTmp += E_kin;
			tTmp += T;
			pTmp += P;	
		}
		if(s % STEPS_OUT == 0 && s > STEPS_O){
			writes++;
			out<<s*dt<<"\t"<<V<<"\t"<<energyTmp/STEPS_OUT<<"\t"<<(potentialTmp+energyTmp)/STEPS_OUT<<"\t"<<tTmp/STEPS_OUT<<"\t"<<pTmp/STEPS_OUT<<"\t"<<std::endl;
			potentialTmp = 0;
			energyTmp = 0;
			pTmp = 0;
			tTmp = 0;
		}
		if(s % STEPS_XYZ == 0){
			for(unsigned int i=0;i<N;i++){
				xyz<<rVectors[i]<<std::endl;
			}
			xyz<<std::endl<<std::endl;
		}

		if(s % STEPS_O){
			std::cout<< s << std::endl;
			//printf("\r%f%\t\t",(double) s / ((double) STEPS +STEPS_O) * 100.0);
		}
		
		
	}
	std::cout<<std::endl;
	plik.close();
	
	out.close();

	
	if(paint){
		std::ifstream in("output/out" + argList[1] + ".dat");
//std::cout<<"output/out" + argList[1] + ".dat";

		Double_t time[writes], potential[writes], Ek[writes], Ec[writes], T[writes], P[writes];
			Double_t av = 0;
			
			
			for(Int_t ii=0;!in.eof();ii++){	    
				
				

				in>>time[ii]>>potential[ii]>>Ek[ii]>>Ec[ii]>>T[ii]>>P[ii];
				av+=potential[ii];
				

}

		Double_t stdev = 0;
		av/=writes;
		std::cout<< R << "\t" << av <<std::endl;
		
		TApplication* rootapp = new TApplication("Pusta",0, argv);
		
		TCanvas *tc=new TCanvas("c1","potential",0,0,400,400);
		TGraph* graph1=new TGraph(writes-1,time, potential);
		graph1->SetMarkerStyle(19);
		tc->DrawFrame(time[0],getMin(potential, writes),time[writes-1],getMax(potential, writes));
		graph1->SetTitle("Potential energy; Time [ps]; Energy");
		graph1->Draw();

		TCanvas *tc2=new TCanvas("c2","Ek",400,0,400,400);

		TGraph* graph2=new TGraph(writes-1,time, Ek);
		graph2->SetMarkerStyle(19);
		tc2->DrawFrame(time[0],getMin(Ek, writes),time[writes-1],getMax(Ek, writes));
		graph2->SetTitle("Kinetic energy; Time [ps]; Energy");
		graph2->Draw();
		TCanvas *tc3=new TCanvas("c3","Ec",0,400,400,400);
		TGraph* graph3=new TGraph(writes-1,time, Ec);
		graph3->SetTitle("Total energy; Time [ps]; Energy");
		tc3->DrawFrame(time[0],-100+getMin(Ec, writes),time[writes-1],100+getMax(Ec, writes));
		graph3->SetMarkerStyle(19);
		graph3->Draw();
		TCanvas *tc4=new TCanvas("c4","T",400,400,400,400);
		TGraph* graph4=new TGraph(writes-1,time, T);
		tc4->DrawFrame(time[0],getMin(T, writes),time[writes-1],getMax(T, writes));
		graph4->SetMarkerStyle(19);
		graph4->SetTitle("Temperature; Time [ps]; T [K]");
		graph4->Draw();

		TCanvas *tc5=new TCanvas("c5","P",400,400,400,400);
		TGraph* graph5=new TGraph(writes-1,time, P);
		tc5->DrawFrame(time[0],getMin(P, writes),time[writes-1],getMax(P, writes));
		graph5->SetMarkerStyle(19);
		graph5->SetTitle("Temperature; Time [ps]; T [K]");
		graph5->Draw();
	    	rootapp->Run();
		
		
		out.close();
		
		
		
	}
	
}
