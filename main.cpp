#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
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
const double dt = 0.01;
const unsigned int STEPS = 10000;
const unsigned int STEPS_O = 1000;
const unsigned int STEPS_XYZ = 10;
const unsigned int STEPS_OUT = 100;
const unsigned int n = 5;
const unsigned int N = n * n * n;

void calculateForces(std::vector<Vec3>& F_i, std::vector<Vec3>& rVectors, double& V, double& P){
	V = 0;
	P = 0;
	for(unsigned int i=0;i<N;i++){
		F_i[i].zero();
	}
	// GLOWNA PETLA
	for (unsigned int ii=1;ii<N;ii++){
		
		// V - (10) i dodac do V_calk
		double r_i = rVectors[ii].getLength(); // dlugosc wektora r_i
		//std::cout<< r_i << std::endl;
		double temp = r_i < L ? 0.0 : 0.5 * (r_i - L) * (r_i - L);
		V += temp; // 													akumulacja V
		
		// F - (14) i dodac do F_i
		Vec3 tempF = r_i < L ? 0 : f*(L-r_i) * rVectors[ii] / r_i;
		F_i[ii] = F_i[ii] + tempF;
			
		// Cisnienie chwilowe (15)
		
		P += tempF.getLength();
		if(ii>0){
			// WEWNETRZNA PETLA
			// <!!!!!!!!>
			for(unsigned int jj=0;jj<ii;jj++){
				Vec3 len;
				len.setR(rVectors[ii]);
				len = len - rVectors[jj];
				double length = len.getLength();
				// pot. par atomowych (9)
				double V_p = e * (pow((R / length),12.0) - 2.0 *  pow((R / length),6.0));
				V += V_p;// 											akumulacja V
				
				// Sily miedzyatomowe Fi Fj (13)
				Vec3 F_ij;
				F_ij = 12.0 * e *((pow(R/length,12.0) - pow(R/length,6.0))) * ((rVectors[ii] - rVectors[jj]) / (length * length));
				F_i[ii] =  F_i[ii] + F_ij;
				F_i[jj] =  F_i[jj] + F_ij * -1.0;
				
			}
			// <!!!!!!!!>
		}
	}
	P *= 1.0 / (4.0 * PI * L * L);
}

int main(int argc, char** argv){
	std::ofstream xyz("xyz.dat");
	std::ofstream out("out.dat");
	srand48(time(NULL));
	Vec3 b0 (a, 0.0, 0.0);
	Vec3 b1 (a / 2.0, a / 2.0 * sqrt(3), 0.0);
	Vec3 b2 (a / 2.0, a / 6.0 * sqrt(3), a * sqrt(2.0/ 3.0));
	
	ofstream plik("dane");
	
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
	std::vector<Vec3> F_i;
	for(unsigned int i=0;i<N;i++){
		F_i.push_back(Vec3(0,0,0));
	}
	
	calculateForces(F_i, rVectors, V, P);
	for(unsigned int i=0;i<N;i++){
				std::cout<<F_i[i]<<std::endl;
	}
	for(unsigned int s = 1;s< STEPS +STEPS_O; s++){
		
		for(unsigned int i=0;i<N;i++){
			pVectors[i] = pVectors[i] + 0.5 * F_i[i] * dt;
			rVectors[i] = rVectors[i] + 1.0 / m * pVectors[i] * dt;	
		}
		
		calculateForces(F_i, rVectors, V, P);
		
		for(unsigned int i=0;i<N;i++){
			pVectors[i] = pVectors[i] + 0.5 * F_i[i] * dt;	
		}
		
		if(s % STEPS_OUT == 0){
			double E_kin = 0;
			double T = 0;
			for(unsigned int i=0;i<N;i++){
				E_kin += pVectors[i].getLengthSquared() / (2.0 * m);
			}
			T = 2.0 / (3.0 * N * k) * E_kin;
			
			out<<s*dt<<"\t"<<V<<"\t"<<E_kin<<"\t"<<E_kin+V<<"\t"<<T<<"\t"<<P<<"\t"<<std::endl;
		}
		if(s % STEPS_XYZ == 0){
			for(unsigned int i=0;i<N;i++){
				xyz<<rVectors[i]<<std::endl;
			}
			xyz<<std::endl<<std::endl;
		}
		
		if(s % STEPS_O){
			//std::cout<< (double) s / ((double) STEPS +STEPS_O) * 100.0<<std::endl;
		}
		
	}
	
	plik.close();
	
	out.close();
	
}
