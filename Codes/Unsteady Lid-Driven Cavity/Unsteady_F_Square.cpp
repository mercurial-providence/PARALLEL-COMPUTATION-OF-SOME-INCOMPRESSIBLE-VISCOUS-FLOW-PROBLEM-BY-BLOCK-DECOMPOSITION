#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

#define C_ERROR 0.000001
//#define Re 100
#define grid 129

using namespace std;

ofstream File;
int i,j,iterpsi,iteromega;
long double dt, dx, dy, beta, U, h;
long double error_W,error_S;

vector< vector<double> > psi(grid, vector<double>(grid,0));
vector< vector<double> > omega(grid, vector<double>(grid,0));
vector< vector<double> > u(grid, vector<double>(grid,0));
vector< vector<double> > v(grid, vector<double>(grid,0));
vector< vector<double> > temppsi(grid, vector<double>(grid,0));
vector< vector<double> > tempomega(grid, vector<double>(grid,0));

vector<double> a(grid,0);
vector<double> b(grid,0);
vector<double> c(grid,0);
vector<double> d(grid,0);
vector<double> p(grid,0);
vector<double> q(grid,0);

void Initialise(){
	//Change
	dx=(double)1/(grid-1);
	dy=(double)1/(grid-1);

	beta=dx/dy;
	U=1;
	h=(double)1/(grid-1);
	dt=0.001;

	for(i=1;i<=grid-2;i++){
		u[i][grid-1]=1.0;
	}
}

void Print(const vector< vector<double> > &v){
	for(uint i=0;i<v.size();i++){
		for(uint j=0;j<v[i].size();j++)
			cout<<"["<<i<<"]"<<"["<<j<<"]"<<v[i][j]<<"\t";
		cout<<endl;
	}
	cout<<endl;
}

void PrintSingle(const vector< double > &v){
	for(uint i=0;i<v.size();i++){
		cout<<"["<<i<<"]"<<v[i]<<"\t";
	}
	cout<<endl<<endl;
}

void SetBoundaryConditions(){
/*
//Slower Conversion, @ 0.8
		for(i=1;i<=grid-2;i++){
			for(j=1;j<=grid-2;j++){
				if(i==1 || i==grid-2 || j==1 || j==grid-2){
					omega[0][j]=(v[1][j]-v[0][j])/dx;
					omega[grid-1][j]=(v[grid-1][j]-v[grid-2][j])/dx;
					omega[i][0]=-(u[i][1]-u[i][0])/dy;
					omega[i][grid-1]=-(u[i][grid-1]-u[i][grid-2])/dy;
				}
			}
		}
*/
	for(i=1;i<grid-1;i++){
		omega[0][i]=-(2*psi[1][i])/(dx*dx);
		omega[grid-1][i]=-(2*psi[grid-2][i])/(dx*dx);
		omega[i][0]=-(2*psi[i][1])/(dy*dy);
		omega[i][grid-1]=-(2*psi[i][grid-2]+2*dy*U)/(dy*dy);
	}
	
}

void WriteFile(string s,int iteromega){
	ofstream File(s,ios::app | ios::out);
	File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\",\"omega\",\"u\",\"v\"";
	
	File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<grid<<",J="<<grid<<",\tF=POINT\n";
	
	for(uint i=0;i<grid;i++){
		for(uint j=0;j<grid;j++)
			File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<" "<<omega[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<endl;
	}
	
	File.close();
}


void Solve(int Re){

	iteromega=0;
	error_W=100.0;


	do{
		iteromega++;
		error_W=0.0;

		SetBoundaryConditions();

		error_S=100.00;
		iterpsi=0;
		do{
			iterpsi++;
			error_S=0.00;
			for(i=1;i<=grid-2;i++)
				for(j=1;j<=grid-2;j++){
					psi[i][j]=(omega[i][j]*(dx*dx)*(dy*dy)+(psi[i+1][j]+psi[i-1][j])*(dy*dy)+(psi[i][j+1]+psi[i][j-1])*(dx*dx) )/(2.0*(dx*dx)+2.0*(dy*dy));
					error_S=error_S+fabs(psi[i][j]-temppsi[i][j]);
				}
		}while(iterpsi<=30);
		temppsi=psi;

		for(i=1;i<=grid-2;i++)
			for(j=1;j<=grid-2;j++){
				u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2.0*dy);
				v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2.0*dx);
			}
		
		if(iteromega%2==0){
			//X Sweep
			for(j=1;j<=grid-2;j++){
				for(i=1;i<=grid-2;i++){
					a[i]=-(u[i][j]*dt)/(2*dx)-dt/(Re*dx*dx);
					b[i]=1.0+(2*dt)/(Re*dx*dx);
					c[i]=(u[i][j]*dt)/(2*dx)-dt/(Re*dx*dx);
					d[i]=omega[i][j]*(-2*dt/(Re*dy*dy)+1)+omega[i][j-1]*((v[i][j]*dt)/(2*dy)+dt/(Re*dy*dy))+omega[i][j+1]*((-v[i][j]*dt)/(2*dy)+dt/(Re*dy*dy));
				}
				
				d[1]=d[1]-a[1]*omega[0][j]; a[1]=0.0;
				d[grid-2]=d[grid-2]-c[grid-2]*omega[grid-1][j];c[grid-2]=0.0;

				b[0]=1.0;a[0]=0.0;c[0]=0.0;d[0]=omega[0][j];
				p[0]=1.0; q[0]=d[0]/b[0];

				for(i=1;i<=grid-2;i++){
					p[i]=b[i]-(a[i]*c[i-1] )/p[i-1];
					q[i]=d[i]-(q[i-1]*a[i])/p[i-1];
				}

				for(i=grid-2;i>=1;i--){
					omega[i][j]=(q[i]-c[i]*omega[i+1][j])/p[i];
				}
			}
		}
		else{
			//Y Sweep
			for(i=1;i<=grid-2;i++){
				for(j=1;j<=grid-2;j++){
					a[j]=(-v[i][j]*dt)/(2*dy)-(dt/(Re*dy*dy));
					b[j]=(1+(2*dt)/(Re*dy*dy));
					c[j]=(v[i][j]*dt)/(2*dy)-dt/(Re*dy*dy);
					d[j]=omega[i][j]*((1-(2*dt)/(Re*dx*dx)))+omega[i-1][j]*(u[i][j]*dt/(2*dx)+dt/(Re*dx*dx))
						+omega[i+1][j]*((-u[i][j]*dt)/(2*dx)+dt/(Re*dx*dx));
				}
				d[1]=d[1]-a[1]*omega[i][0]; a[1]=0.0;
				d[grid-2]=d[grid-2]-c[grid-2]*omega[i][grid-1];c[grid-2]=0.0;

				a[0]=0.0;b[0]=1.0;c[0]=0.0;d[0]=omega[i][0];
				p[0]=1.0; q[0]=d[0]/b[1];

				for(j=1;j<=grid-2;j++){
					p[j]=b[j]-(a[j]*c[j-1])/p[j-1];
					q[j]=d[j]-(q[j-1]*a[j])/p[j-1];
				}

				for(j=grid-2;j>=1;j--){
					omega[i][j]=(q[j]-c[j]*omega[i][j+1])/p[j];
				}
			}
		}

		for(i=1;i<=grid-2;i++)
			for(j=1;j<=grid-2;j++)
				error_W=error_W+fabs(omega[i][j]-tempomega[i][j]);

		tempomega=omega;

		cout<<"Error in Omega = "<<error_W<<endl;

		if(iteromega%100==0){
		
			ofstream File("UnsteadyRe400ONLYPSI.plt",ios::app | ios::out);
			File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			
			File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<grid<<",J="<<grid<<",\tF=POINT\n";
			
			for(uint i=0;i<grid;i++){
				for(uint j=0;j<grid;j++)
					File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<endl;
			}
			
			File.close();
//			WriteFile("UnsteadyRe400.plt",iteromega);
		}



	}while(error_W>0.01 || iteromega<1000);

			ofstream File("UnsteadyRe400ONLYPSI.plt",ios::app | ios::out);
			File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			
			File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<grid<<",J="<<grid<<",\tF=POINT\n";
			
			for(uint i=0;i<grid;i++){
				for(uint j=0;j<grid;j++)
					File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<endl;
			}
			
			File.close();
	
			WriteFile("UnsteadyRe400.plt",iteromega);

}

int main(){
	Initialise();
	Solve(400);

	
}











