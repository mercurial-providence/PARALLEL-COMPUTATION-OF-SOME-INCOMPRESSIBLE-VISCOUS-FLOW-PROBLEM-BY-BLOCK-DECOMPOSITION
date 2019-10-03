#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

#define C_ERROR 0.001
#define Re 400
#define gridi 100
#define gridj 100
#define MAX(X,Y) (X>Y)?X:Y
using namespace std;

ofstream File1,File2;
int i,j,iterpsi,iteromega;
long double dt, dx, dy, beta, U, h;
long double error_W,error_S;

vector< vector<double> > psi(gridi, vector<double>(gridj,0));
vector< vector<double> > omega(gridi, vector<double>(gridj,0));
vector< vector<double> > u(gridi, vector<double>(gridj,0));
vector< vector<double> > v(gridi, vector<double>(gridj,0));
vector< vector<double> > temppsi(gridi, vector<double>(gridj,0));
vector< vector<double> > tempomega(gridi, vector<double>(gridj,0));


vector<double> a(MAX(gridi,gridj),0);
vector<double> b(MAX(gridi,gridj),0);
vector<double> c(MAX(gridi,gridj),0);
vector<double> d(MAX(gridi,gridj),0);
vector<double> p(MAX(gridi,gridj),0);
vector<double> q(MAX(gridi,gridj),0);

void Initialise(){
	//Change
	dx=(double)1/(gridi-1);
	dy=(double)1/(gridj-1);

	beta=dx/dy;
	dt=0.001;

	for(i=1;i<=gridi-2;i++){
		u[i][gridj-1]=1.0;
		

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

//Slower Conversion, @ 0.8
		for(i=1;i<=gridi-2;i++){
			for(j=1;j<=gridj-2;j++){
				if(i==1 || i==gridi-2 || j==1 || j==gridj-2){
					omega[0][j]=(v[1][j]-v[0][j])/dx;
					omega[gridi-1][j]=(v[gridi-1][j]-v[gridi-2][j])/dx;
					omega[i][0]=-(u[i][1]-u[i][0])/dy;
					omega[i][gridj-1]=-(u[i][gridj-1]-u[i][gridj-2])/dy;
				}
			}
		}

/*
	for(i=1;i<grid-1;i++){
		omega[0][i]=-(2*psi[1][i])/(dx*dx);
		omega[gridi-1][i]=-(2*psi[gridi-2][i])/(dx*dx);
		omega[i][0]=-(2*psi[i][1])/(dy*dy);
		omega[i][gridj-1]=-(2*psi[i][gridj-2]+2*dy*U)/(dy*dy);
	}
*/	
}

void Solve(uint is, uint ie, uint js, uint je){


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
			for(i=is+1;i<=ie-2;i++)
				for(j=js+1;j<=je-2;j++){
					psi[i][j]=(omega[i][j]*(dx*dx)*(dy*dy)+(psi[i+1][j]+psi[i-1][j])*(dy*dy)+(psi[i][j+1]+psi[i][j-1])*(dx*dx) )/(2.0*(dx*dx)+2.0*(dy*dy));
					error_S=error_S+fabs(psi[i][j]-temppsi[i][j]);
				}
		}while(iterpsi<=30);
		temppsi=psi;

		for(i=is+1;i<=ie-2;i++)
			for(j=js+1;j<=je-2;j++){
				u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2.0*dy);
				v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2.0*dx);
			}
		
		if(iteromega%2==0){
			//X Sweep
			for(j=js+1;j<=je-2;j++){
				for(i=is+1;i<=ie-2;i++){
					a[i]=-(u[i][j]*dt)/(2*dx)-dt/(Re*dx*dx);
					b[i]=1.0+(2*dt)/(Re*dx*dx);
					c[i]=(u[i][j]*dt)/(2*dx)-dt/(Re*dx*dx);
					d[i]=omega[i][j]*(-2*dt/(Re*dy*dy)+1)+omega[i][j-1]*((v[i][j]*dt)/(2*dy)+dt/(Re*dy*dy))+omega[i][j+1]*((-v[i][j]*dt)/(2*dy)+dt/(Re*dy*dy));

				}
				d[is+1]=d[is+1]-a[is+1]*omega[is][j]; a[is+1]=0.0;
				d[ie-2]=d[ie-2]-c[ie-2]*omega[ie-1][j];c[ie-2]=0.0;

				b[is]=1.0;a[is]=0.0;c[is]=0.0;d[is]=omega[is][j];
				p[is]=1.0; q[is]=d[is]/b[is];

				for(i=is+1;i<=ie-2;i++){
					p[i]=b[i]-(a[i]*c[i-1] )/p[i-1];
					q[i]=d[i]-(q[i-1]*a[i])/p[i-1];
				}

				for(i=ie-2;i>=is+1;i--){
					omega[i][j]=(q[i]-c[i]*omega[i+1][j])/p[i];
				}
			}
		}
		else{
			//Y Sweep
			for(i=is+1;i<=ie-2;i++){
				for(j=js+1;j<=je-2;j++){
					a[j]=(-v[i][j]*dt)/(2*dy)-(dt/(Re*dy*dy));
					b[j]=(1+(2*dt)/(Re*dy*dy));
					c[j]=(v[i][j]*dt)/(2*dy)-dt/(Re*dy*dy);
					d[j]=omega[i][j]*((1-(2*dt)/(Re*dx*dx)))+omega[i-1][j]*(u[i][j]*dt/(2*dx)+dt/(Re*dx*dx))
						+omega[i+1][j]*((-u[i][j]*dt)/(2*dx)+dt/(Re*dx*dx));

				}
				d[js+1]=d[js+1]-a[js+1]*omega[i][js]; a[js+1]=0.0;
				d[je-2]=d[je-2]-c[je-2]*omega[i][je-1];c[je-2]=0.0;

				a[js]=0.0;b[js]=1.0;c[js]=0.0;d[js]=omega[i][js];
				p[js]=1.0; q[js]=d[js]/b[js];

				for(j=js+1;j<=je-2;j++){
					p[j]=b[j]-(a[j]*c[j-1])/p[j-1];
					q[j]=d[j]-(q[j-1]*a[j])/p[j-1];
				}

				for(j=je-2;j>=js+1;j--){
					omega[i][j]=(q[j]-c[j]*omega[i][j+1])/p[j];
				}
			}
		}

		for(i=is+1;i<=ie-2;i++)
			for(j=js+1;j<=je-2;j++)
				error_W=error_W+fabs(omega[i][j]-tempomega[i][j]);

		tempomega=omega;

		cout<<"Error in Omega = "<<error_W<<endl;

		if(iteromega%100==0){
			File1.open("psi.plt", ios::out);
			File2.open("omega.plt", ios::out);

			File1<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			File2<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"omega\"\n";

			File1<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<ie<<",J="<<je<<",\tF=POINT\n";
			File2<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<ie<<",J="<<je<<",\tF=POINT\n";

			for(j=js;j<=je-1;j++)
				for(i=is;i<=ie-1;i++){
					File1<<i*dx<<"\t"<<j*dx<<"\t"<<psi[i][j]<<endl;
					File2<<i*dx<<"\t"<<j*dx<<"\t"<<omega[i][j]<<endl;
				}
			File1.close();
			File2.close();
		}


	}while(error_W>0.001 || iteromega<1000);

			File1.open("psi.plt", ios::out);
			File2.open("omega.plt", ios::out);

			File1<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			File2<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"omega\"\n";

			File1<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<ie<<",J="<<je<<",\tF=POINT\n";
			File2<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<ie<<",J="<<je<<",\tF=POINT\n";

			for(j=js;j<=je-1;j++)
				for(i=is;i<=ie-1;i++){
					File1<<i*dx<<"\t"<<j*dx<<"\t"<<psi[i][j]<<endl;
					File2<<i*dx<<"\t"<<j*dx<<"\t"<<omega[i][j]<<endl;
				}
			File1.close();
			File2.close();


}

int main(){
	Initialise();
	Solve(0,gridi,0,gridj);
	
}











