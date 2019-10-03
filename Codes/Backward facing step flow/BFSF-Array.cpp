#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

#define C_ERROR 0.001
#define Re 800
#define gridi 3001
#define gridj 101
#define MAX(X,Y) (X>Y)?X:Y
using namespace std;

ofstream File1,File2;
int i,j,iterpsi,iteromega;
long double dt, dx, dy, beta, U, h;
long double error_W,error_S;

double psi[gridi][gridj], omega[gridi][gridj], u[gridi][gridj], v[gridi][gridj], temppsi[gridi][gridj], tempomega[gridi][gridj];

double a[gridi], b[gridi], c[gridi], d[gridi], p[gridi], q[gridi];

void updateWalls(){
double y=0;

	//Left Wall
	for(j=(gridj-1)/2;j<=gridj-1;j++){
		
		y=j*((double)1/(gridj-1)) - 0.5;
		
		u[0][j]=12*y*(1-2*y);
		psi[0][j]=2*y*y*(3-4*y);
		omega[0][j]=12*(4*y-1);
		
		
	}

	//Right Wall
	for(j=0;j<=gridj-1;j++){
		
		y=j*((double)1/(gridj-1)) - 0.5;
		
		u[gridi-1][j]=(3*(1-4*y*y))/4;
		psi[gridi-1][j]=(1+3*y-4*y*y*y)/4;
		omega[gridi-1][j]=6*y;
	}
}

void Initialise(){
	//Change
	dx=(double)30/(gridi-1);
	dy=(double)1/(gridj-1);

	beta=dx/dy;
	dt=0.001;

	
	for(i=0;i<gridi;i++)
		for(j=0;j<gridj;j++){
			omega[i][j]=0;
			psi[i][j]=0;
			u[i][j]=0;
			v[i][j]=0;
			temppsi[i][j]=0;
			tempomega[i][j]=0;
		}
	
	updateWalls();

	for(i=0;i<gridi;i++){
		psi[i][gridj-1]=0.5;
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

void CopyTwoinOne(double a[gridi][gridj], double b[gridi][gridj]){
	for(uint i=0;i<gridi;i++){
		for(uint j=0;j<gridj;j++)
			a[i][j]=b[i][j];
	}
}




void WriteFile(string s,int iteromega){
	ofstream File(s,ios::app | ios::out);
	File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\",\"omega\",\"u\",\"v\"";
	
	File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<gridi<<",J="<<gridj<<",\tF=POINT\n";
	
	for(uint j=0;j<gridj;j++){
		for(uint i=0;i<gridi;i++)
			File<<(i)*dx<<" "<<(j)*dy<<" "<<psi[i][j]<<" "<<omega[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<endl;
	}
	
	File.close();
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
				}
		}while(iterpsi<=30);

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

		CopyTwoinOne(tempomega,omega);

		if(iteromega%5000==0){
		
		cout<<"Error in ("<<iteromega<<") Omega = "<<error_W<<endl;
			
			ofstream File("PSI-Plot.plt",ios::app | ios::out);
			File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			
			File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<gridi<<",J="<<gridj<<",\tF=POINT\n";
			
			for(uint j=0;j<gridj;j++){
				for(uint i=0;i<gridi;i++)
					File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<endl;
			}
			
			File.close();
				
		}


	}while(iteromega<1000 || error_W>C_ERROR);

			ofstream File("PSI-Plot.plt",ios::app | ios::out);
			File<<"TITLE=\"2D\"\nVARIABLES=\"X\",\"Y\",\"psi\"";
			
			File<<"\nZONE\tT=BLOCK"<<iteromega<<"\tI="<<gridi<<",J="<<gridj<<",\tF=POINT\n";
			
			for(uint j=0;j<gridj;j++){
				for(uint i=0;i<gridi;i++)
					File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<endl;
			}
			
			File.close();
			WriteFile("PSI-OMEGA-U-V.plt",iteromega);


}

int main(){
	Initialise();
	Solve(0,gridi,0,gridj);
	
}











