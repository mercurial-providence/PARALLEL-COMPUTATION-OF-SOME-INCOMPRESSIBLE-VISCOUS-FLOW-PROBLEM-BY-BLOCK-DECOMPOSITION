/*
Compile:

g++ filename.cpp -std=c++0x

*/

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<cmath>
#include<fstream>


/*
 * Defining constants used in code
 */
#define C_ERROR 0.0000001
//#define Re 700
//#define rel 0.3
#define grid 257

using namespace std;

/*
 * Defining Global Variables
 */
vector< vector<double> > psi;
vector< vector<double> > omega;
vector< vector<double> > u;
vector< vector<double> > v;

double dx=0, dy=0, beta=0, U=0;

/*
 * This function allocates size and values of Global variables
 */
void Initialise(){
	psi.resize(grid, vector<double>(grid,0));
	omega.resize(grid, vector<double>(grid,0));
	u.resize(grid, vector<double>(grid,0));
	v.resize(grid, vector<double>(grid,0));

	for(uint i=0;i<u[0].size();i++)
		u[0][i]=1;

	dx=dy=(double)1/grid;
	beta=dx/dy;
	U=1;
}

/*
 * This function can be used to Print any Matrix/2D-Vector variable
 */
void Print(vector< vector<double> > &v){
	for(uint i=0;i<v.size();i++){
		for(uint j=0;j<v[i].size();j++)
			cout<<"["<<i<<"]"<<"["<<j<<"]"<<v[i][j]<<"\t";
		cout<<endl;
	}
}

/*
 * This function updates the boundary conditions on each call
 */
void SetBoundaryConditions(){
	uint i=0;

	/*
	//Totally accurate conversion. Do not doubt it again.
	for(i=1;i<omega.size()-1;i++){
		omega[0][i]=-1*(2*psi[1][i]+2*dy*U)/(dy*dy); 	//TOP WALL
		omega[i][grid-1]=-1*(2*psi[i][grid-2])/(dx*dx); //RIGHT WALL
		omega[grid-1][i]=-1*(2*psi[grid-2][i])/(dy*dy); //BOTTOM WALL
		omega[i][0]=-1*(2*psi[i][1])/(dx*dx);			//LEFT WALL
	}
	*/

	for(i=1;i<grid-1;i++){
		omega[0][i]=-(2*psi[1][i])/(dx*dx);
		omega[grid-1][i]=-(2*psi[grid-2][i])/(dx*dx);
		omega[i][0]=-(2*psi[i][1])/(dy*dy);
		omega[i][grid-1]=-(2*psi[i][grid-2]+2*dy*U)/(dy*dy);
	}


}

/*
 * This function calculate the Error and returns it
 */
double Error(vector< vector<double> > a, vector< vector<double> > b){
	double diffSum=0, nSum=0;
	for(uint i=0;i<a.size();i++)
		for(uint j=0;j<b.size();j++) {
		nSum+=fabs(a[i][j]);
		diffSum+=fabs(a[i][j]-b[i][j]);
	}
	return diffSum/nSum;
}

/*
 * This is the main iterator. This function take the Reynolds number as input and computer psi and omega.
 */
void Solve(uint Re,double rel){
	uint i=0, j=0, iter=0;
	double error=0;
	vector< vector<double> > temp;


	do{
	temp=omega;

	SetBoundaryConditions();



	//Computing PSI
	for(uint k=0;k<3;k++)
	for(i=1;i<psi.size()-1;i++)
		for(j=1;j<psi[i].size()-1;j++)
			psi[i][j]=((psi[i+1][j]+psi[i-1][j]+beta*beta*(psi[i][j+1]+psi[i][j-1])+(dx*dx)*omega[i][j]))/(2*(1+beta*beta));

	//Computing OMEGA
	for(uint k=0;k<3;k++)
	for(i=1;i<omega.size()-1;i++)
		for(j=1;j<omega[i].size()-1;j++)
			omega[i][j]=(1-rel)*omega[i][j]+rel*(omega[i+1][j]+omega[i-1][j]+beta*beta*(omega[i][j+1]+omega[i][j-1])-0.25*beta*Re*(omega[i+1][j]-omega[i-1][j])*(psi[i][j+1]-psi[i][j-1])+0.25*beta*Re*(omega[i][j+1]-omega[i][j-1])*(psi[i+1][j]-psi[i-1][j]))/(2*(1+(beta*beta)));



	error=Error(omega,temp);
	if(++iter%100==0) cout<<"Iteration : "<<iter<<"\t Error : "<<error<<endl;
	}while(error>C_ERROR);

	//Solving u and v
	for(i=1;i<grid-1;i++){
		for(j=1;j<grid-1;j++){
			u[i][j]=-(psi[i][j+1]-psi[i][j-1])/(2*dy);
			v[i][j]=(psi[i+1][j]-psi[i-1][j])/(2*dy);
		}
	}

}

/*
 * This function writes the output data in File
 */
void WriteFile(string s){
	ofstream File(s,ios::app | ios::out);
	for(uint i=0;i<grid;i++){
		for(uint j=0;j<grid;j++)
			File<<(i+1)*dx<<" "<<(j+1)*dy<<" "<<psi[i][j]<<" "<<omega[i][j]<<" "<<u[i][j]<<" "<<v[i][j]<<endl;
	}
}

/*
 * Required main function controlling and calling all functions.
 */
int main(){
	Initialise();


	cout<<"Switching to Re=100\n";
	Solve(100,1);
	WriteFile("Re100.plt");
	
	cout<<"Switching to Re=200\n";
	Solve(200,1);
	WriteFile("Re200.plt");
		
	cout<<"Switching to Re=400\n";
	Solve(400,1);
	WriteFile("Re400.plt");
	
	cout<<"Switching to Re=800\n";
	Solve(800,0.8);
	WriteFile("Re800.plt");
	
	cout<<"Switching to Re=2000\n";
	Solve(2000,0.3);
	WriteFile("Re2000.plt");
	
	


}



