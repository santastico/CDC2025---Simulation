#ifndef _WIN32
#include "MinMaxGD-main/include/lminmaxgd.h"
#else
#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\include\lminmaxgd.h"
using namespace std;
#endif
#include <time.h>


// This code provide you the tool to obtain Fopt and Popt
/////////////////////////////////////////////////////////////////////
/* Example : Discrete Event Systems in Dioid Framework : Control Theory
 Laurent Hardouin, Olivier Boutin, Bertrand Cottenceau, Thomas Brunsch, Joerg Raisch
 Springer, DISC Book
This program compute :
1. the optimal control when a reference output is a priori known
2. the optimal input filtering when the reference output is unknown
3. the optimal feedback
*/

#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\src\gd.cpp"
#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\src\poly.cpp"
#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\src\serie.cpp"
#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\src\smatrix.cpp"
#include "C:\Users\santo\Desktop\GET_BO\MinMaxGD-main\src\tools.cpp"

using namespace std;
using namespace mmgd;



int main(void)
{
try{


//This example is from the link http://www.istia.univ-angers.fr/~hardouin/GET_incertain.html
//by adding disturbance signals and only consider the black token, i.e. the lower bound case
int n=4;
// System Matrices
smatrix A(n,n), B(n,2), C(1,n);
smatrix M(n,n), N(n,n);
smatrix Astar;
smatrix CA,CAB,AB;
smatrix P0;
serie s;
poly p,q;
gd r;
int i,j;




// defining A=A_bar= A * \gamma in CDC11 paper
    A(0,1)=gd(2,2);
   A(1,3)=gd(0,4);
    A(1,0)=gd(0,3);

   A(2,3)=gd(3,7);

     A(3,1)=gd(2,2);
   A(3,2)=gd(0,7);

    cout << "A :" << A<< endl;

    //defining B matrix

    B(0,0)=gd(0,1);
     B(2,1)=gd(0,1);

    cout << "B :" << B<< endl;

    //defining C matrix
    C(0,1)=gd(0,1);
    cout << "C :" << C<< endl;






// Calculation of the transfer function matrix

Astar = star(A);

cout<<"Astar\n"<<Astar<<endl;
CA=otimes(C,Astar);
AB=otimes(Astar,B);
CAB=otimes(C,AB);

cout<<"AB\n"<<AB<<endl;

P0=lfrac(CAB,CAB);
cout<<"P0\n"<<P0<<endl;

smatrix ABP0=otimes(AB,P0);
cout<<"ABP0 : "<<ABP0<<endl;



for(int i=0;i<n;i++)
{

    M(i,i)=gd(0,0);
    N(i,i)=gd(0,0);
}
 /* q.init(0,6);
 r.init(3,14);
 s.init(epsilon,q,r);
 N(1,0)=s;
*/
 M(0,1)=gd(0,-6);
 M(2,1)=gd(0,-12);




 cout<<"M!!!!"<<M<<endl;
 cout<<"N!!!!"<<N<<endl;

smatrix MAB=otimes(M,AB);
smatrix NAB=otimes(N,AB);

     cout<<"MAB"<<MAB<<endl;
     cout<<"NAB"<<NAB<<endl;



smatrix Pn1=prcaus(P0);
smatrix Pn,CABPn,ABPn,MABPn,NABPn;
smatrix temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8;

cout<<"P0+"<<Pn1<<endl;
/*
smatrix P1(2,2),ABP1,MABP1,NABP1;
i=0;
for(i=0;i<2;i++)
{
    P1(i,i)=e;
}
ABP1=otimes(AB,P1);
MABP1=otimes(M,ABP1);
     //cout<<"3"<<" i="<<i<<"MABPn"<<MABPn<<endl;
     NABP1=otimes(N,ABP1);
   //  cout<<"4"<<" i="<<i<<"NABPn"<<NABPn<<endl;
     temp=lfrac(NABP1,MAB);
     //cout<<"5"<<" i="<<i<<" NABPn/ABPn="<<temp<<endl;
     temp1=lfrac(MABP1,NAB);
    // cout<<"6"<<" i="<<i<<" MABPn/ABPn="<<temp1<<endl;
     Pn1=inf(temp1,temp);

     //We come back */
     Pn1=P0;
do
{    Pn=Pn1;
 //    cout<<"1"<<" i="<<i<<"Pn"<<Pn<<endl;
     ABPn=otimes(AB,Pn);
   //  cout<<"2"<<" i="<<i<<"ABPn"<<ABPn<<endl;
     MABPn=otimes(M,ABPn);
     //cout<<"3"<<" i="<<i<<"MABPn"<<MABPn<<endl;
     NABPn=otimes(N,ABPn);
   //  cout<<"4"<<" i="<<i<<"NABPn"<<NABPn<<endl;
     temp=lfrac(NABPn,MAB);
     //cout<<"5"<<" i="<<i<<" NABPn/ABPn="<<temp<<endl;
     temp1=lfrac(MABPn,NAB);
    // cout<<"6"<<" i="<<i<<" MABPn/ABPn="<<temp1<<endl;
     temp2=inf(temp1,temp);
    // cout<<"7"<<" i="<<i<<" MABPn/ABPn ^ NABPn/ABPn="<<temp2<<endl;
     Pn1=inf(Pn,temp2);

     Pn1=prcaus(Pn1);
    // cout<<"8"<<" i="<<i<<" Pn1 "<<Pn1<<endl;




    if(Pn1==Pn)
    {
                cout<<"convergence achieved Pn is the greatest smaller than P0 such that it exist a Popt solving constraint"<<endl;
                cout<<"Pn"<<Pn<<endl;
                cout<<"P0"<<P0<<endl;
                ABPn=otimes(AB,Pn);

                MABPn=otimes(M,ABPn);
                NABPn=otimes(N,ABPn);
                if(MABPn==NABPn)
                { cout<<"This Pn respects the constraints  \n"<<Pn<<endl;
                  cout<<" with ABPn \n"<<ABPn<<endl;
                  cout<<" with ABP0 \n"<<ABP0<<endl;

                }
                else
                {
                    cout<<"This Pn doesn't respect the constraints  \n"<<Pn<<endl;

                  cout<<" with ABPn \n"<<ABPn<<endl;

                  cout<<"with MABPn \n"<<MABPn<<endl;
                   cout<<"with NABPn \n"<<MABPn<<endl;

                }

    }
    else
    {
        i++;


         cout<<"This Pn1 doesn't respect the constraints  \n"<<Pn1<<endl;


    }

}while(i<100 && !(Pn1==Pn));

smatrix Fopt;
temp=lfrac(Pn,Pn);
temp1=otimes(CAB,Pn);
Fopt=rfrac(temp,temp1);
Fopt=prcaus(Fopt);
cout<<"Fopt"<<Fopt<<endl;

cout<<"AB "<<AB<<endl;
cout<<"ABPn "<<ABPn<<endl;
temp1=AB(2,0);
temp2=AB(1,0);
temp=lfrac(temp1,temp2);
cout<<"Tau1 "<<temp<<endl;
temp1=AB(2,1);
temp2=AB(1,1);
temp=lfrac(temp1,temp2);
cout<<"Tau2 "<<temp<<endl;


cout<<"ABPn "<<ABPn<<endl;
temp1=ABPn(2,0);
temp2=ABPn(1,0);
temp=lfrac(temp1,temp2);
cout<<"Tau1Pn "<<temp<<endl;
temp1=ABPn(2,1);
temp2=ABPn(1,1);
temp=lfrac(temp1,temp2);
cout<<"Tau2Pn "<<temp<<endl;


smatrix CrossABAB;
CrossABAB=rfrac(AB,AB);
cout<<"CrossABAB "<<CrossABAB<<endl;

}



   catch(mem_limite l)
 {
	 cout<<"Exception : too many coefficent in polynom "<<l.memoire<<endl;
	 return(1);
 }

 catch(taille_incorrecte obj)
 { // 0 : r non causal
   // 1 : tentative d'accès à un element d'une matrice avec un indice incorrect
   // 2 : matrice de taille incompatible pour oplus, inf, otimes, rfrac, lfrac
   // 3 : etoile de matrice carré uniquement
	 cout<<"Exception  "<<obj.erreur<<endl;
	 return(1);
 }

}
