// AnalyzeFullC.cpp -- Analyze C matrix from MC_WangLandau.hpp
//
// October 23, 2014
// Greg Brown (gbrown@fsu.edu,browngrg@comcast.net)


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

// DGEBAL - balance a general real matrix A
#ifdef USE_BLAS
extern "C" void dgeev_(char* jobvl, char* jobvr, int* N, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int *lwork, int* info);
extern "C" void dgemm_(char* transa, char* transb ,int* M, int* N, int* K, const double* alpha, const double* A, int* lda, const double* B, int* ldb, const double* beta, double* C, int* ldc);
extern "C" void dgemv_(char* trans, int* M, int* N, const double* alpha, const double* A, int* lda, const double* X, int* incx, const double* beta, double* y, int* incy);
extern "C" void dgetrf_(int* M, int* N, double* A, int* lda, int* ipivot, int* info);
extern "C" void dgetri_(int* N, double* A, int* lda, int* ipivot, double* work, int* lwork, int* info);
#endif


using namespace std;


void multiply(const vector<double>& mat, vector<double>& vec)
{
   int N = vec.size();
   vector<double> m(N,0);
#  if USE_BLAS
   // y := alpha*A*x + beta*y
   char trans = 'N';
   int incr = 1;
   double alpha = 1;
   double beta  = 0;
   dgemv_(&trans,&N,&N,&alpha,&(mat[0]),&N,&(vec[0]),&incr,&beta,&(m[0]),&incr);
#  else
   for(int jcol=0; jcol<N; jcol++)
      for(int irow=0; irow<N; irow++)
         m[irow] += mat[irow+N*jcol]*vec[jcol];
#  endif
   vec = m;
}

// Returns the the Rayleigh quotient of x with respect to A, r(x) = xT*A*x/(xT*x)
double Rayleigh(const vector<double>& A, const vector<double>& x)
{
#  if USE_BLAS
   // This implementation assumes xT*x = 1
   int N = x.size();
   char trans = 'N';
   int incr = 1;
   double alpha = 1;
   double beta = 0;
   vector<double> Ax(N,0);
   dgemv_(&trans,&N,&N,&alpha,&(A[0]),&N,&(x[0]),&incr,&beta,&(Ax[0]),&incr);
   double xTAx = 0;
   double xTx = 0;
   for(int i=0; i<N; i++) 
   {
      xTAx += x[i]*Ax[i];
      xTx  += x[i]*x[i];
   }
   return xTAx/xTx;
#  else
   int N = x.size();
   double xTAx = 0;
   double xTx  = 0;
   for(int i=0; i<N; i++)
   {
     for(int k=0; i<N; i++)
        xTAx += x[i]*A[i+k*N]*x[k];
     xTx += x[i]*x[i];
   }
   return xTAx/xTx;
#  endif
}


double RayleighQuotientIteration(const vector<double>& mat, vector<double>& evec, int niter=20)
{
#  ifndef USE_BLAS
   std::cout << "RayleighQuotientIteration only works with -DUSE_BLAS" << std::endl;
   return 0;
#  endif
   // Converge to the eigenvector using the Rayleigh quotient method
   // Returns the associated eigenvalue
   // This is known to have problems for non-symmetric matrice
   // Maysum Panju, "Iterativec Methods for Computing Eigenvalues and Eigenvectors,"
   // The Waterloo Mathematics Review 1, 9-18 (2011).
   ofstream fout; if(true) fout.open("RayleighQuotient.csv");
   int N = evec.size();
   double vmag = 0; for(int i=0; i<N; i++) vmag += evec[i]*evec[i];      // Normalize the eigenvector
   vmag = sqrt(vmag);
   for(int i=0; i<N; i++) evec[i] /= vmag;
   double lambda = Rayleigh(mat,evec);
   fout << "#index=" << 0 << " lambda=" << lambda << endl;
   for(int i=0; i<N; i++) fout << i << " " << log(max(1.e-20,fabs(evec[i]))) << " " << evec[i] << endl;
   //cout << "# converge init lambda=" << lambda << std::endl;
   for(int k=0; k<niter; k++)
   {
      // Solve (A-lambda*I)*w = v
      vector<double> matinv = mat;
      for(int i=0; i<N; i++) matinv[i+i*N] -= lambda;                                    // mat = A - lambda*I
      vector<int> ipivot(N);
      int info;
      dgetrf_(&N,&N,&(matinv[0]),&N,&(ipivot[0]),&info);                                 // LU factorization of matrix A=P*L*U
      vector<double> work(10*N);
      int lwork = work.size();
      dgetri_(&N,&(matinv[0]),&N,&(ipivot[0]),&(work[0]),&lwork,&info);                  // Computes inverse matrix of LU factorized matrix
      char trans = 'N';
      int incr = 1;
      double alpha = 1;
      double beta = 0;
      vector<double> w(N);
      dgemv_(&trans,&N,&N,&alpha,&(matinv[0]),&N,&(evec[0]),&incr,&beta,&(w[0]),&incr);  // w = mat^{-1}*v
      double wmag = 0; for(int i=0; i<N; i++) wmag += w[i]*w[i];
      wmag = sqrt(wmag);
      for(int i=0; i<N; i++) evec[i] = w[i]/wmag;
      lambda = Rayleigh(mat,evec);
      if(fout && fout.is_open())
      {
         fout << endl << endl;
         fout << "#index=" << k+1 << " lambda=" << lambda << endl;
         for(int i=0; i<N; i++) fout << i << " " << log(max(1.e-20,fabs(evec[i]))) << " " << evec[i] << endl;
      }
      //cout << "# converge: " << k <<  "  lambda=" << lambda << std::endl;
   }
   if(fout && fout.is_open()) fout.close();
   return lambda;
}


// Normalize dos
void fix(vector<double>& dos)
{
   int nbin = dos.size();
   // starts at zero
   if( dos[0]!=0 ) for(int irow=nbin-1; irow>=0; irow--) dos[irow] -= dos[0];
   // positive
   if( dos[nbin-1]<0 ) for(int irow=0; irow<nbin; irow++) dos[irow] *= -1;
   // L1Norm = 1;
   // double norm=0; 
   // for(int irow=0; irow<nbin; irow++) norm += dos[irow]; 
   // for(int irow=0; irow<nbin; irow++) dos[irow] /= norm;
}



int main(int argc, char* argv[])
{

   // Get basic details
   int nbin;
   cin >> nbin;

   // Read count matrix
   vector<long long> Cread(nbin*nbin);
   for(int jcol=0; jcol<nbin; jcol++)
      for(int irow=0; irow<nbin; irow++)
         cin >> Cread[irow+nbin*jcol];
   
   // look for empty rows and bins
   int minbin = 0;
   long long cts = 0;
   while( minbin<nbin && cts==0 )
   {
      for(int jcol=minbin; jcol<nbin; jcol++)
         cts += Cread[minbin+nbin*jcol];
      minbin++;
   }
   if(cts>0) minbin--;
   cts = 0;
   while( minbin<nbin && cts==0 )
   {
      for(int irow=minbin; irow<nbin; irow++)
         cts += Cread[irow+nbin*minbin];
      minbin++;
   }
   if(cts>0) minbin--;

   // possibly downsize
   int nold = nbin;
   nbin = nbin - minbin;
   vector<long long> C(nbin*nbin);
   for(int jcol=0; jcol<nbin; jcol++)
      for(int irow=0; irow<nbin; irow++)
         C[irow+nbin*jcol] = Cread[(irow+minbin)+nold*(jcol+minbin)];

   // Create normalized transition matrix
   vector<long long> Hsum(nbin);
   vector<double>    Tinf(nbin*nbin);
#  if 1
   // Transpose the transition matrix
   for(int irow=0; irow<nbin; irow++)
      for(int jcol=irow+1; jcol<nbin; jcol++)
         swap(C[irow+nbin*jcol],C[jcol+nbin*irow]);
   // Normalize rows
   for(int irow=0; irow<nbin; irow++)
   {
      Hsum[irow] = 0;
      for(int jcol=0; jcol<nbin; jcol++)
         Hsum[irow] += C[irow+nbin*jcol];
      for(int jcol=0; jcol<nbin; jcol++)
         Tinf[irow+nbin*jcol] = static_cast<double>(C[irow+nbin*jcol])/static_cast<double>(Hsum[irow]);
   }
#  else
   // Normalize columns (should give same results as above)
   for(int jcol=0; jcol<nbin; jcol++)
   {
      Hsum[jcol] = 0;
      for(int irow=0; irow<nbin; irow++)
         Hsum[jcol] += C[irow+nbin*jcol];
      for(int irow=0; irow<nbin; irow++)
         Tinf[irow+nbin*jcol] = static_cast<double>(C[irow+nbin*jcol])/static_cast<double>(Hsum[jcol]);
   }
#  endif

   // Calculate beta from Eq. (3) of PMC de Oliveira, TJP Penna, HJ Herrmann Eur Phys J B 1, 205 (1998)
   // This is idential to lng2 when strictly tri-diagonal
   vector<double> beta(nbin);   // actually deltaE*beta(E)
   vector<double> lng(nbin);    // from all up/down
   lng[0] = 0;
   for(int irow=0; irow<(nbin-1); irow++)
   {
      beta[irow] =  -log(Tinf[irow+nbin*(irow+1)]/Tinf[(irow+1)+nbin*irow]); // minus sign cuz transpose of Tinf in paper
      if( beta[irow]!=beta[irow] ) beta[irow]=0;
      lng[irow+1] = lng[irow] + beta[irow];
   }
   beta[nbin-1] = beta[nbin-2];
   fix(lng);

   // Improve estimate by finding eigenvector of transition matrix
   double gmax = max(lng[0],lng[nbin-1]);
   vector<double> ev(nbin); for(int irow=0; irow<nbin; irow++) ev[irow] = exp(lng[irow]-gmax);
   double lambda = RayleighQuotientIteration(Tinf,ev);

   double lambda_dgeev = 0;
   vector<double> ev_dgeev(nbin,1);
#  ifdef USE_BLAS
   if(true)
   {
      // Find the eigenvectors with a canned library routine
      // Doesn't give particularly good results when try using scipy
      char job   = 'B';            // N=none, S=scale, B=scale and permute
      char side  = 'R';            // R = right eigenvectors, L = left eigenvectors
      int N = nbin;                // extent of matrix
      int lda = nbin;              // leading dimension of matrix
      vector<double> tmat(Tinf);
      int info;
      vector<double> wr(N);  // real and imaginary parts of the eigenvalues
      vector<double> wi(N);
      int ldvl = N;
      int ldvr = N;
      char jobvl = 'V';      // N = eigenvectors are not computed, V = eigenvectors are computed
      char jobvr = 'N';      // N = eigenvectors are not computed, V = eigenvectors are computed
      if( side=='R' )
      {
         jobvl = 'N';
         jobvr = 'V';
      }
      vector<double> vl;     // matrix of left eigenvectors
      vector<double> vr;     // matrix of right eigenvectors
      if( jobvl=='V') vl.resize(ldvl*N);
      if( jobvr=='V') vr.resize(ldvr*N);
      int lwork = 50*N;
      vector<double> work(lwork);
      dgeev_(&jobvl,&jobvr,&N,&(tmat[0]),&lda,&(wr[0]),&(wi[0]),&(vl[0]),&ldvl,&(vr[0]),&ldvr,&(work[0]),&lwork,&info);
      int jcol= 1;
      if( side=='R' )
         for(int irow=0; irow<N; irow++) ev_dgeev[irow] = vr[irow+N*jcol];
      else
         for(int irow=0; irow<N; irow++) ev_dgeev[irow] = vl[irow+N*jcol];
      lambda_dgeev = wr[1];
   }
#  endif

   for(int irow=0; irow<nbin; irow++) ev[irow] = log(fabs(ev[irow]));
   for(int irow=0; irow<nbin; irow++) ev_dgeev[irow] = log(fabs(ev_dgeev[irow]));

   fix(ev);
   fix(ev_dgeev);

   cout << "# Data from Analysis of Full C Matrix" << endl;
   cout << "# orig nbin = " << nold <<" analyze nbin=" << nbin << endl;
   cout << "# lambda = " << lambda << " lambda_dgeev=" << lambda_dgeev << endl;
   cout << "# Column 1: ibin" << endl;
   cout << "# Column 2: Broad-histogram estimate of lng(E) Tinf" << endl;
   cout << "# Column 3: iterative eigenvector estimate" << endl;
   cout << "# Column 4: dgeev calculation of eigenvector" << endl;
   cout << "# Column 5: deltaE*beta" << endl;
   cout << "# Column 6: H total counts in bin" << endl;
   cout << setprecision(16);
   for(int irow=0; irow<nbin; irow++)
   {
      cout << irow << " " << lng[irow] << " " << ev[irow] << " " << ev_dgeev[irow] << " " << beta[irow] << " " << Hsum[irow] << endl;
   }
      
}
