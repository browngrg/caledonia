#ifndef WL_WALKER_RAYLEIGH_HPP
#define WL_WALKER_RAYLEIGH_HPP

#include<cmath>

// DGEBAL - balance a general real matrix A
#ifdef USE_BLAS
extern "C" void dgeev_(char* jobvl, char* jobvr, int* N, double* a, int* lda, double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int *lwork, int* info);
extern "C" void dgemm_(char* transa, char* transb ,int* M, int* N, int* K, const double* alpha, const double* A, int* lda, const double* B, int* ldb, const double* beta, double* C, int* ldc);
extern "C" void dgemv_(char* trans, int* M, int* N, const double* alpha, const double* A, int* lda, const double* X, int* incx, const double* beta, double* y, int* incy);
extern "C" void dgetrf_(int* M, int* N, double* A, int* lda, int* ipivot, int* info);
extern "C" void dgetri_(int* N, double* A, int* lda, int* ipivot, double* work, int* lwork, int* info);
#endif

// Returns the the Rayleigh quotient of x with respect to A, r(x) = xT*A*x/(xT*x)
double Rayleigh(const std::vector<double>& A, const std::vector<double>& x)
{
#  ifndef USE_BLAS
   std::cout << __FILE__ << ":" << __LINE__ << " member function Rayleigh(A,x) requires -DUSE_BLAS" << std::endl;
   return 1;
#  endif
   // This implementation assumes xT*x = 1
   int N = x.size();
   char trans = 'N';
   int incr = 1;
   double alpha = 1;
   double beta = 0;
   std::vector<double> Ax(N,0);
   dgemv_(&trans,&N,&N,&alpha,&(A[0]),&N,&(x[0]),&incr,&beta,&(Ax[0]),&incr);
   double xTAx = 0;
   double xTx = 0;
   for(int i=0; i<N; i++) 
   {
      xTAx += x[i]*Ax[i];
      xTx  += x[i]*x[i];
   }
   return xTAx/xTx;
}

double RayleighQuotientIteration(const std::vector<double>& mat, std::vector<double>& evec, int niter=20)
{
#  ifndef USE_BLAS
   std::cout << __FILE__ << ":" << __LINE__ << " member function RayleighQuotientIteration requires -DUSE_BLAS" << std::endl;
   return 0;
#  endif
   // Converge to the eigenvector using the Rayleigh quotient method
   // Returns the associated eigenvalue
   // This is known to have problems for non-symmetric matrice
   // Maysum Panju, "Iterativec Methods for Computing Eigenvalues and Eigenvectors,"
   // The Waterloo Mathematics Review 1, 9-18 (2011).
   int N = evec.size();
   double vmag = 0; for(int i=0; i<N; i++) vmag += evec[i]*evec[i];      // Normalize the eigenvector
   vmag = sqrt(vmag);
   for(int i=0; i<N; i++) evec[i] /= vmag;
   double lambda = Rayleigh(mat,evec);
   //cout << "# converge init lambda=" << lambda << std::endl;
   for(int k=0; k<niter; k++)
   {
      // Solve (A-lambda*I)*w = v
      std::vector<double> matinv = mat;
      for(int i=0; i<N; i++) matinv[i+i*N] -= lambda;                                    // mat = A - lambda*I
      std::vector<int> ipivot(N);
      int info;
      dgetrf_(&N,&N,&(matinv[0]),&N,&(ipivot[0]),&info);                                 // LU factorization of matrix A=P*L*U
      std::vector<double> work(10*N);
      int lwork = work.size();
      dgetri_(&N,&(matinv[0]),&N,&(ipivot[0]),&(work[0]),&lwork,&info);                  // Computes inverse matrix of LU factorized matrix
      char trans = 'N';
      int incr = 1;
      double alpha = 1;
      double beta = 0;
      std::vector<double> w(N);
      dgemv_(&trans,&N,&N,&alpha,&(matinv[0]),&N,&(evec[0]),&incr,&beta,&(w[0]),&incr);  // w = mat^{-1}*v
      double wmag = 0; for(int i=0; i<N; i++) wmag += w[i]*w[i];
      wmag = sqrt(wmag);
      for(int i=0; i<N; i++) evec[i] = w[i]/wmag;
      lambda = Rayleigh(mat,evec);
   }
   return lambda;
}

#endif
