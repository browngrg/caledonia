#ifndef PSIMAG_VEC_HPP_
#define PSIMAG_VEC_HPP_

#include<iostream>
#include<cmath>
#include<cstddef>

// This was originally Vec.h from the psimag toolset 
// @author Thomas C. Schulthess
  
// 
// Cartesian vector with fixed dimensions. 
// It is implemented as a template to allow different types of elements, and
// the lower dimensions are implemented through partial specialization
// for better efficiency.
// @author T. Schulthess and G. Brown, ORNL, October 2000
//
template < class T, size_t DIM > 
class Vec {
    
    template <class S, size_t D> friend
    std::istream& operator >> (std::istream& is, Vec<S,D>& v);

    template <class S, size_t D> friend
    std::ostream& operator << (std::ostream& os, const Vec<S,D> &v);
    
  public:
    
    // Default constructor, initializes all components to 0.
    Vec() { for(size_t i=0;i<DIM;i++) d[i]=0; }
    
    // Construct a vector with all components set to value of a.
    explicit Vec(const T& a) { for(size_t i=0;i<DIM;i++) d[i]=a; }
    
    // Construct a vector with components \f$v_i = a[ {\rm stride} * i]\f$.
    Vec(const T* const a, size_t stride=1) { for(size_t i=0;i<DIM;i++) d[i]=a[i*stride]; }
    
    // This copy constructor makes deep copies. 
    Vec(const Vec<T,DIM>& v) { for(size_t i=0;i<DIM;i++) d[i]= v.d[i]; }
    
    // Create a null default constructor for safety. 
    ~Vec() {}
 
    // Gives the size of the vector 
    static size_t size () { return DIM; }
    static const size_t static_size=DIM;

    // Return constant reference to i<SUP>th</SUP> component, note zero i=0,1, ,DIM-1. 
    const T& operator [] (size_t i) const { return d[i]; }
    
    // Return reference to i<SUP>th</SUP> component,
    T& operator [] (size_t i) { return d[i]; }
    
    // Assign the value of a to each component.
    Vec<T,DIM>& operator = (const Vec<T,DIM>& v) { for(size_t i=0;i<DIM;i++) d[i] = v.d[i]; return *this; } 

    // Assign the value of a to each component.
    Vec<T,DIM>& operator = (const T& a) { for(size_t i=0;i<DIM;i++) d[i] = a; return *this; }

    // Increment component-wise by vector v.
    Vec<T,DIM>& operator += (const Vec<T,DIM>& v) { for(int i=0; i<DIM; i++) d[i] += v.d[i]; return *this; }

    // Decrement component-wise by vector v.
    Vec<T,DIM>& operator -= (const Vec<T,DIM>& v) { for(int i=0; i<DIM; i++) d[i] -= v.d[i]; return *this; }

    // Scale vector by scalar a.
    Vec<T,DIM>& operator *= (const T& a) { for(int i=0; i<DIM; i++) d[i] *= a; return *this; }

    // Scale vector by scalar 1/a.
    Vec<T,DIM>& operator /= (const T& a) { for(int i=0; i<DIM; i++) d[i] /= a; return *this; }

  private:
    
    T d[DIM];
    
  };
  
  template <class T, size_t DIM>
  inline std::istream& operator >> (std::istream& is, Vec<T,DIM>& v)
  { 
    for(size_t i=0;i<DIM;i++) 
      is >> v.d[i]; 
    return is; 
  }  
  
  template <class T, size_t DIM>
  inline std::ostream& operator << (std::ostream& os, const Vec<T,DIM> &v)
  { 
    os << v.d[0]; 
    for(size_t i=1;i<DIM;i++) 
      os << "\t" << v.d[i]; 
    return os; 
  }
  
  // L1Norm \f$= \sum_{j=0}^{DIM-1} \|x_i\|\f$.
  template < class T, size_t DIM > 
  T L1Norm(const Vec <T,DIM> &a )
  { 
    T val=0;
    for(size_t i=0; i<DIM; i++)
      val += std::abs(a[i]);
    return val;
  }
  
  // Returns scalar product of a and b. 
  template < class T, size_t DIM >
  T operator * (const Vec<T,DIM>& a, 
		const Vec<T,DIM>& b)
  { 
    T val=0;
    for(size_t i=0; i<DIM; i++)
      val += a[i]*b[i];
    return val;
  }
  
  // Left-multiplication of vector by scalar a.
  template <class T, size_t DIM>
  Vec<T,DIM> operator * (const T& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val=a;
    for(size_t i=0; i<DIM; i++)
      val[i] *= a;
    return val;
  }

  // Right-multiplication of vector by scalar a.
  template < class T, size_t DIM >
  Vec<T,DIM> operator * (const Vec<T,DIM>& b, 
			  const T& a)
  { 
    Vec<T,DIM> val=b;
    for(size_t i=0; i<DIM; i++)
      val[i] *= a;
    return val;
  }

  // Division of a vector by a scalar a.
  template < class T, size_t DIM >
  Vec<T,DIM> operator / (const Vec<T,DIM>& b, 
			  const T& a)
  { 
    Vec<T,DIM> val=b;
    for(size_t i=0; i<DIM; i++)
      val[i] /= a;
    return val;
  }

  // Vector addition.
  template < class T, size_t DIM >
  Vec<T,DIM> operator + (const Vec<T,DIM>& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = a[i] + b[i];
    return val;
  }

  // Vector subtraction.
  template < class T, size_t DIM >
  Vec<T,DIM> operator - (const Vec<T,DIM>& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = a[i] - b[i];
    return val;
  }

  // Negation.
  template < class T, size_t DIM >
  Vec<T,DIM> operator - (const Vec<T,DIM>& a)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = -a[i];
    return val;
  }

  // Equivalence operator
  template<class T, size_t DIM>
  bool operator == (const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
     //
     // old code, semantics easy to see
     // bool equal = true;
     // for(size_t i=0; equal && i<DIM; i++) if(a[i]!=b[i]) equal=false;
     // return equal;
     //
     // new code, HKL thinks its faster. Depends on optimizer
     for(size_t i=0; i<DIM; i++) if(a[i]!=b[i]) return false;
     return true;
  }

  // Nonequivalence operator
  template<class T, size_t DIM>
  bool operator != (const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
    //bool nequal = false;
    //for(size_t i=0; !nequal && i<DIM; i++) if(a[i]!=b[i]) nequal=true;
    //return nequal;
    for(size_t i=0; i<DIM; i++) if(a[i]!=b[i]) return true; 
    return false;
  }
 
  // Ordering operator
  template<class T, size_t DIM>
  bool operator < ( const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
      // bool lessthan = false;
      // bool equal = true;
      // for(size_t i=0; equal && !lessthan && i<DIM; i++) {
      //    if(a[i]<b[i]) lessthan=true;
      //    else if(a[i]>b[i]) equal = false;
      // }
      // return lessthan;
      for(size_t i=0; i<DIM; i++)
      {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
      }
      return false;
  }

  template<class T, size_t DIM>
  bool operator > ( const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
      // bool lessthan = false;
      // bool equal = true;
      // for(size_t i=0; equal && !lessthan && i<DIM; i++) {
      //    if(a[i]<b[i]) lessthan=true;
      //    else if(a[i]>b[i]) equal = false;
      // }
      // return lessthan;
      for(size_t i=0; i<DIM; i++)
      {
        if(a[i]>b[i]) return true;
        if(a[i]<b[i]) return false;
      }
      return false;
  }


// --------------------------------------------------------------------------
  
//  Cartesian vector in 3D (specialized template Vec<T,DIM=3>).
//  All methods of Vec<T,DIM> apply to Vec<T,3>.
//  @author T. Schulthess and G. Brown, ORNL, October 2000
template < class T > 
class Vec< T, 3 > {
public:
    template < class S> friend std::istream& operator >> (std::istream& is, Vec<S,3>& v);  
    template < class S> friend std::ostream& operator << (std::ostream& os, const Vec<S,3> &v);
    Vec() { d[0]=d[1]=d[2]=0; }
    explicit Vec(const T& a) { d[0]=d[1]=d[2]=a; }
    Vec(const T& x,const T& y,const T& z) { d[0]=x; d[1]=y; d[2]=z; }
    Vec(const T * const a,size_t stride=1) { d[0]=a[0]; d[1]=a[stride]; d[2]=a[2*stride]; }
    Vec(const Vec<T,3>& v) { d[0]=v.d[0]; d[1]=v.d[1]; d[2]=v.d[2]; }
    ~Vec() {}
    static size_t size () { return 3; }
    static const size_t static_size=3;
    const T& operator [] (size_t i) const { return d[i]; }
    T& operator [] (size_t i) { return d[i]; }
    Vec<T,3>& operator = (const Vec<T,3>& v) { d[0]=v.d[0]; d[1]=v.d[1]; d[2]=v.d[2]; return *this; }
    Vec<T,3>& operator = (const T& a) { d[0]=d[1]=d[2]=a; return *this; }
    Vec<T,3>& operator += (const Vec<T,3>& v) { d[0]+=v.d[0]; d[1]+=v.d[1]; d[2]+=v.d[2]; return *this; }
    Vec<T,3>& operator -= (const Vec<T,3>& v) { d[0]-=v.d[0]; d[1]-=v.d[1]; d[2]-=v.d[2]; return *this; }
    Vec<T,3>& operator *= (const T& a) { d[0]*=a; d[1]*=a; d[2]*=a; return *this; }
  private:
    T d[3];
};

template < class T > 
std::istream& operator >> (std::istream& is, Vec<T,3>& v)
{ 
  for(size_t i=0;i<3;i++) 
    is >> v.d[i]; 
  return is; 
}  

template < class T > 
std::ostream& operator << (std::ostream& os, const Vec<T,3> &v)
{ 
  os << v.d[0]; 
  for(size_t i=1;i<3;i++) 
    os << "\t" << v.d[i]; 
  return os; 
}

template < class T > T L1Norm(const Vec<T,3>& a) { return std::abs(a[0])+std::abs(a[1])+std::abs(a[2]); } 
template < class T > T operator * (const Vec<T,3>& a, const Vec<T,3>& b) { return T(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]); }
template < class T > Vec<T,3> operator * (const T& a, const Vec<T,3>& b) { return Vec<T,3> (a * b[0], a * b[1], a * b[2]); }
template < class T > Vec<T,3> operator * (const Vec<T,3>& b, const T& a) { return Vec<T,3> (a * b[0], a * b[1], a * b[2]); }
template < class T > Vec<T,3> operator / (const Vec<T,3>& b, const T& a) { return Vec<T,3> (b[0]/a, b[1]/a, b[2]/a); }
template < class T > Vec<T,3> operator + (const Vec<T,3>& a, const Vec<T,3>& b) { return Vec<T,3> (a[0]+b[0],a[1]+b[1],a[2]+b[2]); }
template < class T > Vec<T,3> operator - (const Vec<T,3>& a, const Vec<T,3>& b) { return Vec<T,3> (a[0]-b[0],a[1]-b[1],a[2]-b[2]); }
template < class T > Vec<T,3> operator - (const Vec<T,3>& a) { return Vec<T,3> (-a[0],-a[1],-a[2]); }

// --------------------------------------------------------------------------

// Cartesian vector in 2D (specialized template Vec<T,DIM=3>).
// All methods of Vec<T,DIM> apply to Vec<T,2>.
// @author T. Schulthess and G. Brown, ORNL, October 2000

template < class T > 
class Vec< T, 2 > {
public:
    template <class S> friend std::istream& operator >> (std::istream& is, Vec<S,2>& v);
    template <class S> friend std::ostream& operator << (std::ostream& os, const Vec<S,2> &v);
    Vec() { d[0]=d[1]=0; }
    explicit Vec(const T& a) { d[0]=d[1]=a; }
    Vec(const T& x,const T& y) { d[0]=x; d[1]=y; }
    Vec(const T * const a,size_t stride=1) {  d[0]=a[0]; d[1]=a[stride]; }
    Vec(const Vec<T,2>& v) { d[0]=v.d[0]; d[1]=v.d[1]; }
    ~Vec() {}
    static size_t size () { return 2; }
    static const size_t static_size=2;
    const T& operator [] (size_t i) const { return d[i]; }
    T& operator [] (size_t i) { return d[i]; }
    Vec<T,2>& operator = (const Vec<T,2>& v) { d[0]=v.d[0]; d[1]=v.d[1]; return *this; }
    Vec<T,2>& operator = (const T& a) { d[0]=d[1]=a; return *this; }
    Vec<T,2>& operator += (const Vec<T,2>& v) { d[0]+=v.d[0]; d[1]+=v.d[1]; return *this; }
    Vec<T,2>& operator -= (const Vec<T,2>& v) { d[0]-=v.d[0]; d[1]-=v.d[1]; return *this; }
    Vec<T,2>& operator *= (const T& a) { d[0]*=a; d[1]*=a; return *this; }
  private:
    T d[2];
};

template < class T > 
std::istream& operator >> (std::istream& is, Vec<T,2>& v)
{ 
  for(size_t i=0;i<2;i++) 
    is >> v.d[i]; 
  return is; 
}  

template < class T > 
std::ostream& operator << (std::ostream& os, const Vec<T,2> &v)
{ 
  os << v.d[0]; 
  for(size_t i=1;i<2;i++) 
    os << "\t" << v.d[i]; 
  return os; 
}

template < class T > T L1Norm(const Vec<T,2>& a) { return std::abs(a[0])+std::abs(a[1]); }
template < class T > T operator * (const Vec<T,2>& a, const Vec<T,2>& b) { return T(a[0]*b[0]+a[1]*b[1]); }
template < class T > T operator % (const Vec<T,2>& a, const Vec<T,2>& b) { return T(a[0] * b[1] - a[1] * b[0]); }
template < class T > Vec<T,2> operator * (const T& a, const Vec<T,2>& b) { return Vec<T,2> (a * b[0], a * b[1]); }
template < class T > Vec<T,2> operator * (const Vec<T,2>& b, const T& a) { return Vec<T,2> (a * b[0], a * b[1]); }
template < class T > Vec<T,2> operator / (const Vec<T,2>& b, const T& a) { return Vec<T,2> (b[0]/a, b[1]/a); }
template < class T > Vec<T,2> operator + (const Vec<T,2>& a, const Vec<T,2>& b) { return Vec<T,2> (a[0]+b[0],a[1]+b[1]); }
template < class T > Vec<T,2> operator - (const Vec<T,2>& a, const Vec<T,2>& b) { return Vec<T,2> (a[0]-b[0],a[1]-b[1]); }
template < class T > Vec<T,2> operator - (const Vec<T,2>& a) { return Vec<T,2> (-a[0],-a[1]); }

// --------------------------------------------------------------------------

//   These methods are overloaded to provide proper behavior for scalars.
//   Having these methods defined allows many algorithms that work for 
//   Heisenberg models to also work for Ising models.

//  L2Norm for scalar returns the absolute value of the scalar.
//  This function is defined to allow generic algorithms that only need
//  to know the length of a vector work equally well with scalar
//  quantities. This is frequently needed for algorithms that are 
//  equally valid for Heisenberg and Ising models.
//  There is a good implementation question here: are these specializations
//  better replaced by a true templated function that returns the absolute
//  value, so that only the Vec forms are specializations?
//  @author G. Brown, ORNL, September 2004

template<class T> T L2Norm(T a) { return std::abs(a); }

template < class T, size_t DIM > 
T L2Norm(const Vec <T,DIM> &a)
{
  T val=0;
  for(size_t i=0; i<DIM; i++)
    val += a[i]*a[i];
  return sqrt(val);
}
  
template < class T > T L2Norm(const Vec<T,2>& a) { return sqrt(a[0]*a[0]+a[1]*a[1]); }
template < class T > T L2Norm(const Vec<T,3>& a) { return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); }


#endif   // PSIMAG_VEC_HPP_
