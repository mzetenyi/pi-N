#ifndef VECTORS_HPP
#define VECTORS_HPP

#include <iostream>
#include <cmath>

using namespace std;

namespace Vectors {

  bool isIndex3(int);
  bool isIndex4(int);

  /**
     Forward declarations.
  */
  class ThreeVector;
  class ThreeTensor;
  class FourVector;
  class FourTensor;

  /**
     Represent a threevector.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, etc.) 
  */
  class ThreeVector {
  public:
    /// Create a nullvector.
    ThreeVector();
    /// Copy constructor.
    ThreeVector(const ThreeVector&);
    /// Create a vector with the given components.
    ThreeVector(double x1, double x2, double x3);
    /// Destructor. Does nothing.
    ~ThreeVector();
    /// Assignment operator.
    ThreeVector& operator = (const ThreeVector&);
    /// Return vector component.
    double& operator [] (int);
    const double& operator [] (int) const;
    double& operator () (int);
    const double& operator () (int) const;
    /// Positive sign.
    const ThreeVector& operator + () const;
    /// Negative sign.
    const ThreeVector operator - () const;
    /// Vector addition.
    const ThreeVector operator + (const ThreeVector&) const;
    /// Vector subtraction.
    const ThreeVector operator - (const ThreeVector&) const;
    /// Increment operator.
    ThreeVector& operator += (const ThreeVector&);
    /// Decrement operator.
    ThreeVector& operator -= (const ThreeVector&);
    /// Inner product.
    const double operator * (const ThreeVector&) const;
    /// Multiplication with scalar (double).
    const ThreeVector operator * (const double&) const;
    /// Division by scalar.
    const ThreeVector operator / (const double&) const;
    ThreeVector& operator *= (const double&);
    ThreeVector& operator /= (const double&);
    /// Friend operator for scalar*vector.
    friend const ThreeVector operator * (const double&, const ThreeVector&);
    /// Vectorial product
    friend ThreeVector crossprod(const ThreeVector&, const ThreeVector&);
    /// Equality check.
    bool operator == (const ThreeVector&) const;
    /// Unequality check.
    bool operator != (const ThreeVector&) const;
    /// Print the vector to the given stream in the form '(x,y)'.
    friend ostream& operator << (ostream&, const ThreeVector&);
    /// Read the vector from the given stream. It must be in the form '(x,y)'.
    friend istream& operator >> (istream&, ThreeVector&);
    /// Absolute value.
    double abs() const;
    /// Azimuth angle in spherical coordinates.
    double theta() const;
    /// Polar angle in spherical coordinated.
    double phi() const;
    /// Square of the vector.
    double square() const;
    /// Decide if any component is NaN.
    bool isNaN() const;
    /// Null vector.
    static const ThreeVector nullVector;
    /// Basis vectors.
    static const ThreeVector e[];
  private:
    /// Variables to store the vector components. (x[0] is not used)
    double x[4];
  };

  /**
     Represent a threetensor.
     Has the usual functionality (addition, multiplication with scalar, 
     tensor product, etc.) 
  */
  class ThreeTensor {
  public:
    /// Create a nulltensor.
    ThreeTensor();
    /// Copy constructor.
    ThreeTensor(const ThreeTensor&);
    /// Create a tensor from components.
    ThreeTensor(double,double,double,double,double,double,double,double,double);
    /// Create a diagonal tensor with the given components.
    ThreeTensor(double,double,double);
    /// Destructor. Does nothing.
    ~ThreeTensor();
    /// Assignment operator.
    ThreeTensor& operator = (const ThreeTensor&);
    /// Return tensor component.
    double& operator () (int,int);
    const double& operator () (int,int) const;
    /// Positive sign.
    const ThreeTensor& operator + () const;
    /// Negative sign.
    const ThreeTensor operator - () const;
    /// Tensor addition.
    const ThreeTensor operator + (const ThreeTensor&) const;
    /// Tensor subtraction.
    const ThreeTensor operator - (const ThreeTensor&) const;
    /// Increment operator.
    ThreeTensor& operator += (const ThreeTensor&);
    /// Decrement operator.
    ThreeTensor& operator -= (const ThreeTensor&);
    /// Tensor product.
    const ThreeTensor operator * (const ThreeTensor&) const;
    /// Tensor*vector product.
    const ThreeVector operator * (const ThreeVector&) const;
    /// Vector*tensor product.
    friend const ThreeVector operator * (const ThreeVector&, const ThreeTensor&);
    /// Multiplication with scalar (double).
    const ThreeTensor operator * (const double&) const;
    /// Division by scalar.
    const ThreeTensor operator / (const double&) const;
    ThreeTensor& operator *= (const double&);
    ThreeTensor& operator /= (const double&);
    /// Friend operator for scalar*tensor.
    friend const ThreeTensor operator * (const double&, const ThreeTensor&);
    /// Equality check.
    bool operator == (const ThreeTensor&) const;
    /// Unequality check.
    bool operator != (const ThreeTensor&) const;
    /// Print the tensor to the given stream in the form '((x11,x12,x13),(...),(...))'.
    friend ostream& operator << (ostream&, const ThreeTensor&);
    /// Print to a stream in the form of a table.
    void output(ostream&, int fieldWidth=10) const;
    /// Trace.
    double trace() const;
    /// create a diagonal tensor with the components of the given vector
    static ThreeTensor diagonalTensor(const ThreeVector&);
    /// Null tensor.
    static const ThreeTensor nullTensor;
    /// Unit tensor.
    static const ThreeTensor unitTensor;
  private:
    /// Variables to store the tensor components. (x[0] is not used)
    double x[4][4];
  };

  /// Diadic product.
  const ThreeTensor diad(const ThreeVector&, const ThreeVector&);

  /**
     Represent a fourvector.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, etc.) 

     \todo spatial() could return a reference (FourVector and
     ThreeVector have exactly equal representations.
  */
  class FourVector {
  public:
    /// Create a nullvector.
    FourVector();
    /// Copy constructor.
    FourVector(const FourVector&);
    /// Create a vector with the given components.
    FourVector(double x0, double x1, double x2, double x3);
    /// Create a vector from a three scalar and a threevector.
    FourVector(double x0, const ThreeVector& xx);
    /// Destructor. Does nothing.
    ~FourVector();
    /// Assignment operator.
    FourVector& operator = (const FourVector&);
    /// Return contravariant vector component.
    double& operator [] (int);
    /// Return contravariant vector component.
    const double& operator [] (int) const;
    /// Return covariant vector component.
    const double operator () (int) const;
    /// Positive sign.
    const FourVector& operator + () const;
    /// Negative sign.
    const FourVector operator - () const;
    /// Vector addition.
    const FourVector operator + (const FourVector&) const;
    /// Vector subtraction.
    const FourVector operator - (const FourVector&) const;
    /// Increment operator.
    FourVector& operator += (const FourVector&);
    /// Decrement operator.
    FourVector& operator -= (const FourVector&);
    /// Inner product.
    const double operator * (const FourVector&) const;
    /// Vector*tensor product.
    friend const FourVector operator * (const FourVector&, const FourTensor&);
    /// Diadic product.
    friend const FourTensor diad(const FourVector&, const FourVector&);
    /// Multiplication with scalar (double).
    const FourVector operator * (const double&) const;
    /// Division by scalar.
    const FourVector operator / (const double&) const;
    FourVector& operator *= (const double&);
    FourVector& operator /= (const double&);
    /// Friend operator for scalar*vector.
    friend const FourVector operator * (const double&, const FourVector&);
    /// Equality check.
    bool operator == (const FourVector&) const;
    /// Unequality check.
    bool operator != (const FourVector&) const;
    /// Print the vector to the given stream.
    friend ostream& operator << (ostream&, const FourVector&);
    /// Read in components of the vector from the given stream.
    friend istream& operator >> (istream&, FourVector&);
    /// Square of the vector.
    double square() const;
    /// Return the spacelike part as a ThreeVector.
    ThreeVector spacial() const;
    /// Reverse the spacial part of the FourVector. (Same as reflect())
    const FourVector reverse() const;
    /// Perform a space reflection on the FourVector.
    const FourVector reflect() const;
    /// Perform a space reflection on the FourVector.
    friend FourVector reflect(const FourVector&);
    /// Decide if the FourVector is timelike.
    bool timelike() const;
    /// Decide if the FourVector is timelike and future pointing.
    bool future() const;
    /// Decide if the FourVector is timelike and past pointing.
    bool past() const;
    /// Decide if the FourVector is spacelike.
    bool spacelike() const;
    /// Decide if any component is NaN.
    bool isNaN() const;
    /// Null vector.
    static const FourVector nullVector;
    /// Basis vectors.
    static const FourVector e[];
    /// Declare FourTensor as friend.
    friend class FourTensor;
  private:
    /// Variables to store the vector components.
    double x[4];
  };

  /**
     Represent a fourtensor.
     Has the usual functionality (addition, multiplication with scalar, 
     inner product, etc.) 
     The components are stored in the array x[0..4][0..4].

     Convention:

     In the storage:
     FIRST index is UPPER, SECOND index is LOWER.

     In the constructors:
     BOTH indices are UPPER.
  */
  class FourTensor {
  public:
    /// Create a nulltensor.
    FourTensor();
    /// Copy constructor.
    FourTensor(const FourTensor&);
    /// Create a tensor from components. (Indices are upper.)
    FourTensor(double,double,double,double,double,double,double,double,
               double,double,double,double,double,double,double,double);
    /// Create a diagonal tensor with the given components. (Indices are upper.)
    FourTensor(double x0, double x1, double x2, double x3);
    /// Create a Fourtensor from blocks. (Indices are upper.)
    FourTensor(double x00, const ThreeVector& x0j,  const ThreeVector& xi0, const ThreeTensor& xij);
    /// Destructor. Does nothing.
    ~FourTensor();
    /// Assignment operator.
    FourTensor& operator = (const FourTensor&);
    /// Return lower index (covariant) tensor component.
    double operator () (int,int) const;
    /// Return contravariant and covariant components.
    const FourVector operator [] (int) const;
    const FourVector operator () (int) const;
    /// Positive sign.
    const FourTensor& operator + () const;
    /// Negative sign.
    const FourTensor operator - () const;
    /// Tensor addition.
    const FourTensor operator + (const FourTensor&) const;
    /// Tensor subtraction.
    const FourTensor operator - (const FourTensor&) const;
    /// Increment operator.
    FourTensor& operator += (const FourTensor&);
    /// Decrement operator.
    FourTensor& operator -= (const FourTensor&);
    /// Tensor product.
    const FourTensor operator * (const FourTensor&) const;
    /// Tensor*vector product.
    const FourVector operator * (const FourVector&) const;
    /// Vector*tensor product.
    friend const FourVector operator * (const FourVector&, const FourTensor&);
    /// Diadic product.
    friend const FourTensor diad(const FourVector&, const FourVector&);
    /// Multiplication with scalar (double).
    const FourTensor operator * (const double&) const;
    /// Division by scalar.
    const FourTensor operator / (const double&) const;
    FourTensor& operator *= (const double&);
    FourTensor& operator /= (const double&);
    /// Friend operator for scalar*tensor.
    friend const FourTensor operator * (const double&, const FourTensor&);
    /// Equality check.
    bool operator == (const FourTensor&) const;
    /// Unequality check.
    bool operator != (const FourTensor&) const;
    /// Print the tensor to the given stream in the form '(x00,x01,x02,x03),(..),(..),(..))'.
    friend ostream& operator << (ostream&, const FourTensor&);
    /// Print to a stream in the form of a table.
    void output(ostream&, int fieldWidth=10) const;
    /// Trace.
    double trace() const;
    /// Null tensor.
    static const FourTensor nullTensor;
    /// Unit tensor.
    static const FourTensor unitTensor;
    /// Declare FourVector as friend.
    friend class FourVector;
  private:
    /// Variables to store the tensor components.
    double x[4][4];
  };

  /**
     Represents a Lorentz transformation that transforms fourvectors into
     the rest frame of the fourvector a.
  */
  class LorentzTransformation {
  public:
    /**
       Constructor that creates a LorentzTransformation transforming to
       the rest frame of a.
    */
    LorentzTransformation(const FourVector& a);
    /// Transform the fourvector x to the rest frame of a.
    FourVector operator()(const FourVector& x) const;
    /// Transform the fourtensor T to the rest frame of a.
    FourTensor operator()(const FourTensor& T) const;
    /**
       Transform the quantity x to the rest frame of a. This operator
       allows the use of operator() to perform the Lorentz
       transfromation on an object of any type T, that has a member
       function transform(const LorentzTransformation&) - this member
       implements the transformation.
    */
    template <typename T>
    T operator()(const T& x) const {
      return x.transform(*this);
    }
    /// The inverse Lorentz transformation.
    LorentzTransformation inverse() const;

    /**
       Print components of the Lorentz tranformation matrix.
    */
    friend ostream& operator<<(ostream&, const LorentzTransformation&);
    /**
       Generate a LorentzTransformation that transforms to the rest
       frame of the given FourVector a.
    */
    static LorentzTransformation toRestFrame(const FourVector& a);
    /**
       Generate a LorentzTransformation that transforms a FourVector at
       rest to a FourVector parallel to the given FourVector a.
    */
    static LorentzTransformation fromRestFrame(const FourVector& a);
  private:
    const FourVector a;
    FourTensor L;
  };

  /**
     Signature of the metric.
  */
  double sign_(int);

  /**
     Metric tensor.
  */
  extern const FourTensor g_;

  /**
     Totally antisymmetric epsilon tensor. EpsilonTensor(0,1,2,3) = 1,
     i.e. indices are upper.
  */
  int EpsilonTensor(int mu, int nu, int ro, int si);
} // end of namespace Vectors

#endif
