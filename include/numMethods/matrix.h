/**
 *  @file      matrix.h
 *  @brief     Vector and Matrix classes.
 *  @author    Mikica Kocic, Press et al
 *  @copyright TODO
 */

#ifndef _MATRIX_H_INCLUDED
#define _MATRIX_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g11 Numerical Methods                                                    */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
#ifdef _DEBUG
    #define assertBounds(v,e)  \
        do if(!(v)) { std::cerr << e << std::endl; throw(e); } while(0)
#else
    #define assertBounds(v,e)
#endif

/////////////////////////////////////////////////////////////////////////////////////////
/** A vector (or an array) of elements.
 */
template <class T>
class Vector
{
    Int nn; //!< The size of the array. The upper index is `nn-1`.
    T* v;   //!< The elements of the array.

public:
    typedef T value_type; //!< Makes `T` available externally

    /** Constructs the empty vector.
     */
    Vector ()
        : nn( 0 )
        , v( NULL )
    {}

    /** Destructs the vector, freeing allocated resources.
     */
    ~Vector ()
    {
        if( v != NULL ) {
            delete[] v;
        }
    }

    /** Constructs the zero-based vector.
     */
    explicit Vector( Int n )
        : nn( n )
        , v( n > 0 ? new T[n] : NULL )
    {}

    /** Constructs the vector initialized to a constant value.
     */
    Vector( Int n, const T& c )
        : nn( n )
        , v( n > 0 ? new T[n] : NULL )
    {
        for( Int i = 0; i < n; ++i ) {
            v[i] = c;
        }
    }

    /** Constructs the vector from an array.
     */
    Vector( Int n, const T* ap )
        : nn( n )
        , v( n > 0 ? new T[n] : NULL )
    {
        for( Int i = 0; i < n; ++i ) {
            v[i] = *ap++;
        }
    }

    Vector( const Vector<T>& rhs )
        : nn( rhs.nn )
        , v( nn > 0 ? new T[nn] : NULL )
    {
        for( Int i = 0; i < nn; ++i ) {
            v[i] = rhs[i];
        }
    }

    /** Assigns an existing vector. Postcondition: normal assignment via copying has
     *  been performed;  if vector and rhs were different sizes, vector has been resized
     *  to match the size of rhs.
     */
    Vector<T>& operator=( const Vector<T>& rhs )
    {
        if( this != &rhs )
        {
            if( nn != rhs.nn ) {
                if( v != NULL ) {
                    delete[] v;
                }
                nn = rhs.nn;
                v= nn > 0 ? new T[nn] : NULL;
            }
            for( Int i = 0; i < nn; ++i ) {
                v[i] = rhs[i];
            }
        }
        return *this;
    }

    T& operator[]( const Int i )
    {
        assertBounds( i >= 0 && i < nn, "Vector subscript out of bounds" );
        return v[i];
    }

    const T& operator[](const Int i) const
    {
        assertBounds( i >= 0 && i < nn, "Vector subscript out of bounds" );
        return v[i];
    }

    /** Returns the size.
     */
    Int size () const {
        return nn;
    }

    /** Resize the array.
     */
    void resize( Int newn )
    {
        if( newn != nn )
        {
            if( v != NULL ) {
                delete[] v;
            }
            nn = newn;
            v = nn > 0 ? new T[nn] : NULL;
        }
    }

    /** Resize and assign a constant value.
     */
    void assign( Int newn, const T& c )
    {
        if( newn != nn ) {
            if( v != NULL ) {
                delete[] v;
            }
            nn = newn;
            v = nn > 0 ? new T[nn] : NULL;
        }
        for( Int i = 0; i < nn; ++i ) {
            v[i] = c;
        }
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
/** A matrix (a two-dimensional array) of elements.
 */
template <class T>
class Matrix
{
    Int nn;  //!< The number of rows.
    Int mm;  //!< The number of columns.
    T** v;   //!< The elements of the two-dimensional array.

public:
    typedef T value_type; //!< Makes `T` available externally.

    Matrix ()
        : nn(0), mm(0)
        , v(NULL)
    {}

    ~Matrix ()
    {
        if( v != NULL ) {
            delete[] v[0];
            delete[] v;
        }
    }

    Matrix( Int n, Int m )
        : nn(n), mm(m)
        , v( n > 0 ? new T*[n] : NULL )
    {
        Int nel = m * n;
        if( v ) {
            v[0] = nel > 0 ? new T[nel] : NULL;
        }
        for( Int i = 1; i < n; ++i ) {
            v[i] = v[i-1] + m;
        }
    }

    Matrix(Int n, Int m, const T& c )
        : nn(n), mm(m)
        , v( n > 0 ? new T*[n] : NULL )
    {
        Int nel = m * n;
        if( v ) {
            v[0] = nel > 0 ? new T[nel] : NULL;
        }
        for( Int i = 1; i < n; ++i ) {
            v[i] = v[i-1] + m;
        }
        for( Int i = 0; i < n; ++i ) {
            for( Int j = 0; j < m; ++j ) {
                v[i][j] = c;
            }
        }
    }

    Matrix(Int n, Int m, const T* ap )
        : nn(n), mm(m)
        , v( n > 0 ? new T*[n] : NULL )
    {
        Int nel = m * n;
        if( v ) {
            v[0] = nel > 0 ? new T[nel] : NULL;
        }
        for( Int i = 1; i < n; ++i ) {
            v[i] = v[i-1] + m;
        }
        for( Int i = 0; i < n; ++i ) {
            for( Int j = 0; j < m; ++j ) {
                v[i][j] = *ap++;
            }
        }
    }

    Matrix( const Matrix& rhs )
        : nn(rhs.nn), mm(rhs.mm)
        , v( nn > 0 ? new T*[nn] : NULL )
    {
        Int nel = mm * nn;
        if( v ) {
            v[0] = nel > 0 ? new T[nel] : NULL;
        }
        for( Int i = 1; i< nn; ++i ) {
            v[i] = v[i-1] + mm;
        }
        for( Int i = 0; i< nn; ++i ) {
            for( Int j = 0; j<mm; ++j ) {
                v[i][j] = rhs[i][j];
            }
        }
    }

    /** The matrix assignment. Postcondition: normal assignment via copying has been
     *  performed; If matrix and rhs were different sizes, matrix has been resized
     * to match the size of rhs
     */
    Matrix<T>& operator=( const Matrix<T>& rhs )
    {
        if( this != &rhs ) {
            if (nn != rhs.nn || mm != rhs.mm) {
                if (v != NULL) {
                    delete[] v[0];
                    delete[] v;
                }
                nn = rhs.nn;
                mm = rhs.mm;
                v = nn > 0 ? new T*[nn] : NULL;
                Int nel = mm * nn;
                if( v ) {
                    v[0] = nel > 0 ? new T[nel] : NULL;
                }
                for( Int i = 1; i < nn; ++i )
                    v[i] = v[i-1] + mm;
            }
            for( Int i = 0; i< nn; ++i )
                for( Int j = 0; j<mm; ++j )
                    v[i][j] = rhs[i][j];
        }
        return *this;
    }

    /** Returns the pointer to a row.
     */
    T* operator[]( const Int i )
    {
        assertBounds( i >= 0 && i < nn, "Matrix subscript out of bounds" );
        return v[i];
    }

    const T* operator[]( const Int i ) const
    {
        assertBounds( i >= 0 && i < nn, "Matrix subscript out of bounds" );
        return v[i];
    }

    Int nrows () const {
        return nn;
    }

    Int ncols () const {
        return mm;
    }

    void resize( Int newn, Int newm )
    {
        if( newn != nn || newm != mm ) {
            if( v != NULL ) {
                delete[] v[0];
                delete[] v;
            }
            nn = newn;
            mm = newm;
            v = nn > 0 ? new T*[nn] : NULL;
            Int nel = mm * nn;
            if( v ) {
                v[0] = nel > 0 ? new T[nel] : NULL;
            }
            for( Int i = 1; i < nn; ++i ) {
                v[i] = v[i-1] + mm;
            }
        }
    }


    void assign( Int newn, Int newm, const T& c )
    {
        if( newn != nn || newm != mm ) {
            if( v != NULL ) {
                delete[] v[0];
                delete[] v;
            }
            nn = newn;
            mm = newm;
            v = nn > 0 ? new T*[nn] : NULL;
            Int nel = mm * nn;
            if( v ) v[0] = nel > 0 ? new T[nel] : NULL;
            for( Int i = 1; i < nn; ++i )
                v[i] = v[i-1] + mm;
        }
        for( Int i = 0; i < nn; ++i )
            for( Int j = 0; j < mm; ++j )
                v[i][j] = c;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////
// Datatypes needed by banmul and Bandec (compatible with NR3).

typedef const Vector<Real> VecReal_I;
typedef Vector<Real> VecReal, VecReal_O, VecReal_IO;

typedef Vector<Int> VecInt;

typedef const Matrix<Real> MatReal_I;
typedef Matrix<Real> MatReal, MatReal_O, MatReal_IO;
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _MATRIX_H_INCLUDED
