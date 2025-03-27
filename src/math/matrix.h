#ifndef MATRIXH
#define MATRIXH

#include <cstring>
#include <span>
#include "../math/vectypes.h"
#include "../math/simd.h"
#include "../math/mathinline.h"
#include <optional>

struct alignas(16) Matrix4x4 {
  Matrix4x4() {
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
    m[0][1] = m[0][2] = m[0][3] = m[1][0] =
    m[1][2] = m[1][3] = m[2][0] = m[2][1] = m[2][3] =
    m[3][0] = m[3][1] = m[3][2] = 0.f;
  }
  Matrix4x4(Float mat[4][4]);
#ifdef RAYSIMDVEC
  Matrix4x4(FVec4 mat[4]) {
    m[0] = mat[0];
    m[1] = mat[1];
    m[2] = mat[2];
    m[3] = mat[3];
  }
#endif
  Matrix4x4(Float t00, Float t01, Float t02, Float t03,
            Float t10, Float t11, Float t12, Float t13,
            Float t20, Float t21, Float t22, Float t23,
            Float t30, Float t31, Float t32, Float t33);
  bool operator==(const Matrix4x4 &m2) const {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (m[i][j] != m2.m[i][j]) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(const Matrix4x4 &m2) const {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (m[i][j] != m2.m[i][j]) {
          return true;
        }
      }
    }
    return false;
  }
  friend Float Det(const Matrix4x4 &m1) {
    return(
        m1.m[0][3] * m1.m[1][2] * m1.m[2][1] * m1.m[3][0] - m1.m[0][2] * m1.m[1][3] * m1.m[2][1] * m1.m[3][0] -
        m1.m[0][3] * m1.m[1][1] * m1.m[2][2] * m1.m[3][0] + m1.m[0][1] * m1.m[1][3] * m1.m[2][2] * m1.m[3][0] +
        m1.m[0][2] * m1.m[1][1] * m1.m[2][3] * m1.m[3][0] - m1.m[0][1] * m1.m[1][2] * m1.m[2][3] * m1.m[3][0] -
        m1.m[0][3] * m1.m[1][2] * m1.m[2][0] * m1.m[3][1] + m1.m[0][2] * m1.m[1][3] * m1.m[2][0] * m1.m[3][1] +
        m1.m[0][3] * m1.m[1][0] * m1.m[2][2] * m1.m[3][1] - m1.m[0][0] * m1.m[1][3] * m1.m[2][2] * m1.m[3][1] -
        m1.m[0][2] * m1.m[1][0] * m1.m[2][3] * m1.m[3][1] + m1.m[0][0] * m1.m[1][2] * m1.m[2][3] * m1.m[3][1] +
        m1.m[0][3] * m1.m[1][1] * m1.m[2][0] * m1.m[3][2] - m1.m[0][1] * m1.m[1][3] * m1.m[2][0] * m1.m[3][2] -
        m1.m[0][3] * m1.m[1][0] * m1.m[2][1] * m1.m[3][2] + m1.m[0][0] * m1.m[1][3] * m1.m[2][1] * m1.m[3][2] +
        m1.m[0][1] * m1.m[1][0] * m1.m[2][3] * m1.m[3][2] - m1.m[0][0] * m1.m[1][1] * m1.m[2][3] * m1.m[3][2] -
        m1.m[0][2] * m1.m[1][1] * m1.m[2][0] * m1.m[3][3] + m1.m[0][1] * m1.m[1][2] * m1.m[2][0] * m1.m[3][3] +
        m1.m[0][2] * m1.m[1][0] * m1.m[2][1] * m1.m[3][3] - m1.m[0][0] * m1.m[1][2] * m1.m[2][1] * m1.m[3][3] -
        m1.m[0][1] * m1.m[1][0] * m1.m[2][2] * m1.m[3][3] + m1.m[0][0] * m1.m[1][1] * m1.m[2][2] * m1.m[3][3]
    );
  }
  friend std::ostream& operator<<(std::ostream& o, Matrix4x4 const& m1)  {
    return (o << m1.m[0][0] << " " << m1.m[0][1] << " " << m1.m[0][2] << " " << m1.m[0][3] << "\n" <<
      m1.m[1][0] << " " << m1.m[1][1] << " " << m1.m[1][2] << " " << m1.m[1][3] << "\n" <<
      m1.m[2][0] << " " << m1.m[2][1] << " " << m1.m[2][2] << " " << m1.m[2][3] << "\n" << 
      m1.m[3][0] << " " << m1.m[3][1] << " " << m1.m[3][2] << " " << m1.m[3][3] << "\n");
  }
  friend Matrix4x4 Transpose(const Matrix4x4 &);

  static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
#ifdef RAYSIMDVEC
        FVec4 col{m2.m[0][j], m2.m[1][j], m2.m[2][j], m2.m[3][j]};
        r.m[i][j] = simd_dot(m1.m[i], col);
#else
        r.m[i][j] = m1.m[i][0] * m2.m[0][j] + 
                    m1.m[i][1] * m2.m[1][j] + 
                    m1.m[i][2] * m2.m[2][j] + 
                    m1.m[i][3] * m2.m[3][j];
#endif
      }
    }
    return r;
  }
  
  friend Matrix4x4 Inverse(const Matrix4x4 &);
  
#ifdef RAYSIMDVEC
  FVec4 m[4];
#else
  Float m[4][4];
#endif
};


template <int N>
inline void init(Float m[N][N], int i, int j) {}

template <int N, typename... Args>
inline void init(Float m[N][N], int i, int j, Float v, Args... args) {
    m[i][j] = v;
    if (++j == N) {
        ++i;
        j = 0;
    }
    init<N>(m, i, j, args...);
}

template <int N>
inline void initDiag(Float m[N][N], int i) {}

template <int N, typename... Args>
inline void initDiag(Float m[N][N], int i, Float v, Args... args) {
    m[i][i] = v;
    initDiag<N>(m, i + 1, args...);
}


// SquareMatrix Definition
template <int N>
class SquareMatrix {
  public:
    // SquareMatrix Public Methods
    
    static SquareMatrix Zero() {
        SquareMatrix m;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m.m[i][j] = 0;
        return m;
    }

    
    SquareMatrix() {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = (i == j) ? 1 : 0;
    }
    
    SquareMatrix(const Float mat[N][N]) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = mat[i][j];
    }
    
    SquareMatrix(std::span<const Float> t);
    template <typename... Args>
    SquareMatrix(Float v, Args... args) {
        static_assert(1 + sizeof...(Args) == N * N,
                      "Incorrect number of values provided to SquareMatrix constructor");
        init<N>(m, 0, 0, v, args...);
    }
    template <typename... Args>
    static SquareMatrix Diag(Float v, Args... args) {
        static_assert(1 + sizeof...(Args) == N,
                      "Incorrect number of values provided to SquareMatrix::Diag");
        SquareMatrix m;
        initDiag<N>(m.m, 0, v, args...);
        return m;
    }

    
    SquareMatrix operator+(const SquareMatrix &m) const {
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] += m.m[i][j];
        return r;
    }

    
    SquareMatrix operator*(Float s) const {
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] *= s;
        return r;
    }
    
    SquareMatrix operator/(Float s) const {
        //DCHECK_NE(s, 0);
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] /= s;
        return r;
    }

    
    bool operator==(const SquareMatrix<N> &m2) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (m[i][j] != m2.m[i][j])
                    return false;
        return true;
    }

    
    bool operator!=(const SquareMatrix<N> &m2) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (m[i][j] != m2.m[i][j])
                    return true;
        return false;
    }

    
    bool operator<(const SquareMatrix<N> &m2) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                if (m[i][j] < m2.m[i][j])
                    return true;
                if (m[i][j] > m2.m[i][j])
                    return false;
            }
        return false;
    }

    
    bool IsIdentity() const;

    std::string ToString() const;

    
    std::span<const Float> operator[](int i) const { return m[i]; }
    
    std::span<Float> operator[](int i) { return std::span<Float>(m[i]); }

  private:
    Float m[N][N];
};

// SquareMatrix Inline Methods
template <int N>
inline bool SquareMatrix<N>::IsIdentity() const {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                if (m[i][j] != 1)
                    return false;
            } else if (m[i][j] != 0)
                return false;
        }
    return true;
}

// SquareMatrix Inline Functions
template <int N>
inline SquareMatrix<N> operator*(Float s, const SquareMatrix<N> &m) {
    return m * s;
}

template <typename Tresult, int N, typename T>
inline Tresult Mul(const SquareMatrix<N> &m, const T &v) {
    Tresult result;
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j)
            result[i] += m[i][j] * v[j];
    }
    return result;
}

template <int N>
Float Determinant(const SquareMatrix<N> &m);

template <>
inline Float Determinant(const SquareMatrix<3> &m) {
    Float minor12 = DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
    Float minor02 = DifferenceOfProducts(m[1][0], m[2][2], m[1][2], m[2][0]);
    Float minor01 = DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
    return std::fma(m[0][2], 
                  minor01,
                  DifferenceOfProducts(m[0][0], minor12, m[0][1], minor02));
}

template <int N>
inline SquareMatrix<N> Transpose(const SquareMatrix<N> &m);
template <int N>
std::optional<SquareMatrix<N>> Inverse(const SquareMatrix<N> &);

template <int N>
SquareMatrix<N> InvertOrExit(const SquareMatrix<N> &m) {
    std::optional<SquareMatrix<N>> inv = Inverse(m);
    // CHECK(inv.has_value());
    return *inv;
}

template <int N>
inline SquareMatrix<N> Transpose(const SquareMatrix<N> &m) {
    SquareMatrix<N> r;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i][j] = m[j][i];
    return r;
}

template <>
inline std::optional<SquareMatrix<3>> Inverse(const SquareMatrix<3> &m) {
    Float det = Determinant(m);
    if (det == 0)
        return {};
    Float invDet = 1 / det;

    SquareMatrix<3> r;

    r[0][0] = invDet * DifferenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
    r[1][0] = invDet * DifferenceOfProducts(m[1][2], m[2][0], m[1][0], m[2][2]);
    r[2][0] = invDet * DifferenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
    r[0][1] = invDet * DifferenceOfProducts(m[0][2], m[2][1], m[0][1], m[2][2]);
    r[1][1] = invDet * DifferenceOfProducts(m[0][0], m[2][2], m[0][2], m[2][0]);
    r[2][1] = invDet * DifferenceOfProducts(m[0][1], m[2][0], m[0][0], m[2][1]);
    r[0][2] = invDet * DifferenceOfProducts(m[0][1], m[1][2], m[0][2], m[1][1]);
    r[1][2] = invDet * DifferenceOfProducts(m[0][2], m[1][0], m[0][0], m[1][2]);
    r[2][2] = invDet * DifferenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);

    return r;
}

template <int N, typename T>
inline T operator*(const SquareMatrix<N> &m, const T &v) {
    return Mul<T>(m, v);
}

template <>
inline SquareMatrix<4> operator*(const SquareMatrix<4> &m1,
                                 const SquareMatrix<4> &m2) {
    SquareMatrix<4> r;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            r[i][j] = InnerProduct(m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2],
                                   m2[2][j], m1[i][3], m2[3][j]);
    return r;
}

template <>
inline SquareMatrix<3> operator*(const SquareMatrix<3> &m1,
                                              const SquareMatrix<3> &m2) {
    SquareMatrix<3> r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[i][j] =
                InnerProduct(m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2], m2[2][j]);
    return r;
}

template <int N>
inline SquareMatrix<N> operator*(const SquareMatrix<N> &m1,
                                              const SquareMatrix<N> &m2) {
    SquareMatrix<N> r;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            r[i][j] = 0;
            for (int k = 0; k < N; ++k)
                r[i][j] = FMA(m1[i][k], m2[k][j], r[i][j]);
        }
    return r;
}

template <int N>
inline SquareMatrix<N>::SquareMatrix(std::span<const Float> t) {
    //CHECK_EQ(N * N, t.size());
    for (int i = 0; i < N * N; ++i)
        m[i / N][i % N] = t[i];
}

template <int N>
SquareMatrix<N> operator*(const SquareMatrix<N> &m1,
                                       const SquareMatrix<N> &m2);

template <>
inline Float Determinant(const SquareMatrix<1> &m) {
    return m[0][0];
}

template <>
inline Float Determinant(const SquareMatrix<2> &m) {
    return DifferenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);
}

template <>
inline Float Determinant(const SquareMatrix<4> &m) {
    Float s0 = DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
    Float s1 = DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
    Float s2 = DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

    Float s3 = DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
    Float s4 = DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
    Float s5 = DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

    Float c0 = DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
    Float c1 = DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
    Float c2 = DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

    Float c3 = DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
    Float c4 = DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
    Float c5 = DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

    return (DifferenceOfProducts(s0, c5, s1, c4) + DifferenceOfProducts(s2, c3, -s3, c2) +
            DifferenceOfProducts(s5, c0, s4, c1));
}

template <int N>
inline Float Determinant(const SquareMatrix<N> &m) {
    SquareMatrix<N - 1> sub;
    Float det = 0;
    // Inefficient, but we don't currently use N>4 anyway..
    for (int i = 0; i < N; ++i) {
        // Sub-matrix without row 0 and column i
        for (int j = 0; j < N - 1; ++j)
            for (int k = 0; k < N - 1; ++k)
                sub[j][k] = m[j + 1][k < i ? k : k + 1];

        Float sign = (i & 1) ? -1 : 1;
        det += sign * m[0][i] * Determinant(sub);
    }
    return det;
}

template <>
inline std::optional<SquareMatrix<4>> Inverse(const SquareMatrix<4> &m) {
    // Via: https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
    // (c) Google, Apache license.

    // For 4x4 do not compute the adjugate as the transpose of the cofactor
    // matrix, because this results in extra work. Several calculations can be
    // shared across the sub-determinants.
    //
    // This approach is explained in David Eberly's Geometric Tools book,
    // excerpted here:
    //   http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
    Float s0 = DifferenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
    Float s1 = DifferenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
    Float s2 = DifferenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

    Float s3 = DifferenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
    Float s4 = DifferenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
    Float s5 = DifferenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

    Float c0 = DifferenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
    Float c1 = DifferenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
    Float c2 = DifferenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

    Float c3 = DifferenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
    Float c4 = DifferenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
    Float c5 = DifferenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

    Float determinant = InnerProduct(s0, c5, -s1, c4, s2, c3, s3, c2, s5, c0, -s4, c1);
    if (determinant == 0)
        return {};
    Float s = 1 / determinant;

    Float inv[4][4] = {{s * InnerProduct(m[1][1], c5, m[1][3], c3, -m[1][2], c4),
                        s * InnerProduct(-m[0][1], c5, m[0][2], c4, -m[0][3], c3),
                        s * InnerProduct(m[3][1], s5, m[3][3], s3, -m[3][2], s4),
                        s * InnerProduct(-m[2][1], s5, m[2][2], s4, -m[2][3], s3)},

                       {s * InnerProduct(-m[1][0], c5, m[1][2], c2, -m[1][3], c1),
                        s * InnerProduct(m[0][0], c5, m[0][3], c1, -m[0][2], c2),
                        s * InnerProduct(-m[3][0], s5, m[3][2], s2, -m[3][3], s1),
                        s * InnerProduct(m[2][0], s5, m[2][3], s1, -m[2][2], s2)},

                       {s * InnerProduct(m[1][0], c4, m[1][3], c0, -m[1][1], c2),
                        s * InnerProduct(-m[0][0], c4, m[0][1], c2, -m[0][3], c0),
                        s * InnerProduct(m[3][0], s4, m[3][3], s0, -m[3][1], s2),
                        s * InnerProduct(-m[2][0], s4, m[2][1], s2, -m[2][3], s0)},

                       {s * InnerProduct(-m[1][0], c3, m[1][1], c1, -m[1][2], c0),
                        s * InnerProduct(m[0][0], c3, m[0][2], c0, -m[0][1], c1),
                        s * InnerProduct(-m[3][0], s3, m[3][1], s1, -m[3][2], s0),
                        s * InnerProduct(m[2][0], s3, m[2][2], s0, -m[2][1], s1)}};

    return SquareMatrix<4>(inv);
}

extern template class SquareMatrix<2>;
extern template class SquareMatrix<3>;
extern template class SquareMatrix<4>;

#endif