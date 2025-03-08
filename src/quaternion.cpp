#include "quaternion.h"
#include "transform.h"

// Quaternion Method Definitions
Transform Quaternion::ToTransform() const {
  Float xx = v.x * v.x, yy = v.y * v.y, zz = v.z * v.z;
  Float xy = v.x * v.y, xz = v.x * v.z, yz = v.y * v.z;
  Float wx = v.x * w, wy = v.y * w, wz = v.z * w;
  
  Matrix4x4 m;
  m.m[0][0] = 1 - 2 * (yy + zz);
  m.m[0][1] = 2 * (xy + wz);
  m.m[0][2] = 2 * (xz - wy);
  m.m[1][0] = 2 * (xy - wz);
  m.m[1][1] = 1 - 2 * (xx + zz);
  m.m[1][2] = 2 * (yz + wx);
  m.m[2][0] = 2 * (xz + wy);
  m.m[2][1] = 2 * (yz - wx);
  m.m[2][2] = 1 - 2 * (xx + yy);
  
  // Transpose since we are left-handed.  Ugh.
  return Transform(Transpose(m), m);
}

Quaternion::Quaternion(const Transform &t) {
  const Matrix4x4 &m = t.m;
  Float trace = m.m[0][0] + m.m[1][1] + m.m[2][2];
  if (trace > 0.f) {
    // Compute w from matrix trace, then xyz
    // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
    Float s = std::sqrt(trace + static_cast<Float>(1));
    w = s / 2.0f;
    s = 0.5f / s;
    v[0] = (m.m[2][1] - m.m[1][2]) * s;
    v[1] = (m.m[0][2] - m.m[2][0]) * s;
    v[2] = (m.m[1][0] - m.m[0][1]) * s;
  } else {
    // Compute largest of $x$, $y$, or $z$, then remaining components
    const int nxt[3] = {1, 2, 0};
    Float q[3];
    int i = 0;
    if (m.m[1][1] > m.m[0][0]) i = 1;
    if (m.m[2][2] > m.m[i][i]) i = 2;
    int j = nxt[i];
    int k = nxt[j];
    Float s = std::sqrt((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + static_cast<Float>(1));
    q[i] = s * 0.5f;
    if (s != 0.f) s = 0.5f / s;
    w = (m.m[k][j] - m.m[j][k]) * s;
    q[j] = (m.m[j][i] + m.m[i][j]) * s;
    q[k] = (m.m[k][i] + m.m[i][k]) * s;
    v[0] = q[0];
    v[1] = q[1];
    v[2] = q[2];
  }
}


Quaternion::Quaternion(const Matrix4x4 &m) {
  Float trace = m.m[0][0] + m.m[1][1] + m.m[2][2];
  if (trace > 0.f) {
    // Compute w from matrix trace, then xyz
    // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
    Float s = std::sqrt(trace + static_cast<Float>(1));
    w = s / 2.0f;
    s = 0.5f / s;
    v[0] = (m.m[2][1] - m.m[1][2]) * s;
    v[1] = (m.m[0][2] - m.m[2][0]) * s;
    v[2] = (m.m[1][0] - m.m[0][1]) * s;
  } else {
    // Compute largest of $x$, $y$, or $z$, then remaining components
    const int nxt[3] = {1, 2, 0};
    Float q[3];
    int i = 0;
    if (m.m[1][1] > m.m[0][0]) i = 1;
    if (m.m[2][2] > m.m[i][i]) i = 2;
    int j = nxt[i];
    int k = nxt[j];
    Float s = std::sqrt((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + static_cast<Float>(1));
    q[i] = s * 0.5f;
    if (s != 0.f) s = 0.5f / s;
    w = (m.m[k][j] - m.m[j][k]) * s;
    q[j] = (m.m[j][i] + m.m[i][j]) * s;
    q[k] = (m.m[k][i] + m.m[i][k]) * s;
    v[0] = q[0];
    v[1] = q[1];
    v[2] = q[2];
  }
}

Quaternion Slerp(Float t, const Quaternion &q1, const Quaternion &q2) {
  Float cosTheta = dot(q1, q2);
  if (cosTheta > .9995f) {
    return Normalize((1 - t) * q1 + t * q2);
  } else {
    Float theta = std::acos(clamp(cosTheta, -1, 1));
    Float thetap = theta * t;
    Quaternion qperp = Normalize(q2 - q1 * cosTheta);
    return q1 * std::cos(thetap) + qperp * std::sin(thetap);
  }
}