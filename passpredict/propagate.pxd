cdef extern from "ext/ast2body/ast2Body.h" namespace "ast2Body":
    void pkepler(double r1[3], double v1[3], double& dtsec, double& ndot, double& nddot, double r2[3], double v2[3])
