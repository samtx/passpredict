#ifndef OMM_H
#define OMM_H

namespace passpredict{
    struct Omm
    {
        char classification;
        char satnum[9];
        int revnum;
        int elnum;
        int ephtype;
        double jdsatepoch;
        double jdsatepochF;
        double ndot;  // rad/s
        double nddot;  // rad/s^2
        double bstar;
        double inclo; // deg
        double nodeo; // deg
        double ecco;
        double argpo; // deg
        double mo;    // deg
        double no_kozai;
    };
};

#endif