#ifndef OMM_H
#define OMM_H

namespace passpredict{
    struct Omm
    {
        char satnum[9];
        double jdsatepoch;
        double jdsatepochF;
        double bstar;
        double inclo; // deg
        double nodeo; // deg
        double ecco;
        double argpo; // deg
        double mo;    // deg
        double no_kozai;
        int revnum;
        int elnum;
        char classification;
        int ephtype;
    };
};

#endif