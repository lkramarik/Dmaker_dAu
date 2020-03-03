#ifndef StAnaCuts_h
#define StAnaCuts_h

using namespace std;

namespace vars {
//    const bool ratioHists = true;
    const bool ratioHists = false;
//    const bool dcaHists = false;
    const bool dcaHists = true;

//histo stuff
    const int m_nPtTOF = 19;
    float const m_PtTOFedge[m_nPtTOF+1] = {0.15,0.3,0.4,0.5,0.6,0.8,1.,1.25,1.5,2.,3.0,4.,5.,6.,7,8,9,10,11,12.};

//    const int m_nZdc = 10;
//    float const m_zdcEdge[m_nZdc+1] = {0,50,70,90,110,130,150,170,190,210,250};

//    const int m_nZdc = 5;
//    float const m_zdcEdge[m_nZdc+1] = {0,50,90,130,170,210};

    const int m_nZdc = 2;
    float const m_zdcEdge[m_nZdc+1] = {0,150,210};

    const int m_nZdcDCA = 2;
    float const m_zdcEdgeDCA[m_nZdcDCA+1] = {0,150,210};

    const int m_nParticles = 2;
    const TString m_ParticleName[m_nParticles] = {"Pion", "Kaon"};

    const int m_nEtasDca = 3;
    float const m_EtaEdgeDca[m_nEtasDca+1] = {0, 0.4, 0.8, 1.0};

    const int m_nPhisDca = 11;
    float const m_PhiEdgeDca[m_nPhisDca + 1] = {-3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};

    const int m_nVzsDca = 4;
    float const m_VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};

//    const int m_nPtsDca = 19;
//    float const m_PtEdgeDca[m_nPtsDca + 1] = {0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 12.0};

    const int m_nPtsDca = 12;
    float const m_PtEdgeDca[m_nPtsDca + 1] = {0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1., 1.25, 1.5, 2., 3., 5., 12.};

    const int m_nEtasRatio = 10;
    float const m_EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0};

    const int m_nPhisRatio = 11;
    float const m_PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};

    const int m_nVzsRatio = 6;
    float const m_VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};

//    const int m_nmultEdge = 1; //7
//    float const m_multEdge[m_nmultEdge+1] = {0, 200}; //currently not used in dca

    const int m_nmultEdge = 7;
    float const m_multEdge[m_nmultEdge+1] = {0, 4, 8, 12, 16, 20, 24, 200}; //currently not used in dca
//    int multEdge[nmultEdge + 1] = {0, 4, 8, 12, 16, 20, 24, 200};


    const int m_nPtsRatio = 15;
    float const m_PtEdgeRatio[m_nPtsRatio + 1] =
            {
                    0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0, 12.0
            };

//    const int m_nPtsRatio = 36;
//    float const m_PtEdgeRatio[m_nPtsRatio + 1] =  {
//            0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
//            1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
//            2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
//            3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
//            6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11,  12.0 };

    const int m_nDcasDca = 148;
    float const  m_DcaEdgeDca[m_nDcasDca + 1] =   {
            -1.5, -1.2, -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
            -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
            -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
            0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
            0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1, 1.2, 1.5
    };

}


#endif //SIMINPUTS_STANACUTS_H