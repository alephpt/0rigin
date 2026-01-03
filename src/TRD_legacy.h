#pragma once
#include <Nova/Nova.h>

// Forward declaration of TRDEngine
class TRDEngine;

class TRD {
public:
    TRD();
    ~TRD();

    static TRD* manifest();   // Singleton

    static void materialize();      // Draw UI and simulation
    void actualize();               // Run the application

private:
    NovaConfig _config;
    Nova* _actuality;
    TRDEngine* _msftEngine;

    float dt = 0.01f;
    float K = 1.0f;
    float Delta = 2.5f;
    float chiral_angle = 0.0f;
    float R_avg = 0.0f;

    TRD* realize();            // Init
};
