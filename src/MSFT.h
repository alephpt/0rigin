#pragma once
#include <Nova/Nova.h>

// Forward declaration of MSFTEngine
class MSFTEngine;

class MSFT {
public:
    MSFT();
    ~MSFT();

    static MSFT* manifest();   // Singleton

    static void materialize();      // Draw UI and simulation
    void actualize();               // Run the application

private:
    NovaConfig _config;
    Nova* _actuality;
    MSFTEngine* _msftEngine;

    float dt = 0.01f;
    float K = 1.0f;
    float Delta = 2.5f;
    float chiral_angle = 0.0f;
    float R_avg = 0.0f;

    MSFT* realize();            // Init
};
