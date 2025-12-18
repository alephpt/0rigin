#pragma once
#include <Nova/Nova.h>

// Forward declaration of SMFTEngine
class SMFTEngine;

class SMFTCore {
public:
    SMFTCore();
    ~SMFTCore();

    static SMFTCore* manifest();   // Singleton

    static void materialize();      // Draw UI and simulation
    void actualize();               // Run the application

private:
    NovaConfig _config;
    Nova* _actuality;
    SMFTEngine* _smftEngine;

    float dt = 0.01f;
    float K = 1.0f;
    float Delta = 2.5f;
    float chiral_angle = 0.0f;
    float R_avg = 0.0f;

    SMFTCore* realize();            // Init
};
