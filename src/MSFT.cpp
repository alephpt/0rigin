#include "MSFT.h"
#include "MSFTEngine.h" // Assuming this will be created
#include <imgui.h>

MSFT* _essence = nullptr;

MSFT::MSFT() {
    //report(LOGGER::INFO, "MSFT - Constructing Nova ..");
    assert(_essence == nullptr);
    _actuality = nullptr;
    _msftEngine = nullptr;
}

MSFT::~MSFT() {
    //report(LOGGER::INFO, "MSFT - The End is Nigh ..");
    if (_essence != nullptr) {
        delete _msftEngine;
        delete _actuality;
        _essence = nullptr;
    }
}

MSFT* MSFT::manifest() {
    //report(LOGGER::INFO, "MSFT - Manifesting ..");
    if (_essence == nullptr) {
        MSFT* essence = new MSFT();
        _essence = essence->realize();
    }
    return _essence;
}

MSFT* MSFT::realize() {
    _config = {
        .name = "MSFT Nova Demo",
        .screen = {1920, 1080},
        .debug_level = "debug",
        .dimensions = "3D",
        .camera_type = "orbit",
        .compute = true,
    };
    //report(LOGGER::INFO, "MSFT - Configuring Nova with debug level: %s", _config.debug_level.c_str());
    _actuality = new Nova(_config);
    _actuality->initialized = true;

    // Now that the Nova engine is initialized, we can initialize the MSFT engine.
    _msftEngine = new MSFTEngine(_actuality);
    
    const uint32_t Nx = 128;
    const uint32_t Ny = 128;
    const float Delta = 2.5f;
    const float chiral_angle = 0.0f;
    _msftEngine->initialize(Nx, Ny, Delta, chiral_angle);

    // set initial conditions
    std::vector<float> theta_init(Nx * Ny);
    std::vector<float> omega_init(Nx * Ny);
    std::srand(42);
    for (uint32_t y = 0; y < Ny; y++) {
        for (uint32_t x = 0; x < Nx; x++) {
            uint32_t idx = y * Nx + x;
            theta_init[idx] = 2.0f * M_PI * (float)std::rand() / RAND_MAX;
            omega_init[idx] = 0.5f * ((float)std::rand() / RAND_MAX - 0.5f);
            float wave_x = std::cos(2.0f * M_PI * x / Nx * 2.0f);
            float wave_y = std::cos(2.0f * M_PI * y / Ny * 2.0f);
            theta_init[idx] += 0.2f * wave_x * wave_y;
        }
    }
    _msftEngine->setInitialPhases(theta_init);
    _msftEngine->setNaturalFrequencies(omega_init);

    return this;
}

void MSFT::materialize() {
    //report(LOGGER::VERBOSE, "MSFT - Materializing ..");

    _essence->_msftEngine->step(_essence->dt, _essence->K, 0.1f);
    
    std::vector<float> R_field = _essence->_msftEngine->getSyncField();
    float r_sum = 0.0f;
    for (float r : R_field) {
        r_sum += r;
    }
    _essence->R_avg = r_sum / R_field.size();

    // ImGui window
    ImGui::Begin("MSFT Controls");
    ImGui::Text("Simulation Parameters");
    ImGui::SliderFloat("dt", &_essence->dt, 0.001f, 0.1f);
    ImGui::SliderFloat("K", &_essence->K, 0.1f, 10.0f);
    ImGui::SliderFloat("Delta", &_essence->Delta, 0.1f, 10.0f);
    ImGui::Text("Average R: %f", _essence->R_avg);
    ImGui::End();
}

void MSFT::actualize() {
    //report(LOGGER::INFO, "MSFT - Actualizing ..");

    // The materialize function will be called inside the main loop
    // Nova's illuminate() doesn't take parameters - materialize will be called elsewhere
    _actuality->illuminate();
}
