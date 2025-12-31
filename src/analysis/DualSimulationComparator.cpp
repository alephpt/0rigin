/**
 * DualSimulationComparator.cpp
 *
 * Implementation of dual simulation comparison for time dilation measurement
 */

#include "analysis/DualSimulationComparator.h"
#include "Nova.h"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>

DualSimulationComparator::DualSimulationComparator(uint32_t Nx, uint32_t Ny, float Delta)
    : _Nx(Nx), _Ny(Ny), _Delta(Delta) {

    // Create Nova instance for both engines (CPU mode only)
    NovaConfig nova_config = {
        .name = "DualSimComparator",
        .screen = {800, 600},
        .debug_level = "error",
        .dimensions = "2D",
        .camera_type = "fixed",
        .compute = true
    };
    _nova = std::make_unique<Nova>(nova_config);

    // Create dual SMFT engines
    _engine_standard = std::make_unique<SMFTEngine>(_nova.get());
    _engine_timedilation = std::make_unique<SMFTEngine>(_nova.get());

    std::cout << "[DualSimulationComparator] Initialized " << Nx << "x" << Ny
              << " dual engine system (Δ=" << Delta << ")" << std::endl;
}

void DualSimulationComparator::initialize(
    const std::vector<float>& theta_field,
    const std::vector<float>& omega_field,
    float x0, float y0, float sigma) {

    // Initialize both engines with identical parameters
    float chiral_angle = 0.0f;  // Standard chiral coupling

    _engine_standard->initialize(_Nx, _Ny, _Delta, chiral_angle);
    _engine_timedilation->initialize(_Nx, _Ny, _Delta, chiral_angle);

    // Set identical Kuramoto initial conditions
    _engine_standard->setInitialPhases(theta_field);
    _engine_standard->setNaturalFrequencies(omega_field);

    _engine_timedilation->setInitialPhases(theta_field);
    _engine_timedilation->setNaturalFrequencies(omega_field);

    // Initialize identical Dirac wavepackets
    float amplitude = 1.0f;  // Normalized automatically
    _engine_standard->initializeDiracField(x0, y0, sigma, amplitude);
    _engine_timedilation->initializeDiracField(x0, y0, sigma, amplitude);

    // Enable time-dilation mode ONLY for second engine
    auto* dirac_standard = _engine_standard->getDiracEvolutionNonConst();
    auto* dirac_timedilation = _engine_timedilation->getDiracEvolutionNonConst();

    if (dirac_standard && dirac_timedilation) {
        dirac_standard->setTimeDilationMode(false);       // Standard evolution
        dirac_timedilation->setTimeDilationMode(true);    // Time-dilation evolution

        std::cout << "[DualSimulationComparator] Initialized dual Dirac fields:" << std::endl;
        std::cout << "  Engine A: Standard evolution (dτ = dt)" << std::endl;
        std::cout << "  Engine B: Time-dilation evolution (dτ = R·dt)" << std::endl;
    } else {
        std::cerr << "[DualSimulationComparator] ERROR: Failed to initialize Dirac fields" << std::endl;
    }
}

void DualSimulationComparator::step(float dt, float K, float damping, float lambda_coupling) {
    // Step both engines with identical parameters
    // The only difference is time-dilation mode in engine B

    _engine_standard->stepWithDirac(dt, lambda_coupling, 1, K, damping);
    _engine_timedilation->stepWithDirac(dt, lambda_coupling, 1, K, damping);
}

DualSimulationComparator::DualObservables
DualSimulationComparator::computeObservables(double time) const {
    DualObservables obs;
    obs.time = time;

    // Get Dirac evolution objects
    auto* dirac_std = _engine_standard->getDiracEvolution();
    auto* dirac_td = _engine_timedilation->getDiracEvolution();

    if (!dirac_std || !dirac_td) {
        std::cerr << "[DualSimulationComparator] ERROR: Dirac fields not initialized" << std::endl;
        return obs;
    }

    // Standard evolution observables
    obs.norm_standard = dirac_std->getNorm();

    std::vector<float> mass_std = _engine_standard->getMassField();
    float KE_std, PE_std;
    obs.energy_standard = dirac_std->getEnergy(mass_std, KE_std, PE_std);

    obs.phase_standard = computePhaseAccumulation(dirac_std);

    float x_std, y_std;
    dirac_std->getCenterOfMass(x_std, y_std);
    obs.pos_x_standard = x_std;
    obs.pos_y_standard = y_std;

    // Time-dilation evolution observables
    obs.norm_timedilation = dirac_td->getNorm();

    std::vector<float> mass_td = _engine_timedilation->getMassField();
    float KE_td, PE_td;
    obs.energy_timedilation = dirac_td->getEnergy(mass_td, KE_td, PE_td);

    obs.phase_timedilation = computePhaseAccumulation(dirac_td);

    float x_td, y_td;
    dirac_td->getCenterOfMass(x_td, y_td);
    obs.pos_x_timedilation = x_td;
    obs.pos_y_timedilation = y_td;

    // Phase accumulation difference (key observable)
    obs.delta_phase = obs.phase_standard - obs.phase_timedilation;

    // R-field at particle centers
    std::vector<float> R_std = _engine_standard->getSyncField();
    std::vector<float> R_td = _engine_timedilation->getSyncField();

    obs.R_center_standard = dirac_std->getRFieldAtPosition(R_std, x_std, y_std);
    obs.R_center_timedilation = dirac_td->getRFieldAtPosition(R_td, x_td, y_td);

    return obs;
}

double DualSimulationComparator::computePhaseAccumulation(const DiracEvolution* dirac) const {
    /**
     * Compute global phase φ = arg(⟨Ψ|Ψ⟩) = arg(∫Ψ*Ψ dx)
     *
     * For 4-component spinor:
     * ⟨Ψ|Ψ⟩ = Σ_x Σ_α ψ*_α(x) ψ_α(x)
     *
     * This gives the overall phase accumulation of the wavepacket
     */
    if (!dirac) return 0.0;

    std::complex<double> global_amplitude(0.0, 0.0);

    // Sum over all grid points and spinor components
    for (int c = 0; c < 4; c++) {
        const auto& psi_c = dirac->getComponent(c);

        for (size_t i = 0; i < psi_c.size(); i++) {
            // Accumulate ψ*_α · ψ_α
            std::complex<float> psi = psi_c[i];
            global_amplitude += std::complex<double>(
                psi.real() * psi.real() + psi.imag() * psi.imag(),
                2.0 * psi.real() * psi.imag()
            );
        }
    }

    // Extract phase: φ = arg(⟨Ψ|Ψ⟩)
    double phase = std::atan2(global_amplitude.imag(), global_amplitude.real());

    return phase;
}

std::string DualSimulationComparator::getCSVHeader() {
    return "time,norm_std,norm_td,energy_std,energy_td,"
           "phase_std,phase_td,delta_phase,"
           "pos_x_std,pos_y_std,pos_x_td,pos_y_td,"
           "R_center_std,R_center_td";
}

std::string DualSimulationComparator::toCSVLine(const DualObservables& obs) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(8);

    oss << obs.time << ","
        << obs.norm_standard << "," << obs.norm_timedilation << ","
        << obs.energy_standard << "," << obs.energy_timedilation << ","
        << obs.phase_standard << "," << obs.phase_timedilation << ","
        << obs.delta_phase << ","
        << obs.pos_x_standard << "," << obs.pos_y_standard << ","
        << obs.pos_x_timedilation << "," << obs.pos_y_timedilation << ","
        << obs.R_center_standard << "," << obs.R_center_timedilation;

    return oss.str();
}
