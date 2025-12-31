// include/physics/GaugeTheory.h
#pragma once
#include "PhysicsTypes.h"
#include <string>

namespace physics {

class GaugeTheory {
public:
    enum class Mechanism { PROCA, STUCKELBERG, MAXWELL };

    virtual ~GaugeTheory() = default;

    // Core interface
    virtual void computePotentials(const float* theta_field, const float* R_field,
                                   int nx, int ny, float dx, float dt) = 0;
    virtual void computeFieldStrengths() = 0;
    virtual FieldTensor getFieldAt(int i, int j) const = 0;

    // Properties
    virtual bool isGaugeInvariant() const = 0;
    virtual Mechanism getMechanism() const = 0;
    virtual std::string getName() const = 0;

    // Energy
    virtual float computeFieldEnergy() const = 0;
};

} // namespace physics
