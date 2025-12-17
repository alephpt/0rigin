#include "MSFTPipelineFactory.h"
#include <fstream>
#include <iostream>
#include <algorithm>

MSFTPipelineFactory::MSFTPipelineFactory(VkDevice device)
    : _device(device) {
    if (device == VK_NULL_HANDLE) {
        throw std::invalid_argument("MSFTPipelineFactory: Invalid device handle");
    }
}

MSFTPipelineFactory::~MSFTPipelineFactory() {
    // Clean up any remaining pipelines
    destroyAllPipelines();
}

std::vector<uint32_t> MSFTPipelineFactory::loadShaderFile(const std::string& path) {
    // Open shader file in binary mode, positioned at end
    std::ifstream file(path, std::ios::ate | std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "MSFTPipelineFactory: Failed to open shader file: " << path << std::endl;
        return {};
    }

    // Get file size from current position (at end)
    size_t fileSize = static_cast<size_t>(file.tellg());
    if (fileSize == 0 || fileSize % sizeof(uint32_t) != 0) {
        std::cerr << "MSFTPipelineFactory: Invalid shader file size: " << fileSize << std::endl;
        return {};
    }

    // Allocate buffer for SPIR-V code
    std::vector<uint32_t> buffer(fileSize / sizeof(uint32_t));

    // Read file from beginning
    file.seekg(0);
    file.read(reinterpret_cast<char*>(buffer.data()), fileSize);
    file.close();

    // Validate SPIR-V magic number (0x07230203)
    if (!buffer.empty() && buffer[0] != 0x07230203) {
        std::cerr << "MSFTPipelineFactory: Invalid SPIR-V magic number in: " << path << std::endl;
        return {};
    }

    return buffer;
}

VkShaderModule MSFTPipelineFactory::createShaderModule(const std::vector<uint32_t>& code) {
    if (code.empty()) {
        return VK_NULL_HANDLE;
    }

    VkShaderModuleCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    createInfo.codeSize = code.size() * sizeof(uint32_t);
    createInfo.pCode = code.data();

    VkShaderModule shaderModule;
    if (vkCreateShaderModule(_device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS) {
        std::cerr << "MSFTPipelineFactory: Failed to create shader module" << std::endl;
        return VK_NULL_HANDLE;
    }

    return shaderModule;
}

VkPipeline MSFTPipelineFactory::createComputePipeline(VkShaderModule shaderModule,
                                                      VkPipelineLayout pipelineLayout) {
    if (shaderModule == VK_NULL_HANDLE || pipelineLayout == VK_NULL_HANDLE) {
        return VK_NULL_HANDLE;
    }

    // Configure shader stage
    VkPipelineShaderStageCreateInfo shaderStageInfo{};
    shaderStageInfo.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    shaderStageInfo.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    shaderStageInfo.module = shaderModule;
    shaderStageInfo.pName = "main";  // Entry point name in shader

    // Create compute pipeline
    VkComputePipelineCreateInfo pipelineInfo{};
    pipelineInfo.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    pipelineInfo.stage = shaderStageInfo;
    pipelineInfo.layout = pipelineLayout;

    VkPipeline pipeline;
    VkResult result = vkCreateComputePipelines(_device, VK_NULL_HANDLE, 1,
                                               &pipelineInfo, nullptr, &pipeline);

    if (result != VK_SUCCESS) {
        std::cerr << "MSFTPipelineFactory: Failed to create compute pipeline, error: "
                  << result << std::endl;
        return VK_NULL_HANDLE;
    }

    // Track the pipeline for cleanup
    _createdPipelines.push_back(pipeline);

    return pipeline;
}

VkPipeline MSFTPipelineFactory::createPipelineFromShader(const std::string& shaderPath,
                                                         VkPipelineLayout pipelineLayout,
                                                         const std::string& pipelineName) {
    // Load shader bytecode
    auto shaderCode = loadShaderFile(shaderPath);
    if (shaderCode.empty()) {
        std::cerr << "MSFTPipelineFactory: Failed to load shader for "
                  << pipelineName << ": " << shaderPath << std::endl;
        return VK_NULL_HANDLE;
    }

    // Create shader module
    VkShaderModule shaderModule = createShaderModule(shaderCode);
    if (shaderModule == VK_NULL_HANDLE) {
        std::cerr << "MSFTPipelineFactory: Failed to create shader module for "
                  << pipelineName << std::endl;
        return VK_NULL_HANDLE;
    }

    // Create pipeline
    VkPipeline pipeline = createComputePipeline(shaderModule, pipelineLayout);

    // Clean up shader module (pipeline keeps internal reference)
    vkDestroyShaderModule(_device, shaderModule, nullptr);

    if (pipeline == VK_NULL_HANDLE) {
        std::cerr << "MSFTPipelineFactory: Failed to create pipeline for "
                  << pipelineName << std::endl;
    }

    return pipeline;
}

VkPipeline MSFTPipelineFactory::createKuramotoPipeline(const std::string& shaderPath,
                                                       VkPipelineLayout pipelineLayout) {
    /**
     * Kuramoto Pipeline - Phase evolution using Kuramoto dynamics
     *
     * GPU SAFETY: ✅ SAFE - 9 transcendentals per workgroup
     * - Expected shader: kuramoto_step.comp
     * - Workload: Well under 20 Tflops budget
     * - Timeout risk: None
     *
     * Shader bindings:
     * - Binding 0: theta_buffer (input phases)
     * - Binding 1: theta_out_buffer (output phases)
     * - Binding 2: omega_buffer (natural frequencies)
     * - Binding 3: spinor_density_buffer (quantum feedback)
     *
     * Push constants: dt, K, damping, Delta, chiral_angle, Nx, Ny, N_total, neighborhood_radius
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "Kuramoto");
}

VkPipeline MSFTPipelineFactory::createSyncFieldPipeline(const std::string& shaderPath,
                                                        VkPipelineLayout pipelineLayout) {
    /**
     * Sync Field Pipeline - Calculate synchronization field R(x)
     *
     * GPU SAFETY: ✅ SAFE - 37 transcendentals (using sync_field_simple.comp)
     * - RECOMMENDED: sync_field_simple.comp (37 transcendentals)
     * - AVOID: sync_field.comp (36+ transcendentals + Kahan summation = borderline timeout)
     * - Workload: Within 20 Tflops budget
     * - Timeout risk: None with simple version
     *
     * Shader bindings:
     * - Binding 0: theta_buffer (input phases)
     * - Binding 1: R_field_buffer (output sync field)
     *
     * The sync field R represents local phase coherence and acts as
     * the effective mass in the MSFT framework (R = √(ℏc/G) * mass density)
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "SyncField");
}

VkPipeline MSFTPipelineFactory::createGravityFieldPipeline(const std::string& shaderPath,
                                                           VkPipelineLayout pipelineLayout) {
    /**
     * Gravity Field Pipeline - Calculate gravitational field g = -Δ·∇R
     *
     * GPU SAFETY: ✅ SAFE - 0 transcendentals (pure arithmetic)
     * - Expected shader: gravity_field.comp
     * - Workload: Minimal - only gradient computation
     * - Timeout risk: None
     *
     * Shader bindings:
     * - Binding 0: R_field_buffer (input sync field)
     * - Binding 1: gravity_x_buffer (output x-component)
     * - Binding 2: gravity_y_buffer (output y-component)
     *
     * Gravity emerges from spatial gradients in the synchronization field,
     * implementing the MSFT principle that gravity IS sync field geometry
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "GravityField");
}

VkPipeline MSFTPipelineFactory::createKuramotoStochasticPipeline(const std::string& shaderPath,
                                                                 VkPipelineLayout pipelineLayout) {
    /**
     * Stochastic Kuramoto Pipeline - Phase evolution with thermal noise
     *
     * GPU SAFETY: ✅ SAFE - 12-14 transcendentals per workgroup
     * - Expected shader: kuramoto_stochastic.comp
     * - Workload: Well under 20 Tflops budget (includes PRNG overhead)
     * - Timeout risk: None
     *
     * Implements Martin-Siggia-Rose (MSR) formalism for stochastic dynamics.
     * Adds white noise to phase evolution for thermal fluctuations.
     *
     * Additional push constants include noise strength and PRNG seed.
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "KuramotoStochastic");
}

VkPipeline MSFTPipelineFactory::createDiracPipeline(const std::string& shaderPath,
                                                    VkPipelineLayout pipelineLayout) {
    /**
     * Dirac Evolution Pipeline - Quantum spinor field evolution
     *
     * GPU SAFETY: ❌ DANGEROUS - ~3000 FLOPs per workgroup (10× over budget)
     * - Expected shader: dirac_rk4.comp
     * - Workload: RK4 integration with 4 stages, complex arithmetic, gamma matrices
     * - Timeout risk: HIGH - consistent 2+ second timeouts on GTX 1650
     * - STATUS: DISABLED in MSFTEngine::createPipelines()
     *
     * Evolves 4-component Dirac spinors according to:
     * (iγ^μ∂_μ)Ψ = [√(ℏc/G)] · R(x) · e^(iθ(x)γ⁵) Ψ
     *
     * Shader bindings include spinor field buffers and mass field R(x).
     *
     * RECOMMENDATION: Use CPU-based implementation or simplify to Euler integration.
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "Dirac");
}

VkPipeline MSFTPipelineFactory::createDiracStochasticPipeline(const std::string& shaderPath,
                                                              VkPipelineLayout pipelineLayout) {
    /**
     * Stochastic Dirac Pipeline - Quantum evolution with vacuum fluctuations
     *
     * GPU SAFETY: ❌ DANGEROUS - 50-80 transcendentals per workgroup (4× over budget)
     * - Expected shader: dirac_stochastic.comp
     * - Workload: Dirac RK4 + PRNG + noise injection + complex arithmetic
     * - Timeout risk: CRITICAL - 4× transcendental budget exceeded
     * - STATUS: DISABLED in MSFTEngine::createPipelines()
     *
     * Adds stochastic noise to Dirac evolution representing quantum
     * vacuum fluctuations and thermal effects in the spinor field.
     *
     * Implements MSR formalism for quantum stochastic dynamics.
     *
     * RECOMMENDATION: CPU implementation mandatory. GPU version requires major simplification.
     */
    return createPipelineFromShader(shaderPath, pipelineLayout, "DiracStochastic");
}

void MSFTPipelineFactory::destroyPipeline(VkPipeline pipeline) {
    if (pipeline != VK_NULL_HANDLE) {
        vkDestroyPipeline(_device, pipeline, nullptr);

        // Remove from tracking list
        auto it = std::find(_createdPipelines.begin(), _createdPipelines.end(), pipeline);
        if (it != _createdPipelines.end()) {
            _createdPipelines.erase(it);
        }
    }
}

void MSFTPipelineFactory::destroyAllPipelines() {
    for (VkPipeline pipeline : _createdPipelines) {
        if (pipeline != VK_NULL_HANDLE) {
            vkDestroyPipeline(_device, pipeline, nullptr);
        }
    }
    _createdPipelines.clear();
}