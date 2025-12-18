# Dirac Coupling Implementation Guide

**Status**: Headers updated ✅, Implementation pending ❌

## What Was Just Completed

### Header Updates (MSFTEngine.h)
Added declarations for:
- `VkDescriptorSetLayout _dirac_descriptor_layout`
- `VkDescriptorSetLayout _spinor_feedback_descriptor_layout`  
- `VkPipelineLayout _dirac_pipeline_layout`
- `VkPipelineLayout _spinor_feedback_pipeline_layout`
- `VkDescriptorSet _dirac_descriptor_set`
- `VkDescriptorSet _spinor_feedback_descriptor_set`
- `VkPipeline _spinor_feedback_pipeline`

### Constructor Updates (MSFTEngine.cpp)
Initialized all new members to `VK_NULL_HANDLE` in constructor initialization list

### Compilation
✅ Code compiles successfully

---

## Implementation Roadmap

This is a **multi-day effort** (~1-2 weeks as estimated). Below is the step-by-step implementation plan.

---

## Phase 1: Create Dirac Descriptor Set Layout (Day 1)

### File: src/MSFTEngine.cpp
### Location: In `initVulkan()` function, after gravity descriptor layout creation

### Code to Add:

```cpp
// ============================================================================
// DIRAC DESCRIPTOR SET LAYOUT
// ============================================================================
// Bindings for dirac_rk4.comp:
// - binding 0: psi buffer (read/write)
// - binding 1: R_field buffer (readonly)

VkDescriptorSetLayoutBinding dirac_bindings[2] = {};

// Binding 0: Spinor field Ψ (read/write)
dirac_bindings[0].binding = 0;
dirac_bindings[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
dirac_bindings[0].descriptorCount = 1;
dirac_bindings[0].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

// Binding 1: R field (readonly)
dirac_bindings[1].binding = 1;
dirac_bindings[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
dirac_bindings[1].descriptorCount = 1;
dirac_bindings[1].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

VkDescriptorSetLayoutCreateInfo dirac_layout_info = {};
dirac_layout_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
dirac_layout_info.bindingCount = 2;
dirac_layout_info.pBindings = dirac_bindings;

result = vkCreateDescriptorSetLayout(device, &dirac_layout_info, nullptr, &_dirac_descriptor_layout);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create Dirac descriptor set layout!" << std::endl;
    return;
}
```

---

## Phase 2: Create Spinor Feedback Descriptor Set Layout (Day 1)

### Code to Add (after Dirac layout):

```cpp
// ============================================================================
// SPINOR FEEDBACK DESCRIPTOR SET LAYOUT
// ============================================================================
// Bindings for spinor_feedback.comp:
// - binding 0: psi buffer (readonly)
// - binding 1: spinor_density buffer (writeonly)

VkDescriptorSetLayoutBinding feedback_bindings[2] = {};

// Binding 0: Spinor field Ψ (readonly)
feedback_bindings[0].binding = 0;
feedback_bindings[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
feedback_bindings[0].descriptorCount = 1;
feedback_bindings[0].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

// Binding 1: Spinor density ρ (writeonly)
feedback_bindings[1].binding = 1;
feedback_bindings[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
feedback_bindings[1].descriptorCount = 1;
feedback_bindings[1].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

VkDescriptorSetLayoutCreateInfo feedback_layout_info = {};
feedback_layout_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
feedback_layout_info.bindingCount = 2;
feedback_layout_info.pBindings = feedback_bindings;

result = vkCreateDescriptorSetLayout(device, &feedback_layout_info, nullptr, &_spinor_feedback_descriptor_layout);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create spinor feedback descriptor set layout!" << std::endl;
    return;
}
```

---

## Phase 3: Create Pipeline Layouts (Day 1)

### Code to Add (after descriptor layouts):

```cpp
// ============================================================================
// DIRAC PIPELINE LAYOUT
// ============================================================================
VkPushConstantRange dirac_push_constant_range = {};
dirac_push_constant_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
dirac_push_constant_range.offset = 0;
dirac_push_constant_range.size = sizeof(PushConstants);  // Reuse existing struct

VkPipelineLayoutCreateInfo dirac_pipeline_layout_info = {};
dirac_pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
dirac_pipeline_layout_info.setLayoutCount = 1;
dirac_pipeline_layout_info.pSetLayouts = &_dirac_descriptor_layout;
dirac_pipeline_layout_info.pushConstantRangeCount = 1;
dirac_pipeline_layout_info.pPushConstantRanges = &dirac_push_constant_range;

result = vkCreatePipelineLayout(device, &dirac_pipeline_layout_info, nullptr, &_dirac_pipeline_layout);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create Dirac pipeline layout!" << std::endl;
    return;
}

// ============================================================================
// SPINOR FEEDBACK PIPELINE LAYOUT
// ============================================================================
VkPushConstantRange feedback_push_constant_range = {};
feedback_push_constant_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
feedback_push_constant_range.offset = 0;
feedback_push_constant_range.size = sizeof(PushConstants);

VkPipelineLayoutCreateInfo feedback_pipeline_layout_info = {};
feedback_pipeline_layout_info.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
feedback_pipeline_layout_info.setLayoutCount = 1;
feedback_pipeline_layout_info.pSetLayouts = &_spinor_feedback_descriptor_layout;
feedback_pipeline_layout_info.pushConstantRangeCount = 1;
feedback_pipeline_layout_info.pPushConstantRanges = &feedback_push_constant_range;

result = vkCreatePipelineLayout(device, &feedback_pipeline_layout_info, nullptr, &_spinor_feedback_pipeline_layout);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create spinor feedback pipeline layout!" << std::endl;
    return;
}
```

---

## Phase 4: Load Shaders and Create Pipelines (Day 2)

### Code to Add (in pipeline creation section):

```cpp
// ============================================================================
// DIRAC PIPELINE
// ============================================================================
// Load dirac_rk4.comp.spv
std::vector<char> dirac_code = _nova->loadShaderSPV("shaders/smft/dirac_rk4.comp.spv");
VkShaderModule dirac_module = _nova->createShaderModule(dirac_code);

VkPipelineShaderStageCreateInfo dirac_shader_stage = {};
dirac_shader_stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
dirac_shader_stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
dirac_shader_stage.module = dirac_module;
dirac_shader_stage.pName = "main";

VkComputePipelineCreateInfo dirac_pipeline_info = {};
dirac_pipeline_info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
dirac_pipeline_info.stage = dirac_shader_stage;
dirac_pipeline_info.layout = _dirac_pipeline_layout;

result = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &dirac_pipeline_info, nullptr, &_dirac_pipeline);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create Dirac pipeline!" << std::endl;
    return;
}

vkDestroyShaderModule(device, dirac_module, nullptr);

// ============================================================================
// SPINOR FEEDBACK PIPELINE
// ============================================================================
std::vector<char> feedback_code = _nova->loadShaderSPV("shaders/smft/spinor_feedback.comp.spv");
VkShaderModule feedback_module = _nova->createShaderModule(feedback_code);

VkPipelineShaderStageCreateInfo feedback_shader_stage = {};
feedback_shader_stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
feedback_shader_stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
feedback_shader_stage.module = feedback_module;
feedback_shader_stage.pName = "main";

VkComputePipelineCreateInfo feedback_pipeline_info = {};
feedback_pipeline_info.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
feedback_pipeline_info.stage = feedback_shader_stage;
feedback_pipeline_info.layout = _spinor_feedback_pipeline_layout;

result = vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &feedback_pipeline_info, nullptr, &_spinor_feedback_pipeline);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to create spinor feedback pipeline!" << std::endl;
    return;
}

vkDestroyShaderModule(device, feedback_module, nullptr);
```

---

## Phase 5: Update Descriptor Pool (Day 2)

### Find and Update:
Search for `VkDescriptorPoolCreateInfo` and increase pool sizes to account for 2 new descriptor sets.

Current: 3 descriptor sets (kuramoto, sync, gravity)
New: 5 descriptor sets (add dirac, spinor_feedback)

```cpp
// Update poolSize count
VkDescriptorPoolSize poolSize = {};
poolSize.type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
poolSize.descriptorCount = 15;  // Increased from 11 to 15
                                 // kuramoto: 4, sync: 2, gravity: 3, dirac: 2, feedback: 2

VkDescriptorPoolCreateInfo poolInfo = {};
poolInfo.maxSets = 5;  // Increased from 3 to 5
```

---

## Phase 6: Allocate and Write Descriptor Sets (Day 3)

### Code to Add (after gravity descriptor set allocation):

```cpp
// ============================================================================
// ALLOCATE DIRAC DESCRIPTOR SET
// ============================================================================
VkDescriptorSetAllocateInfo dirac_alloc_info = {};
dirac_alloc_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
dirac_alloc_info.descriptorPool = _descriptor_pool;
dirac_alloc_info.descriptorSetCount = 1;
dirac_alloc_info.pSetLayouts = &_dirac_descriptor_layout;

result = vkAllocateDescriptorSets(device, &dirac_alloc_info, &_dirac_descriptor_set);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to allocate Dirac descriptor set!" << std::endl;
    return;
}

// Write descriptor set
VkDescriptorBufferInfo dirac_buffer_infos[2] = {};

// Binding 0: Spinor buffer
dirac_buffer_infos[0].buffer = _spinor_buffer;
dirac_buffer_infos[0].offset = 0;
dirac_buffer_infos[0].range = VK_WHOLE_SIZE;

// Binding 1: R field buffer
dirac_buffer_infos[1].buffer = _R_field_buffer;
dirac_buffer_infos[1].offset = 0;
dirac_buffer_infos[1].range = VK_WHOLE_SIZE;

VkWriteDescriptorSet dirac_writes[2] = {};

dirac_writes[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
dirac_writes[0].dstSet = _dirac_descriptor_set;
dirac_writes[0].dstBinding = 0;
dirac_writes[0].descriptorCount = 1;
dirac_writes[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
dirac_writes[0].pBufferInfo = &dirac_buffer_infos[0];

dirac_writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
dirac_writes[1].dstSet = _dirac_descriptor_set;
dirac_writes[1].dstBinding = 1;
dirac_writes[1].descriptorCount = 1;
dirac_writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
dirac_writes[1].pBufferInfo = &dirac_buffer_infos[1];

vkUpdateDescriptorSets(device, 2, dirac_writes, 0, nullptr);

// ============================================================================
// ALLOCATE SPINOR FEEDBACK DESCRIPTOR SET
// ============================================================================
VkDescriptorSetAllocateInfo feedback_alloc_info = {};
feedback_alloc_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
feedback_alloc_info.descriptorPool = _descriptor_pool;
feedback_alloc_info.descriptorSetCount = 1;
feedback_alloc_info.pSetLayouts = &_spinor_feedback_descriptor_layout;

result = vkAllocateDescriptorSets(device, &feedback_alloc_info, &_spinor_feedback_descriptor_set);
if (result != VK_SUCCESS) {
    std::cerr << "ERROR: Failed to allocate spinor feedback descriptor set!" << std::endl;
    return;
}

VkDescriptorBufferInfo feedback_buffer_infos[2] = {};

// Binding 0: Spinor buffer (readonly)
feedback_buffer_infos[0].buffer = _spinor_buffer;
feedback_buffer_infos[0].offset = 0;
feedback_buffer_infos[0].range = VK_WHOLE_SIZE;

// Binding 1: Spinor density buffer (writeonly)
feedback_buffer_infos[1].buffer = _spinor_density_buffer;
feedback_buffer_infos[1].offset = 0;
feedback_buffer_infos[1].range = VK_WHOLE_SIZE;

VkWriteDescriptorSet feedback_writes[2] = {};

feedback_writes[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
feedback_writes[0].dstSet = _spinor_feedback_descriptor_set;
feedback_writes[0].dstBinding = 0;
feedback_writes[0].descriptorCount = 1;
feedback_writes[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
feedback_writes[0].pBufferInfo = &feedback_buffer_infos[0];

feedback_writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
feedback_writes[1].dstSet = _spinor_feedback_descriptor_set;
feedback_writes[1].dstBinding = 1;
feedback_writes[1].descriptorCount = 1;
feedback_writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
feedback_writes[1].pBufferInfo = &feedback_buffer_infos[1];

vkUpdateDescriptorSets(device, 2, feedback_writes, 0, nullptr);
```

---

## Phase 7: Initialize Spinor Field (Day 3)

### File: src/MSFTEngine.cpp
### Location: In `initVulkan()` after buffer creation

### Code to Add:

```cpp
// ============================================================================
// INITIALIZE SPINOR FIELD Ψ(x,y,0)
// ============================================================================
// Initialize with Gaussian wavepacket centered at grid center

float x0 = _Nx / 2.0f;
float y0 = _Ny / 2.0f;
float sigma = _Nx / 8.0f;  // Width of wavepacket

for (uint32_t y = 0; y < _Ny; y++) {
    for (uint32_t x = 0; x < _Nx; x++) {
        float dx = x - x0;
        float dy = y - y0;
        float r2 = dx*dx + dy*dy;
        
        // Gaussian envelope
        float amplitude = std::exp(-r2 / (2.0f * sigma * sigma));
        
        // Normalize to unit probability
        amplitude /= std::sqrt(2.0f * M_PI * sigma * sigma);
        
        // 4-component spinor: start with spin-up state (1,0,0,0)
        // This is eigenstate of σ³ with eigenvalue +1
        for (int component = 0; component < 4; component++) {
            uint32_t spinor_idx = (y * _Nx + x) * 4 + component;
            
            if (component == 0) {
                _spinor_field[spinor_idx] = amplitude;  // ψ₁ = amplitude
            } else {
                _spinor_field[spinor_idx] = 0.0f;  // ψ₂,₃,₄ = 0
            }
        }
    }
}

// Upload to GPU
void* spinor_data;
vkMapMemory(device, _spinor_memory, 0, spinor_buffer_size, 0, &spinor_data);
memcpy(spinor_data, _spinor_field.data(), spinor_buffer_size);
vkUnmapMemory(device, _spinor_memory);

std::cout << "Initialized spinor field with Gaussian wavepacket (σ=" << sigma << ")" << std::endl;
```

---

## Phase 8: Add Dispatch Calls to step() (Day 4)

### File: src/MSFTEngine.cpp  
### Location: Line 497, after gravity_field dispatch, before vkEndCommandBuffer

### Code to Add:

```cpp
// ============================================================================
// 13. DIRAC EVOLUTION: i∂ₜΨ = [cα·p + βm(x)c²]Ψ
// ============================================================================
vkCmdBindPipeline(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, _dirac_pipeline);
vkCmdBindDescriptorSets(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                       _dirac_pipeline_layout, 0, 1, &_dirac_descriptor_set, 0, nullptr);
vkCmdPushConstants(_compute_command_buffer, _dirac_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                  0, sizeof(pushConstants), &pushConstants);
vkCmdDispatch(_compute_command_buffer, workgroupsX, workgroupsY, 1);

// Memory barrier: Ensure dirac_evolution completes before spinor_feedback
VkMemoryBarrier barrier = {};
barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
vkCmdPipelineBarrier(_compute_command_buffer,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    0, 1, &barrier, 0, nullptr, 0, nullptr);

// ============================================================================
// 14. SPINOR FEEDBACK: Compute ρ(x) = Ψ̄Ψ
// ============================================================================
vkCmdBindPipeline(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, _spinor_feedback_pipeline);
vkCmdBindDescriptorSets(_compute_command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE,
                       _spinor_feedback_pipeline_layout, 0, 1, &_spinor_feedback_descriptor_set, 0, nullptr);
vkCmdPushConstants(_compute_command_buffer, _spinor_feedback_pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT,
                  0, sizeof(pushConstants), &pushConstants);
vkCmdDispatch(_compute_command_buffer, workgroupsX, workgroupsY, 1);

// Memory barrier: Ensure spinor_feedback completes before next kuramoto_step
vkCmdPipelineBarrier(_compute_command_buffer,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    0, 1, &barrier, 0, nullptr, 0, nullptr);

// 15. End command buffer recording
```

**IMPORTANT**: Change the comment number from "12. End command buffer" to "15. End command buffer"

---

## Phase 9: Add Cleanup Code (Day 4)

### File: src/MSFTEngine.cpp
### Location: In destructor `~MSFTEngine()`

### Code to Add (before existing cleanup):

```cpp
// Dirac pipeline cleanup
if (_dirac_pipeline != VK_NULL_HANDLE) {
    vkDestroyPipeline(device, _dirac_pipeline, nullptr);
    _dirac_pipeline = VK_NULL_HANDLE;
}

if (_spinor_feedback_pipeline != VK_NULL_HANDLE) {
    vkDestroyPipeline(device, _spinor_feedback_pipeline, nullptr);
    _spinor_feedback_pipeline = VK_NULL_HANDLE;
}

// Dirac pipeline layout cleanup
if (_dirac_pipeline_layout != VK_NULL_HANDLE) {
    vkDestroyPipelineLayout(device, _dirac_pipeline_layout, nullptr);
    _dirac_pipeline_layout = VK_NULL_HANDLE;
}

if (_spinor_feedback_pipeline_layout != VK_NULL_HANDLE) {
    vkDestroyPipelineLayout(device, _spinor_feedback_pipeline_layout, nullptr);
    _spinor_feedback_pipeline_layout = VK_NULL_HANDLE;
}

// Dirac descriptor layout cleanup
if (_dirac_descriptor_layout != VK_NULL_HANDLE) {
    vkDestroyDescriptorSetLayout(device, _dirac_descriptor_layout, nullptr);
    _dirac_descriptor_layout = VK_NULL_HANDLE;
}

if (_spinor_feedback_descriptor_layout != VK_NULL_HANDLE) {
    vkDestroyDescriptorSetLayout(device, _spinor_feedback_descriptor_layout, nullptr);
    _spinor_feedback_descriptor_layout = VK_NULL_HANDLE;
}
```

---

## Phase 10: Testing (Days 5-7)

### Test 1: Compilation
```bash
make MSFT
```
Expected: Clean build with no errors

### Test 2: Initialization
```bash
./bin/MSFT
```
Expected: No Vulkan validation errors, spinor field initialized

### Test 3: Single Step Execution
Add debug output to check if dispatches execute

### Test 4: Spinor Density Verification
Read back `spinor_density_buffer` and verify ρ > 0

### Test 5: Conservation Check
Verify Ψ̄Ψ is approximately conserved over time

---

## Timeline Summary

| Day | Task | Hours |
|-----|------|-------|
| 1 | Phase 1-3: Descriptor layouts, pipeline layouts | 4-6 |
| 2 | Phase 4-5: Shader loading, pipeline creation, pool update | 4-6 |
| 3 | Phase 6-7: Descriptor sets, spinor initialization | 4-6 |
| 4 | Phase 8-9: Dispatch calls, cleanup code | 3-4 |
| 5-7 | Phase 10: Testing, debugging, validation | 12-16 |

**Total**: ~30-40 hours of focused work

---

## Success Criteria

✅ **Minimum**:
- Code compiles without errors
- No Vulkan validation errors
- Spinor density ρ > 0 after first step
- No crashes during 100-step simulation

✅ **Full**:
- Ψ̄Ψ conserved to within 1%
- Feedback visible in Kuramoto dynamics (R changes)
- Zitterbewegung observable in ρ(x,t)
- Stable for 10k steps

---

## Current Status

- ✅ Headers declared
- ✅ Constructor updated
- ✅ Compilation successful
- ❌ Pipelines not created yet
- ❌ Descriptor sets not allocated yet
- ❌ Dispatch calls not added yet

**Next immediate step**: Implement Phase 1 (Dirac descriptor set layout)
