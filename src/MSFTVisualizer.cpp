#include "MSFTVisualizer.h"
#include <cstring>
#include <algorithm>
#include <iostream>

MSFTVisualizer::MSFTVisualizer(Nova* nova) : _nova(nova), _width(0), _height(0) {
    _texture_image = VK_NULL_HANDLE;
    _texture_memory = VK_NULL_HANDLE;
    _texture_image_view = VK_NULL_HANDLE;
    _texture_sampler = VK_NULL_HANDLE;
    _staging_buffer = VK_NULL_HANDLE;
    _staging_memory = VK_NULL_HANDLE;
}

MSFTVisualizer::~MSFTVisualizer() {
    destroyResources();
}

void MSFTVisualizer::initialize(uint32_t width, uint32_t height) {
    _width = width;
    _height = height;

    createStagingBuffer();
    createTextureResources();
}

void MSFTVisualizer::updateRFieldTexture(const std::vector<float>& r_field) {
    if (r_field.size() != _width * _height) {
        std::cerr << "MSFTVisualizer: R field size mismatch!" << std::endl;
        return;
    }

    // Create RGBA texture data from R field
    std::vector<uint8_t> rgba_data(_width * _height * 4);

    for (uint32_t y = 0; y < _height; ++y) {
        for (uint32_t x = 0; x < _width; ++x) {
            uint32_t idx = y * _width + x;
            float r_val = r_field[idx];

            // Clamp to [0,1]
            r_val = std::max(0.0f, std::min(1.0f, r_val));

            // Convert to RGBA color
            uint8_t* pixel = &rgba_data[idx * 4];
            mapRToColor(r_val, pixel);
        }
    }

    // For now, we'll just log that we would update the texture
    // Full implementation would upload rgba_data to GPU texture
    static int frame_count = 0;
    if (frame_count++ % 60 == 0) {  // Log every 60 frames
        float avg_r = 0.0f;
        for (float r : r_field) avg_r += r;
        avg_r /= r_field.size();
        std::cout << "MSFTVisualizer: Average R = " << avg_r << std::endl;
    }
}

void MSFTVisualizer::mapRToColor(float r, uint8_t* rgba) {
    // Heat map: blue (R=0) -> green (R=0.5) -> red (R=1)
    float red, green, blue;

    if (r < 0.5f) {
        // Blue to Green transition
        float t = r * 2.0f;  // Normalize to [0,1]
        red = 0.0f;
        green = t;
        blue = 1.0f - t;
    } else {
        // Green to Red transition
        float t = (r - 0.5f) * 2.0f;  // Normalize to [0,1]
        red = t;
        green = 1.0f - t;
        blue = 0.0f;
    }

    // Convert to uint8_t
    rgba[0] = static_cast<uint8_t>(red * 255.0f);
    rgba[1] = static_cast<uint8_t>(green * 255.0f);
    rgba[2] = static_cast<uint8_t>(blue * 255.0f);
    rgba[3] = 255;  // Full opacity
}

void MSFTVisualizer::createTextureResources() {
    // Simplified for now - full implementation would create Vulkan texture
    std::cout << "MSFTVisualizer: Creating " << _width << "x" << _height << " texture resources" << std::endl;
}

void MSFTVisualizer::createStagingBuffer() {
    // Simplified for now - full implementation would create Vulkan staging buffer
    std::cout << "MSFTVisualizer: Creating staging buffer for " << _width << "x" << _height << " texture" << std::endl;
}

void MSFTVisualizer::destroyResources() {
    // Clean up Vulkan resources (simplified for now)
    if (_texture_image_view != VK_NULL_HANDLE) {
        // vkDestroyImageView(device, _texture_image_view, nullptr);
    }
    if (_texture_image != VK_NULL_HANDLE) {
        // vkDestroyImage(device, _texture_image, nullptr);
    }
    if (_texture_memory != VK_NULL_HANDLE) {
        // vkFreeMemory(device, _texture_memory, nullptr);
    }
    if (_texture_sampler != VK_NULL_HANDLE) {
        // vkDestroySampler(device, _texture_sampler, nullptr);
    }
    if (_staging_buffer != VK_NULL_HANDLE) {
        // vkDestroyBuffer(device, _staging_buffer, nullptr);
    }
    if (_staging_memory != VK_NULL_HANDLE) {
        // vkFreeMemory(device, _staging_memory, nullptr);
    }
}