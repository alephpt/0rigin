#pragma once
#include <Nova/Nova.h>
#include <vulkan/vulkan.h>
#include <vector>
#include <cstdint>

/**
 * SMFTVisualizer - Handles visualization of SMFT physics fields
 *
 * Responsible for:
 * - Converting R field data to texture
 * - Creating color mapped visualization
 * - Managing texture updates each frame
 */
class SMFTVisualizer {
public:
    SMFTVisualizer(Nova* nova);
    ~SMFTVisualizer();

    /**
     * Initialize visualization resources
     * @param width Grid width
     * @param height Grid height
     */
    void initialize(uint32_t width, uint32_t height);

    /**
     * Update visualization texture with new R field data
     * @param r_field Synchronization field values [0,1]
     */
    void updateRFieldTexture(const std::vector<float>& r_field);

    /**
     * Get the Vulkan texture image for rendering
     */
    VkImage getTextureImage() const { return _texture_image; }

    /**
     * Get the texture image view for binding
     */
    VkImageView getTextureImageView() const { return _texture_image_view; }

private:
    Nova* _nova;
    uint32_t _width, _height;

    // Vulkan resources
    VkImage _texture_image;
    VkDeviceMemory _texture_memory;
    VkImageView _texture_image_view;
    VkSampler _texture_sampler;

    // Staging buffer for texture updates
    VkBuffer _staging_buffer;
    VkDeviceMemory _staging_memory;

    /**
     * Convert R value to RGBA color using heat map
     * R=0 (blue) -> R=0.5 (green) -> R=1 (red)
     */
    void mapRToColor(float r, uint8_t* rgba);

    void createTextureResources();
    void createStagingBuffer();
    void destroyResources();
};