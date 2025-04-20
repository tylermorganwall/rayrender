#include <Rcpp.h>
#include "RcppThread.h"

#include "math/RayMatrix.h"
#include "utils/assert.h"

#include "math/mathinline.h"
#include "ext/pbrt/util/colorspace.h"
#include "ext/pbrt/util/parallel.h"
#include "ext/pbrt/util/image.h"
#include "ext/pbrt/options.h"

#include "ext/tinyobj/tinyexr.h"


// These headers are C-based, so import them via extern "C":
extern "C" {
#include "ext/sky/ArHosekSkyModel.h"
// #include "sky/ArHosekSkyModelData_Spectral.h"
// #include "sky/ArHosekSkyModelData_RGB.h"
// #include "sky/ArHosekSkyModelData_CIEXYZ.h"
}

int saveEXRImage(const std::string& outfile, const RayMatrix& sky_image, 
                 int width, int height, const pbrt::RGBColorSpace* colorSpace,
                 float albedo, float elevation, float turbidity) {
    // Create EXR header
    EXRHeader header;
    InitEXRHeader(&header);

    // Set compression type
    if ((width < 16) && (height < 16)) {
        // No compression for small image
        header.compression_type = TINYEXR_COMPRESSIONTYPE_NONE;
    } else {
        header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;
    }

    // Create EXR image
    EXRImage image;
    InitEXRImage(&image);

    // 3 channels: RGB
    image.num_channels = 3;
    
    // Split RGB into separate channels
    std::vector<float> R(width * height);
    std::vector<float> G(width * height);
    std::vector<float> B(width * height);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;
            R[idx] = sky_image.value(x, y, 0);
            G[idx] = sky_image.value(x, y, 1);
            B[idx] = sky_image.value(x, y, 2);
        }
    }

    // Setup image pointers
    float* image_ptr[3];
    image_ptr[0] = &R[0]; // B
    image_ptr[1] = &G[0]; // G
    image_ptr[2] = &B[0]; // R

    image.images = reinterpret_cast<unsigned char**>(image_ptr);
    image.width = width;
    image.height = height;

    // Setup EXR header channels
    header.num_channels = 3;
    header.channels = static_cast<EXRChannelInfo*>(malloc(sizeof(EXRChannelInfo) * 3));
    
    // Channel names must be BGR order for most EXR viewers
    strncpy(header.channels[0].name, "R", 255);
    strncpy(header.channels[1].name, "G", 255);
    strncpy(header.channels[2].name, "B", 255);
    
    header.channels[0].name[strlen("R")] = '\0';
    header.channels[1].name[strlen("G")] = '\0';
    header.channels[2].name[strlen("B")] = '\0';

    // Set pixel types
    header.pixel_types = static_cast<int*>(malloc(sizeof(int) * 3));
    header.requested_pixel_types = static_cast<int*>(malloc(sizeof(int) * 3));
    
    for (int i = 0; i < 3; i++) {
        header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // Input image pixel type
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // Output pixel type
    }

    // Add metadata
    header.num_custom_attributes = 3;
    header.custom_attributes = static_cast<EXRAttribute*>(malloc(sizeof(EXRAttribute) * 3));

    // Albedo attribute
    std::string albedo_str = std::to_string(albedo);
    strncpy(header.custom_attributes[0].name, "skymodel.albedo", 255);
    strncpy(header.custom_attributes[0].type, "string", 255);
    header.custom_attributes[0].size = static_cast<int>(albedo_str.size() + 1);
    header.custom_attributes[0].value = static_cast<unsigned char*>(malloc(header.custom_attributes[0].size));
    memcpy(header.custom_attributes[0].value, albedo_str.c_str(), header.custom_attributes[0].size);

    // Elevation attribute
    std::string elevation_str = std::to_string(elevation);
    strncpy(header.custom_attributes[1].name, "skymodel.elevation", 255);
    strncpy(header.custom_attributes[1].type, "string", 255);
    header.custom_attributes[1].size = static_cast<int>(elevation_str.size() + 1);
    header.custom_attributes[1].value = static_cast<unsigned char*>(malloc(header.custom_attributes[1].size));
    memcpy(header.custom_attributes[1].value, elevation_str.c_str(), header.custom_attributes[1].size);

    // Turbidity attribute
    std::string turbidity_str = std::to_string(turbidity);
    strncpy(header.custom_attributes[2].name, "skymodel.turbidity", 255);
    strncpy(header.custom_attributes[2].type, "string", 255);
    header.custom_attributes[2].size = static_cast<int>(turbidity_str.size() + 1);
    header.custom_attributes[2].value = static_cast<unsigned char*>(malloc(header.custom_attributes[2].size));
    memcpy(header.custom_attributes[2].value, turbidity_str.c_str(), header.custom_attributes[2].size);

    // Write to file
    const char* err = NULL;
    int ret = SaveEXRImageToFile(&image, &header, outfile.c_str(), &err);
    
    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            Rcpp::Rcout << "Error saving EXR file: \n" <<  err << "\n";
            FreeEXRErrorMessage(err);
        }
    }

    // Cleanup
    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
    
    for (int i = 0; i < header.num_custom_attributes; i++) {
        free(header.custom_attributes[i].value);
    }
    free(header.custom_attributes);
    
    return ret;
}

using Allocator = pstd::pmr::polymorphic_allocator<std::byte>;

vec3f latLongSquareToSphere(const vec2f &sample) {
    Float theta = (Float)M_PI * sample.xy.y;     
    Float phi = 2 * (Float)M_PI * sample.xy.x; 
    return { std::sin(theta) * std::cos(phi),
             std::sin(theta) * std::sin(phi),
             std::cos(theta)};
}

// [[Rcpp::export]]
int makesky(std::string outfile, 
            double albedo = 0.5,
            double turbidity = 3.,
            double elevation = 10, //measured from the horizon
            unsigned int resolution = 2048, 
            unsigned int numbercores = 1,
            bool square_projection = false) {
    ASSERT(!(albedo < 0. || albedo > 1.));
    ASSERT(!(turbidity < 1.7 || turbidity > 10.));
    ASSERT(!(elevation < 0. || elevation > 90.));
    ASSERT(!(resolution < 1));
    int nThreads = 1;
    // int nThreads = pbrt::Options->nThreads != 0 ? pbrt::Options->nThreads : pbrt::AvailableCores();
    pbrt::ParallelInit(nThreads);  
    pbrt::ColorEncoding::Init(Allocator{});
    // Before RGBColorSpace::Init!
    pbrt::Spectra::Init(Allocator{});
    pbrt::RGBToSpectrumTable::Init(Allocator{});

    pbrt::RGBColorSpace::Init(Allocator{});

    elevation = Radians(elevation);

    vec3f sunDir(0., std::cos(elevation), std::sin(elevation));

    // They assert wavelengths are in this range...
    int nLambda = 1 + (720 - 320) / 32;
    std::vector<Float> lambda(nLambda, Float(0));
    for (int i = 0; i < nLambda; ++i) {
        lambda[i] = Lerp(i / Float(nLambda - 1), 320, 720);
    }

    // Assume a uniform spectral albedo
    ArHosekSkyModelState *skymodel_state =
        arhosekskymodelstate_alloc_init(elevation, turbidity, albedo);

    const pbrt::RGBColorSpace *colorSpace = pbrt::RGBColorSpace::ACES2065_1;
    RcppThread::ThreadPool pool(numbercores);
    RayMatrix sky_image(resolution, resolution, 3);
    pbrt::Image img(pbrt::PixelFormat::Float, 
                    {static_cast<int>(resolution), static_cast<int>(resolution)},
                     {"R", "G", "B"});
    
    // auto worker = [&lambda, resolution, sunDir, 
    //                skymodel_state, colorSpace, &sky_image, square_projection] (int iy) {
    // auto worker = [&] (int iy) {
    pbrt::ParallelFor(0, resolution, [&](int64_t start, int64_t end) {
        std::vector<Float> skyv(lambda.size());
        for (int64_t iy = start; iy < end; ++iy) {
            Float y = ((Float)iy + 0.5f) / resolution;
            for (int ix = 0; ix < resolution; ++ix) {
                Float x = ((Float)ix + 0.5f) / resolution;
                vec3f v;
                if(square_projection) {
                    v = EqualAreaSquareToSphere({x, y});
                } else {
                    v = latLongSquareToSphere({x, y});
                }
                if (v.xyz.z <= 0) {
                    // downward hemisphere
                    continue;
                }

                Float theta = SphericalTheta(v);

                // Compute the angle between the pixel's direction and the sun
                // direction.
                Float gamma = pbrt::SafeACos(dot(v, sunDir));
                // DCHECK(gamma >= 0 && gamma <= Pi);

                for (int i = 0; i < lambda.size(); ++i) {
                    skyv[i] = arhosekskymodel_solar_radiance(skymodel_state, theta, gamma,
                                                                lambda[i]);
                }

                pbrt::PiecewiseLinearSpectrum spec(pstd::MakeConstSpan(lambda),
                                                    pstd::MakeConstSpan(skyv));
                pbrt::XYZ xyz = SpectrumToXYZ(&spec);
                pbrt::RGB rgb = colorSpace->ToRGB(xyz);

                for (int c = 0; c < 3; ++c)
                    img.SetChannel({ix, int(iy)}, c, rgb[c]);
                // for (int c = 0; c < 3; ++c) {
                //     sky_image(ix, int(iy), c) =  rgb[c];
                // }
            }
        }
    });
     pbrt::ImageMetadata metadata;
    metadata.colorSpace = colorSpace;
    std::map<std::string, std::vector<std::string>> stringVectors;
    stringVectors["makesky.albedo"] = {std::to_string(albedo)};
    stringVectors["makesky.elevation"] = {std::to_string(elevation)};
    stringVectors["makesky.turbidity"] = {std::to_string(turbidity)};
    metadata.stringVectors = stringVectors;
    CHECK(img.Write(outfile, metadata));
    // img.Write(outfile);
    // for(size_t j = 0; j < resolution; j++) {
    //     pool.push(worker, j);
    // }
    // pool.join();

    // // Save the generated image to EXR file
    // int ret = saveEXRImage(outfile, sky_image, resolution, resolution, colorSpace, 
    //                       albedo, elevation, turbidity);
    int ret = 0;
    // Free resources
    arhosekskymodelstate_free(skymodel_state);

    return ret;
}
