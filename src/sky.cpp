#include <Rcpp.h>
#include "RcppThread.h"

#include "ext/pbrt/util/buffercache.h"
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

using Allocator = pstd::pmr::polymorphic_allocator<std::byte>;

vec3f latLongSquareToSphere(const vec2f &sample) {
    Float theta = static_cast<Float>(M_PI) * sample.xy.y;     
    Float phi = 2 * static_cast<Float>(M_PI) * sample.xy.x; 
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
    int nThreads = numbercores;
    // int nThreads = pbrt::Options->nThreads != 0 ? pbrt::Options->nThreads : pbrt::AvailableCores();
    pbrt::ParallelInit(nThreads);  
    pbrt::ColorEncoding::Init(Allocator{});
    pbrt::Spectra::Init(Allocator{});
    pbrt::RGBToSpectrumTable::Init(Allocator{});

    pbrt::RGBColorSpace::Init(Allocator{});
    // pbrt::InitBufferCaches();  
    // pbrt::PBRTOptions opt;
    // pbrt::InitPBRT(opt);

    Float elevation_rad = Radians((Float)elevation);

    vec3f sunDir(0., std::cos(elevation_rad), std::sin(elevation_rad));

    // They assert wavelengths are in this range...
    int nLambda = 1 + (720 - 320) / 32;
    std::vector<Float> lambda(nLambda, Float(0));
    for (int i = 0; i < nLambda; ++i) {
        lambda[i] = Lerp(i / Float(nLambda - 1), static_cast<Float>(320), static_cast<Float>(720));
    }

    // Assume a uniform spectral albedo
    ArHosekSkyModelState *skymodel_state =
        arhosekskymodelstate_alloc_init(elevation_rad, turbidity, albedo);

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
                if(!(gamma >= 0 && gamma <= M_PI)) {
                    throw std::runtime_error("gamma out of range");
                }

                for (int i = 0; i < lambda.size(); ++i) {
                    skyv[i] = (Float)arhosekskymodel_solar_radiance(skymodel_state, theta, gamma,
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
    int ret = 0;
    // Free resources
    arhosekskymodelstate_free(skymodel_state);
    pbrt::ParallelCleanup();
    // pbrt::CleanupPBRT();
    return ret;
}
