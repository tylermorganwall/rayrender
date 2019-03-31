// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// render_scene_rcpp
List render_scene_rcpp(int nx, int ny, int ns, float fov, bool ambient_light, NumericVector lookfromvec, NumericVector lookatvec, float aperture, NumericVector camera_up, IntegerVector type, NumericVector radius, IntegerVector shape, NumericVector x, NumericVector y, NumericVector z, List properties, List velocity, LogicalVector moving, int n, NumericVector& bghigh, NumericVector& bglow, float shutteropen, float shutterclose, LogicalVector ischeckered, List checkercolors, NumericVector noise, LogicalVector isnoise, NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist, List& angle, LogicalVector& isimage, CharacterVector& filelocation, LogicalVector& islight, NumericVector& lightintensity, LogicalVector& isflipped, float focus_distance, LogicalVector& isvolume, NumericVector& voldensity, bool parallel, LogicalVector& implicit_sample, List& order_rotation_list, float clampval, LogicalVector& isgrouped, List& group_pivot, List& group_translate, List& group_angle, List& group_order_rotation);
RcppExport SEXP _rayrender_render_scene_rcpp(SEXP nxSEXP, SEXP nySEXP, SEXP nsSEXP, SEXP fovSEXP, SEXP ambient_lightSEXP, SEXP lookfromvecSEXP, SEXP lookatvecSEXP, SEXP apertureSEXP, SEXP camera_upSEXP, SEXP typeSEXP, SEXP radiusSEXP, SEXP shapeSEXP, SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP propertiesSEXP, SEXP velocitySEXP, SEXP movingSEXP, SEXP nSEXP, SEXP bghighSEXP, SEXP bglowSEXP, SEXP shutteropenSEXP, SEXP shuttercloseSEXP, SEXP ischeckeredSEXP, SEXP checkercolorsSEXP, SEXP noiseSEXP, SEXP isnoiseSEXP, SEXP noisephaseSEXP, SEXP noiseintensitySEXP, SEXP noisecolorlistSEXP, SEXP angleSEXP, SEXP isimageSEXP, SEXP filelocationSEXP, SEXP islightSEXP, SEXP lightintensitySEXP, SEXP isflippedSEXP, SEXP focus_distanceSEXP, SEXP isvolumeSEXP, SEXP voldensitySEXP, SEXP parallelSEXP, SEXP implicit_sampleSEXP, SEXP order_rotation_listSEXP, SEXP clampvalSEXP, SEXP isgroupedSEXP, SEXP group_pivotSEXP, SEXP group_translateSEXP, SEXP group_angleSEXP, SEXP group_order_rotationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ny(nySEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< float >::type fov(fovSEXP);
    Rcpp::traits::input_parameter< bool >::type ambient_light(ambient_lightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lookfromvec(lookfromvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lookatvec(lookatvecSEXP);
    Rcpp::traits::input_parameter< float >::type aperture(apertureSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type camera_up(camera_upSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< List >::type properties(propertiesSEXP);
    Rcpp::traits::input_parameter< List >::type velocity(velocitySEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type moving(movingSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type bghigh(bghighSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type bglow(bglowSEXP);
    Rcpp::traits::input_parameter< float >::type shutteropen(shutteropenSEXP);
    Rcpp::traits::input_parameter< float >::type shutterclose(shuttercloseSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ischeckered(ischeckeredSEXP);
    Rcpp::traits::input_parameter< List >::type checkercolors(checkercolorsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type isnoise(isnoiseSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type noisephase(noisephaseSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type noiseintensity(noiseintensitySEXP);
    Rcpp::traits::input_parameter< List >::type noisecolorlist(noisecolorlistSEXP);
    Rcpp::traits::input_parameter< List& >::type angle(angleSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type isimage(isimageSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type filelocation(filelocationSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type islight(islightSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lightintensity(lightintensitySEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type isflipped(isflippedSEXP);
    Rcpp::traits::input_parameter< float >::type focus_distance(focus_distanceSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type isvolume(isvolumeSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type voldensity(voldensitySEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type implicit_sample(implicit_sampleSEXP);
    Rcpp::traits::input_parameter< List& >::type order_rotation_list(order_rotation_listSEXP);
    Rcpp::traits::input_parameter< float >::type clampval(clampvalSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type isgrouped(isgroupedSEXP);
    Rcpp::traits::input_parameter< List& >::type group_pivot(group_pivotSEXP);
    Rcpp::traits::input_parameter< List& >::type group_translate(group_translateSEXP);
    Rcpp::traits::input_parameter< List& >::type group_angle(group_angleSEXP);
    Rcpp::traits::input_parameter< List& >::type group_order_rotation(group_order_rotationSEXP);
    rcpp_result_gen = Rcpp::wrap(render_scene_rcpp(nx, ny, ns, fov, ambient_light, lookfromvec, lookatvec, aperture, camera_up, type, radius, shape, x, y, z, properties, velocity, moving, n, bghigh, bglow, shutteropen, shutterclose, ischeckered, checkercolors, noise, isnoise, noisephase, noiseintensity, noisecolorlist, angle, isimage, filelocation, islight, lightintensity, isflipped, focus_distance, isvolume, voldensity, parallel, implicit_sample, order_rotation_list, clampval, isgrouped, group_pivot, group_translate, group_angle, group_order_rotation));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rayrender_render_scene_rcpp", (DL_FUNC) &_rayrender_render_scene_rcpp, 48},
    {NULL, NULL, 0}
};

RcppExport void R_init_rayrender(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
