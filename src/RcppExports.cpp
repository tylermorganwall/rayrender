// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// has_gui_capability
bool has_gui_capability();
RcppExport SEXP _rayrender_has_gui_capability() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(has_gui_capability());
    return rcpp_result_gen;
END_RCPP
}
// render_animation_rcpp
void render_animation_rcpp(List scene, List camera_info, List scene_info, List camera_movement, int start_frame, int end_frame, CharacterVector filenames, Function post_process_frame, int toneval, bool bloom, bool write_image, bool transparent_background);
RcppExport SEXP _rayrender_render_animation_rcpp(SEXP sceneSEXP, SEXP camera_infoSEXP, SEXP scene_infoSEXP, SEXP camera_movementSEXP, SEXP start_frameSEXP, SEXP end_frameSEXP, SEXP filenamesSEXP, SEXP post_process_frameSEXP, SEXP tonevalSEXP, SEXP bloomSEXP, SEXP write_imageSEXP, SEXP transparent_backgroundSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type scene(sceneSEXP);
    Rcpp::traits::input_parameter< List >::type camera_info(camera_infoSEXP);
    Rcpp::traits::input_parameter< List >::type scene_info(scene_infoSEXP);
    Rcpp::traits::input_parameter< List >::type camera_movement(camera_movementSEXP);
    Rcpp::traits::input_parameter< int >::type start_frame(start_frameSEXP);
    Rcpp::traits::input_parameter< int >::type end_frame(end_frameSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type filenames(filenamesSEXP);
    Rcpp::traits::input_parameter< Function >::type post_process_frame(post_process_frameSEXP);
    Rcpp::traits::input_parameter< int >::type toneval(tonevalSEXP);
    Rcpp::traits::input_parameter< bool >::type bloom(bloomSEXP);
    Rcpp::traits::input_parameter< bool >::type write_image(write_imageSEXP);
    Rcpp::traits::input_parameter< bool >::type transparent_background(transparent_backgroundSEXP);
    render_animation_rcpp(scene, camera_info, scene_info, camera_movement, start_frame, end_frame, filenames, post_process_frame, toneval, bloom, write_image, transparent_background);
    return R_NilValue;
END_RCPP
}
// render_scene_rcpp
List render_scene_rcpp(List scene, List camera_info, List scene_info);
RcppExport SEXP _rayrender_render_scene_rcpp(SEXP sceneSEXP, SEXP camera_infoSEXP, SEXP scene_infoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type scene(sceneSEXP);
    Rcpp::traits::input_parameter< List >::type camera_info(camera_infoSEXP);
    Rcpp::traits::input_parameter< List >::type scene_info(scene_infoSEXP);
    rcpp_result_gen = Rcpp::wrap(render_scene_rcpp(scene, camera_info, scene_info));
    return rcpp_result_gen;
END_RCPP
}
// PrintClassSizes
void PrintClassSizes();
RcppExport SEXP _rayrender_PrintClassSizes() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    PrintClassSizes();
    return R_NilValue;
END_RCPP
}
// tonemap_image
Rcpp::List tonemap_image(Rcpp::NumericMatrix routput, Rcpp::NumericMatrix goutput, Rcpp::NumericMatrix boutput, int toneval);
RcppExport SEXP _rayrender_tonemap_image(SEXP routputSEXP, SEXP goutputSEXP, SEXP boutputSEXP, SEXP tonevalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type routput(routputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type goutput(goutputSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type boutput(boutputSEXP);
    Rcpp::traits::input_parameter< int >::type toneval(tonevalSEXP);
    rcpp_result_gen = Rcpp::wrap(tonemap_image(routput, goutput, boutput, toneval));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rayrender_has_gui_capability", (DL_FUNC) &_rayrender_has_gui_capability, 0},
    {"_rayrender_render_animation_rcpp", (DL_FUNC) &_rayrender_render_animation_rcpp, 12},
    {"_rayrender_render_scene_rcpp", (DL_FUNC) &_rayrender_render_scene_rcpp, 3},
    {"_rayrender_PrintClassSizes", (DL_FUNC) &_rayrender_PrintClassSizes, 0},
    {"_rayrender_tonemap_image", (DL_FUNC) &_rayrender_tonemap_image, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_rayrender(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
