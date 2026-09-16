#ifdef NOT_CRAN
#include "PreviewDisplay.h"
#include "preview_sky.h"
#include <testthat.h>

namespace {
std::unique_ptr<PreviewDisplay> NativeTestDisplay(RayCamera& cam,Transform& object,Transform& world) {
#ifdef HAS_OIDN
  return std::unique_ptr<PreviewDisplay>(new PreviewDisplay(4,4,false,true,false,10,&cam,
    &object,&world,nullptr,nullptr,nullptr,false,false));
#else
  return std::unique_ptr<PreviewDisplay>(new PreviewDisplay(4,4,false,true,false,10,&cam,&object,&world,false));
#endif
}
}
context("Native preview controls") {
  test_that("viewport input queues bounded actions and respects focus and modifiers") {
    RayrenderGui gui;gui.can_edit=true;
    rimgui_viewport_v1 view{sizeof(view)};view.width=100;view.height=200;view.mouse_x=25;view.mouse_y=50;
    rimgui_input_v1 input{sizeof(input)};input.keyboard_available=1;input.keys_repeated=RIMGUI_KEY_BIT(RIMGUI_KEY_W);
    for(int i=0;i<1000;++i)gui.collect_input(input,view);
    expect_true((gui.keys.size()==1 && gui.keys.front().count==8));
    input.keyboard_available=0;gui.collect_input(input,view);expect_true((gui.keys.empty()));
    input.keyboard_available=1;input.modifiers=RIMGUI_CTRL;gui.collect_input(input,view);expect_true((gui.keys.empty()));
    input.mouse_available=1;input.mouse_clicked=RIMGUI_MOUSE_RIGHT;gui.collect_input(input,view);
    expect_true((gui.pick_pending && !gui.pick_focus));
    expect_true((gui.pick_u==.75f && gui.pick_v==.75f));
    gui.pick_pending=false;input.mouse_available=0;gui.collect_input(input,view);expect_false(gui.pick_pending);
  }
  test_that("standard camera movement, orbit, pitch and reset use the renderer camera") {
    Transform object,world;
    camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,5,0,1,1);
    auto display=NativeTestDisplay(cam,object,world);RayrenderGui gui;
    display->AttachNativeGui(&gui,true,true);
    auto key=[&](unsigned code,unsigned modifiers=0) { gui.keys.push_back({code,modifiers,1});return display->ApplyNativeControls(nullptr); };
    expect_true((key(RIMGUI_KEY_W)));expect_true((std::abs(cam.get_origin()[2]+9.5f)<1e-6));
    expect_true((key(RIMGUI_KEY_S)));expect_true((std::abs(cam.get_origin()[2]+10)<1e-6));
    expect_false(key(RIMGUI_KEY_E));expect_true((gui.movement_speed==2));
    expect_true((key(RIMGUI_KEY_Q)));expect_true((std::abs(cam.get_origin()[1]-10/std::sqrt(101.0))<1e-6 && std::abs(cam.get_origin().length()-10)<1e-6));
    expect_true((key(RIMGUI_KEY_R)));expect_true((gui.movement_speed==1 && cam.get_origin()[1]==0));
    const auto original_direction=cam.get_w();
    expect_true((key(RIMGUI_KEY_W,RIMGUI_SHIFT)));
    expect_true(((cam.get_w()-original_direction).length()>0));
    expect_true((cam.get_origin()[2]==-10));
    key(RIMGUI_KEY_R);key(RIMGUI_KEY_TAB);
    expect_false(gui.orbit);
    const auto direction=cam.get_w();key(RIMGUI_KEY_A);
    expect_true(((cam.get_w()-direction).length()<1e-6));
    expect_true((std::abs(cam.get_origin()[0])>.1));
    key(RIMGUI_KEY_R);key(RIMGUI_KEY_TAB);gui.movement_speed=128;
    key(RIMGUI_KEY_W);expect_true((cam.get_origin()[2]==-10));
    key(RIMGUI_KEY_F);expect_true((display->write_fast_output && gui.fast_preview));
    gui.fast_preview=0;gui.fast_pending=true;
    expect_true((display->ApplyNativeControls(nullptr)));
    expect_false(display->write_fast_output);expect_false(gui.fast_pending);
    gui.fast_preview=1;gui.fast_pending=true;
    expect_true((display->ApplyNativeControls(nullptr)));expect_true((display->write_fast_output));
    key(RIMGUI_KEY_F);expect_false(display->write_fast_output);expect_false(gui.fast_preview);
    expect_false(display->render_requested);key(RIMGUI_KEY_ENTER);expect_true((display->render_requested));
    display->interactive=false;const auto before=cam.get_origin();key(RIMGUI_KEY_Q);
    expect_true(((cam.get_origin()-before).length()==0));
  }
  test_that("date components retain typed values and clamp to real calendar dates") {
    RayrenderGui gui;gui.set_datetime("2024-01-31 21:42:13");
    expect_true((gui.date[0]==2024 && gui.date[1]==1 && gui.date[2]==31));
    expect_true((gui.time[0]==21 && gui.time[1]==42 && gui.time[2]==13));
    gui.date[1]=2;gui.format_datetime();
    expect_true((std::string(gui.datetime)=="2024-02-29 21:42:13"));
    gui.date[0]=2025;gui.format_datetime();
    expect_true((std::string(gui.datetime)=="2025-02-28 21:42:13"));
    gui.date[0]=1900;gui.date[2]=29;gui.format_datetime();expect_true((gui.date[2]==28));
    gui.date[0]=2000;gui.date[2]=29;gui.format_datetime();expect_true((gui.date[2]==29));
    gui.time[0]=99;gui.time[1]=-1;gui.time[2]=60;gui.format_datetime();
    expect_true((std::string(gui.datetime)=="2000-02-29 23:00:59"));
    gui.date[1]=13;gui.date[2]=0;gui.format_datetime();
    expect_true((std::string(gui.datetime)=="2000-12-01 23:00:59"));
  }
  test_that("sky changes wait for release and a safe checkpoint; invalid dates keep rendering") {
    Transform object,world;
    camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,5,0,1,1);
    auto display=NativeTestDisplay(cam,object,world);RayrenderGui gui;
    display->AttachNativeGui(&gui,true,false);
    int updates=0;double elevation=0;
    display->SetSunControls(30,180,[&](double e,double) { ++updates;elevation=e; });
    gui.sun_pending=true;gui.sun_editing=true;gui.sun_elevation=45;gui.sun_azimuth=90;
    display->ApplyNativeSkyControls();expect_true((updates==0));expect_false(display->ConsumeAtmosphereChange());
    gui.sun_elevation=50;gui.sun_editing=false;display->ApplyNativeSkyControls();
    expect_true((updates==1 && elevation==50));expect_true((display->ConsumeAtmosphereChange()));
    display->ApplyNativeSkyControls();expect_true((updates==1));
    display->SetSkyControls(0,0,"bad",[](double,double,const std::string&) { return std::string("Invalid date"); });
    gui.location_pending=true;display->ApplyNativeSkyControls();
    expect_true((gui.sky_error=="Invalid date"));expect_false(display->ConsumeAtmosphereChange());
    display->SetSunControls(50,90,[](double,double) { throw std::runtime_error("Sky failure"); });
    gui.sun_pending=true;bool threw=false;
    try { display->ApplyNativeSkyControls(); }catch(const std::runtime_error&) { threw=true; }
    expect_true((threw));expect_false(display->ConsumeAtmosphereChange());
  }
  test_that("manual sun direction updates the sky and solar disk without changing other lights") {
    Rcpp::List original=Rcpp::List::create(
      Rcpp::List::create(Rcpp::Named("type")="prague",Rcpp::Named("elevation")=30.,Rcpp::Named("azimuth")=180.,Rcpp::Named("rotation")=25.),
      Rcpp::List::create(Rcpp::Named("type")="disk",Rcpp::Named("radiance_spectrum")="sun",Rcpp::Named("direction")=Rcpp::NumericVector::create(0,1,0),Rcpp::Named("rotation")=0.),
      Rcpp::List::create(Rcpp::Named("type")="disk",Rcpp::Named("radiance_spectrum")="moon",Rcpp::Named("direction")=Rcpp::NumericVector::create(0,1,0)));
    Rcpp::List updated=PreviewSunDescriptions(original,0,0,90);
    Rcpp::List sky=updated[0],sun=updated[1],moon=updated[2],before=original[0];
    auto direction=Rcpp::as<Rcpp::NumericVector>(sun["direction"]);
    expect_true((std::abs(direction[0]+1)<1e-12 && std::abs(direction[1])<1e-12 && std::abs(direction[2])<1e-12));
    expect_true((Rcpp::as<double>(sky["azimuth"])==90 && Rcpp::as<double>(sun["rotation"])==25));
    expect_true((Rcpp::as<Rcpp::NumericVector>(moon["direction"])[1]==1));
    expect_true((Rcpp::as<double>(before["azimuth"])==180));
  }
}
#endif
