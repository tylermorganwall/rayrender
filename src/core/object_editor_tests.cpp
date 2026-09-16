#ifdef NOT_CRAN
#include "preview_scene.h"
#include "buildscene.h"
#include "PreviewDisplay.h"
#include "../hitables/hitablelist.h"
#include "../materials/material.h"
#include "../materials/texturecache.h"
#include "../math/transformcache.h"
#include "../volumes/boundary.h"
#include "../volumes/lights.h"
#include <testthat.h>

namespace {
Rcpp::List EditorFixture(const char* code) {
  Rcpp::Function parse=Rcpp::Environment::base_env()["parse"],eval=Rcpp::Environment::base_env()["eval"];
  Rcpp::Function newenv=Rcpp::Environment::base_env()["new.env"];
  Rcpp::Environment environment=newenv(Rcpp::_["parent"]=Rcpp::Environment::namespace_env("rayrender"));
  return Rcpp::as<Rcpp::List>(eval(parse(Rcpp::_["text"]=code),environment));
}
struct EditorBuild {
  TransformCache transforms;TextureCache texture_cache;
  std::vector<Float*> textures;
  std::vector<unsigned char*> alpha_textures,bump_textures,roughness_textures;
  std::vector<std::shared_ptr<material>> materials;
  std::vector<std::shared_ptr<alpha_texture>> alpha;
  std::vector<std::shared_ptr<bump_texture>> bump;
  std::vector<std::shared_ptr<roughness_texture>> roughness;
  std::vector<std::shared_ptr<hitable>> instances;
  std::vector<std::shared_ptr<hitable_list>> instance_lights;
  std::vector<int> texture_idx;
  hitable_list lights;
  std::shared_ptr<hitable> root;
  EditorBuild(Rcpp::List input,PreviewScene& editor) {
    Rcpp::IntegerVector shape=input["shape"];random_gen rng(1);
    lights.volume_scene=std::make_shared<VolumeScene>();
    root=build_scene(input,shape,0,1,textures,alpha_textures,bump_textures,roughness_textures,&materials,
      alpha,bump,roughness,0,transforms,texture_cache,lights,instances,instance_lights,texture_idx,false,rng,&editor);
    lights.volume_scene->light_sampler=std::make_shared<VolumeLightSampler>(lights);
  }
};
point3f EditorColor(PreviewScene& scene,size_t root) {
  auto fields=PreviewMaterialFields(PreviewMaterials(scene.roots[root].object.get())[0]);
  return point3f(fields[0].field.values[0],fields[0].field.values[1],fields[0].field.values[2]);
}
}
context("Native object editor") {
  test_that("top-level groups and instance placements pick and rebuild independently") {
    auto input=EditorFixture(R"(
      native_editor_scene(process_scene(add_object(
        group_objects(add_object(sphere(x=-3),sphere(x=-1)), translate=c(0,0,0)),
        create_instances(sphere(material=diffuse(color='red')), x=c(2,5))
      ))$scene)
    )");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    editor.prepare_rebuild=[&](PreviewScene& next) { auto staged=std::make_shared<EditorBuild>(input,next);return [&,staged]() { live=staged; }; };
    expect_true((editor.roots.size()==4));
    expect_true((editor.roots[0].id==editor.roots[1].id && editor.roots[2].id!=editor.roots[3].id));
    PreviewObjectState ui;
    expect_true((editor.Pick(Ray(point3f(-3,0,-8),vec3f(0,0,1),.5),ui)));
    expect_true((ui.id==editor.roots[0].id && ui.materials.size()==2));
    ui.translation[0]+=1;ui.numeric_transform=ui.apply_transform=true;
    expect_true((editor.Apply(ui)));expect_true((editor.settings.transforms.size()==2));
    aabb bound;editor.roots[0].object->bounding_box(0,1,bound);
    expect_true((std::abs(bound.Centroid()[0]+2)<1e-5));
    expect_true((editor.Pick(Ray(point3f(2,0,-8),vec3f(0,0,1),.5),ui)));
    expect_true((ui.id==editor.roots[2].id && ui.materials.size()==1));
    expect_true((PreviewMaterials(editor.roots[2].object.get())[0]==PreviewMaterials(editor.roots[3].object.get())[0]));
    ui.materials[0].fields[0].values={0,1,0};ui.materials[0].changed=true;ui.apply_material=true;
    expect_true((editor.Apply(ui)));
    expect_true((PreviewMaterials(editor.roots[2].object.get())[0]!=PreviewMaterials(editor.roots[3].object.get())[0]));
    expect_true(((EditorColor(editor,2)-point3f(0,1,0)).length()<1e-6));
    expect_true(((EditorColor(editor,3)-point3f(1,0,0)).length()<1e-6));
    expect_true((editor.ExportEdits().size()==3));
    ui.revert=true;expect_true((editor.Apply(ui)));
    expect_true(((EditorColor(editor,2)-point3f(1,0,0)).length()<1e-6));
    expect_true((editor.ExportEdits().size()==2));
    expect_false(editor.Pick(Ray(point3f(100,0,-8),vec3f(0,0,1),.5),ui));
    expect_false(ui.selected);
  }
  test_that("invalid transforms and failed rebuilds leave the live scene intact") {
    auto input=EditorFixture("native_editor_scene(process_scene(sphere())$scene)");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    PreviewObjectState ui;editor.Describe(editor.roots[0].id,ui);
    bool published=false;
    editor.prepare_rebuild=[&](PreviewScene& next)->std::function<void()> { throw std::runtime_error("rebuild failed"); };
    auto original=editor.roots[0].object;
    ui.scale[0]=0;ui.numeric_transform=ui.apply_transform=true;
    expect_false(editor.Apply(ui));expect_true((!ui.error.empty() && editor.settings.transforms.empty()));
    ui.translation[0]=7;ui.numeric_transform=ui.apply_transform=true;
    expect_false(editor.Apply(ui));expect_true((ui.error=="rebuild failed" && editor.roots[0].object==original));
    ui.model[0]=std::numeric_limits<float>::quiet_NaN();ui.apply_transform=true;
    expect_false(editor.Apply(ui));expect_true((editor.roots[0].object==original));
    ui.model[12]=3;ui.transform_pending=ui.cancel_transform=true;editor.CancelTransform(ui);
    expect_true((ui.model[12]==0 && !ui.transform_pending && !ui.cancel_transform));
    expect_false(published);
  }
  test_that("mesh selection keeps its root and light edits rebuild light sampling") {
    auto input=EditorFixture(R"(
      path = tempfile(fileext='.obj')
      writeLines(c('v -1 -1 0','v 1 -1 0','v 1 1 0','v -1 1 0','f 1 2 3','f 1 3 4'),path)
      native_editor_scene(process_scene(add_object(
        obj_model(path, material=diffuse(color='blue'),load_material=FALSE),
        sphere(x=4,material=light(color='red',intensity=2))
      ))$scene)
    )");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    editor.prepare_rebuild=[&](PreviewScene& next) { auto staged=std::make_shared<EditorBuild>(input,next);return [&,staged]() { live=staged; }; };
    expect_true((editor.roots.size()==2));
    PreviewObjectState ui;expect_true((editor.Pick(Ray(point3f(0,0,-10),vec3f(0,0,1),.5),ui)));
    expect_true((ui.id==editor.roots[0].id));
    editor.Describe(editor.roots[1].id,ui);
    auto old_sampler=live->lights.volume_scene->light_sampler;
    ui.translation[0]=6;ui.numeric_transform=ui.apply_transform=true;
    ui.materials[0].fields.back().values[0]=5;ui.materials[0].changed=ui.apply_material=true;
    expect_true((editor.Apply(ui)));
    expect_true((live->lights.volume_scene->light_sampler!=old_sampler));
    aabb bound;live->lights.objects[0]->bounding_box(0,1,bound);expect_true((std::abs(bound.Centroid()[0]-6)<1e-5));
    expect_true((ui.materials[0].fields.back().values[0]==5));
  }
  test_that("Shift-click queues object picking without changing the camera target") {
    RayrenderGui gui;gui.can_edit=gui.object.enabled=true;
    rimgui_viewport_v1 view{sizeof(view)};view.width=100;view.height=200;view.mouse_x=25;view.mouse_y=50;
    rimgui_input_v1 input{sizeof(input)};input.mouse_available=1;input.mouse_clicked=RIMGUI_MOUSE_LEFT;input.modifiers=RIMGUI_SHIFT;
    gui.collect_input(input,view);
    expect_true((gui.object.pick_pending && !gui.pick_pending));
    expect_true((gui.object.pick_u==.75f && gui.object.pick_v==.75f));
    gui.object.pick_pending=false;input.mouse_available=0;gui.collect_input(input,view);expect_false(gui.object.pick_pending);
    PreviewObjectState ui;ui.translation={1,2,3};ui.rotation={10,20,30};ui.scale={-2,3,4};
    auto transform=PreviewNumericTransform(ui);PreviewMatrix(transform,ui.model);PreviewDecompose(ui);
    auto again=PreviewNumericTransform(ui);const auto& a=transform.GetMatrix();const auto& b=again.GetMatrix();
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)expect_true((std::abs(a.m[i][j]-b.m[i][j])<1e-5));
  }
  test_that("gizmo projection agrees with perspective and orthographic film rays") {
    auto input=EditorFixture("native_editor_scene(process_scene(sphere())$scene)");PreviewScene editor(input);
    for(int ortho=0;ortho<2;++ortho) {
      Transform object,world;std::unique_ptr<RayCamera> cam;
      if(ortho)cam=std::make_unique<ortho_camera>(point3f(4,3,-10),point3f(0),vec3f(0,1,0),8,6,0,1,1);
      else cam=std::make_unique<camera>(point3f(4,3,-10),point3f(0),vec3f(0,1,0),60,4.f/3,0,10,0,1,1);
#ifdef HAS_OIDN
      PreviewDisplay display(80,60,false,true,false,10,cam.get(),&object,&world,nullptr,nullptr,nullptr,false,false);
#else
      PreviewDisplay display(80,60,false,true,false,10,cam.get(),&object,&world,false);
#endif
      RayrenderGui gui;gui.width=80;gui.height=60;display.AttachNativeGui(&gui,true,true);display.scene_editor=&editor;
      display.UpdateNativeObjectCamera();expect_true(gui.object.projection_valid);
      auto multiply=[](const std::array<float,16>& m,const std::array<double,4>& v) {
        std::array<double,4> r{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)r[i]+=m[i+4*j]*v[j];return r;
      };
      const float x=.2,y=.7;auto ray=cam->get_ray(1-x,1-y,point3f(0),.5);auto p=ray(1);
      auto clip=multiply(gui.object.projection,multiply(gui.object.view,{p[0],p[1],p[2],1}));
      expect_true((std::abs((clip[0]/clip[3]+1)/2-x)<1e-5));
      expect_true((std::abs((1-clip[1]/clip[3])/2-y)<1e-5));
    }
  }
  test_that("material panel defaults remain valid for each supported material class") {
    auto input=EditorFixture(R"(
      mats=list(diffuse(),diffuse(sigma=30),metal(),dielectric(),microfacet(),
        microfacet(transmission=TRUE),glossy(),light(),light(spotlight_focus=c(0,0,0)),hair())
      scene=sphere(material=mats[[1]])
      for(i in 2:length(mats))scene=add_object(scene,sphere(x=i*3,material=mats[[i]]))
      native_editor_scene(process_scene(scene)$scene)
    )");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    for(auto& root:editor.roots)for(auto mat:PreviewMaterials(root.object.get())) {
      auto fields=PreviewMaterialFields(mat);PreviewMaterialEdit edit;edit.type=mat->GetName();
      for(auto& binding:fields)edit.fields.push_back(binding.field);
      expect_true((!edit.fields.empty()));
      bool success=true;try { PreviewApplyMaterial(mat,edit); }catch(...) { success=false; }
      expect_true(success);expect_true((mat->GetName()==edit.type));
    }
  }

  test_that("selected materials and gizmos draw through the provider's headless backend") {
    Rcpp::Function require=Rcpp::Environment::base_env()["requireNamespace"];
    if(!Rcpp::as<bool>(require("rimgui",Rcpp::_["quietly"]=true)))return;
    Rcpp::Function acquire=Rcpp::Environment::namespace_env("rimgui")["acquire_api"];
    Rcpp::RObject handle=acquire();RayrenderGui gui;
    expect_true((rimgui_api_from_R_v1(handle,sizeof(rimgui_api_v1),RIMGUI_CAP_HEADLESS,&gui.api)==0));
    if(gui.api->header.abi_minor<2)return;
    rimgui_session_desc_v1 desc{sizeof(desc),RIMGUI_HEADLESS,"Object editor test",18,1100,760};
    expect_true((gui.api->open(&desc,&gui.session,&gui.error)==0));
    rimgui_callback_v1 callback{sizeof(callback),1,RayrenderGui::draw,&gui};
    expect_true((gui.api->register_callback(gui.session,&callback,&gui.callback,&gui.error)==0));
    gui.width=80;gui.height=60;gui.pixels.assign(80*60*4,128);gui.publish();
    auto input=EditorFixture("native_editor_scene(process_scene(sphere(material=glossy()))$scene)");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    Transform object,world;camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,4.f/3,0,10,0,1,1);
#ifdef HAS_OIDN
    PreviewDisplay display(80,60,false,true,false,10,&cam,&object,&world,nullptr,nullptr,nullptr,false,false);
#else
    PreviewDisplay display(80,60,false,true,false,10,&cam,&object,&world,false);
#endif
    display.AttachNativeGui(&gui,true,true);display.scene_editor=&editor;display.UpdateNativeObjectCamera();
    editor.Describe(editor.roots[0].id,gui.object);
    rimgui_step_v1 step{sizeof(step)};
    for(int i=0;i<3;++i) { gui.object.operation=i;expect_true((gui.api->step(gui.session,&step,&gui.error)==0)); }
    gui.close();
  }
  test_that("rebuilding edited instances retains later alpha textures and volume boundaries") {
    auto input=EditorFixture(R"(
      native_editor_scene(process_scene(add_object(
        create_instances(sphere(),x=c(-4,-2)),
        add_object(sphere(x=2,material=diffuse(alpha_texture=matrix(1,2,2))),
                   set_medium(sphere(x=5),homogeneous_medium(sigma_s=.1)))
      ))$scene)
    )");
    PreviewScene editor(input);auto live=std::make_shared<EditorBuild>(input,editor);
    editor.prepare_rebuild=[&](PreviewScene& next) { auto staged=std::make_shared<EditorBuild>(input,next);return [&,staged]() { live=staged; }; };
    PreviewObjectState ui;editor.Describe(editor.roots[0].id,ui);
    ui.materials[0].fields[0].values={0,1,0};ui.materials[0].changed=ui.apply_material=true;
    expect_true((editor.Apply(ui)));
    expect_true((editor.Pick(Ray(point3f(2,0,-8),vec3f(0,0,1),.5),ui)));
    expect_true((ui.id==editor.roots[2].id));
    auto volume=live->lights.volume_scene;
    expect_true((volume->BoundaryCount()==1));
    editor.Describe(editor.roots[3].id,ui);ui.translation[0]=7;ui.numeric_transform=ui.apply_transform=true;
    expect_true((editor.Apply(ui)));
    aabb bound;live->lights.volume_scene->boundary_bvh->bounding_box(0,1,bound);
    expect_true((std::abs(bound.Centroid()[0]-7)<1e-5));
    auto next_sampler=live->lights.volume_scene->light_sampler;
    volume->ReplacePreparedScene(*live->lights.volume_scene);
    expect_true((volume->light_sampler==next_sampler && volume->BoundaryCount()==1));
    volume->boundary_bvh->bounding_box(0,1,bound);expect_true((std::abs(bound.Centroid()[0]-7)<1e-5));
  }

}
#endif
