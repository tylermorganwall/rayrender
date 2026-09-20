#ifdef NOT_CRAN
#include "preview_scene.h"
#include "preview_surface_maps.h"
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
  Rcpp::Function parse = Rcpp::Environment::base_env()["parse"],
                 eval = Rcpp::Environment::base_env()["eval"];
  Rcpp::Function newenv = Rcpp::Environment::base_env()["new.env"];
  Rcpp::Environment environment =
      newenv(Rcpp::_["parent"] = Rcpp::Environment::namespace_env("rayrender"));
  return Rcpp::as<Rcpp::List>(eval(parse(Rcpp::_["text"] = code), environment));
}
struct EditorBuild {
  TransformCache transforms;
  TextureCache texture_cache;
  std::vector<Float*> textures;
  std::vector<unsigned char*> alpha_textures, bump_textures, roughness_textures;
  std::vector<std::shared_ptr<material>> materials;
  std::vector<std::shared_ptr<alpha_texture>> alpha;
  std::vector<std::shared_ptr<bump_texture>> bump;
  std::vector<std::shared_ptr<roughness_texture>> roughness;
  std::vector<std::shared_ptr<hitable>> instances;
  std::vector<std::shared_ptr<hitable_list>> instance_lights;
  std::vector<int> texture_idx;
  hitable_list lights;
  std::shared_ptr<hitable> root;
  EditorBuild(Rcpp::List input, PreviewScene& editor) {
    Rcpp::IntegerVector shape = input["shape"];
    random_gen rng(1);
    lights.volume_scene = std::make_shared<VolumeScene>();
    root = build_scene(input,
                       shape,
                       0,
                       1,
                       textures,
                       alpha_textures,
                       bump_textures,
                       roughness_textures,
                       &materials,
                       alpha,
                       bump,
                       roughness,
                       0,
                       transforms,
                       texture_cache,
                       lights,
                       instances,
                       instance_lights,
                       texture_idx,
                       false,
                       rng,
                       &editor);
    lights.volume_scene->light_sampler = std::make_shared<VolumeLightSampler>(lights);
  }
};
// Count actual visibility rays to distinguish reusing a mask from recomputing
// an identical mask. Render-color changes should issue no visibility rays.
struct EditorCountingCamera : camera {
  using camera::camera;
  size_t visibility_rays = 0;
  Ray get_ray(Float s, Float t, point3f lens, Float time) override {
    ++visibility_rays;
    return camera::get_ray(s, t, lens, time);
  }
};
uint64_t EditorChild(const PreviewScene& scene, uint64_t parent, size_t index = 0) {
  for (const auto& node : scene.Hierarchy()) {
    if (node.parent == parent) {
      if (index == 0) {
        return node.id;
      }
      --index;
    }
  }
  throw std::runtime_error("Missing editor child in test fixture.");
}
point3f EditorColor(PreviewScene& scene, size_t root) {
  auto fields =
      PreviewMaterialFields(PreviewMaterials(scene.roots[root].object.get())[0]);
  return point3f(
      fields[0].field.values[0], fields[0].field.values[1], fields[0].field.values[2]);
}
}
PreviewField& EditorField(PreviewObjectState& ui, const std::string& name) {
  for (auto& field : ui.materials.at(ui.material_slot).fields) {
    if (field.name == name) {
      return field;
    }
  }
  throw std::runtime_error("Missing material field: " + name);
}
context("Native object editor") {
  test_that("history rebuilds transforms and materials and keeps failed undo atomic") {
    auto input = EditorFixture(
        "native_editor_scene(process_scene(rbind(sphere(material=glossy()), sphere(x=5)))$scene)");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged] {
        live = staged;
      };
    };
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
#ifdef HAS_OIDN
    PreviewDisplay display(4,
                           4,
                           false,
                           true,
                           false,
                           10,
                           &cam,
                           &object,
                           &world,
                           nullptr,
                           nullptr,
                           nullptr,
                           false,
                           false);
#else
    PreviewDisplay display(4, 4, false, true, false, 10, &cam, &object, &world, false);
#endif
    RayrenderGui gui;
    display.AttachNativeGui(&gui, true, false);
    display.scene_editor = &editor;
    editor.Describe(editor.roots[0].id, gui.object);
    const auto original_color = gui.object.materials[0].fields[0].values;
    display.BeginNativeHistory();
    ++gui.history_epoch;
    gui.object.translation[0] = 2;
    gui.object.numeric_transform = gui.object.apply_transform = true;
    expect_true(display.CommitNativeEdits(nullptr));
    expect_true(editor.ExportEdits().size() == 1);
    expect_true(editor.ObjectTransform(0, 0).GetMatrix().m[0][3] == 2);
    gui.history_requests.push_back(-1);
    expect_true(display.CommitNativeEdits(nullptr));
    expect_true(editor.settings.Empty());
    expect_true(gui.object.translation[0] == 0);
    gui.history_requests.push_back(1);
    display.CommitNativeEdits(nullptr);
    expect_true(editor.ObjectTransform(0, 0).GetMatrix().m[0][3] == 2);
    ++gui.history_epoch;
    auto& panel = gui.object.materials[0];
    panel.fields[0].values = {.1, .2, .3};
    panel.fields[0].changed = panel.changed = true;
    gui.object.apply_material = true;
    display.CommitNativeEdits(nullptr);
    const auto edited_color = gui.object.materials[0].fields[0].values;
    expect_true(edited_color != original_color);
    gui.history_requests.push_back(-1);
    display.CommitNativeEdits(nullptr);
    expect_true(gui.object.materials[0].fields[0].values == original_color);
    expect_true(editor.ObjectTransform(0, 0).GetMatrix().m[0][3] == 2);
    gui.history_requests.push_back(1);
    display.CommitNativeEdits(nullptr);
    expect_true(gui.object.materials[0].fields[0].values == edited_color);
    auto build = editor.prepare_rebuild;
    editor.prepare_rebuild = [](PreviewScene&) -> std::function<void()> {
      throw std::runtime_error("Texture file unavailable");
    };
    const auto revision = editor.Revision();
    gui.history_requests.push_back(-1);
    display.CommitNativeEdits(nullptr);
    expect_true(editor.Revision() == revision);
    expect_true(gui.object.materials[0].fields[0].values == edited_color);
    expect_true(gui.history_message == "Undo failed: Texture file unavailable");
    expect_false(gui.can_redo);
    editor.prepare_rebuild = build;
    gui.history_requests.push_back(-1);
    display.CommitNativeEdits(nullptr);
    expect_true(gui.object.materials[0].fields[0].values == original_color);
    expect_true(gui.can_redo);
    // Undo a numeric draft after selecting a different sphere. The next Apply
    // must use the restored sphere's pivot, not the newly selected one's pivot.
    ++gui.history_epoch;
    gui.object.translation[0] = 3;
    gui.object.transform_pending = gui.object.numeric_transform = true;
    gui.history_dirty = true;
    display.CommitNativeEdits(nullptr);
    gui.object.select_id = editor.roots[1].id;
    gui.object.select_pending = true;
    display.CommitNativeEdits(nullptr);
    gui.history_requests.push_back(-1);
    display.CommitNativeEdits(nullptr);
    expect_true(gui.object.id == editor.roots[0].id);
    expect_true(gui.object.translation[0] == 2);
    ++gui.history_epoch;
    gui.object.translation[0] = 4;
    gui.object.numeric_transform = gui.object.apply_transform = true;
    display.CommitNativeEdits(nullptr);
    expect_true(editor.ObjectTransform(0, 0).GetMatrix().m[0][3] == 4);
  }

  test_that("selection outlines scale to film resolution with transparent interiors") {
    std::vector<uint8_t> mask(9 * 9, 0), output;
    for (size_t y = 1; y < 8; ++y) {
      for (size_t x = 1; x < 8; ++x) {
        mask[x + 9 * y] = 1;
      }
    }
    // Neighboring film pixels inside the selected object remain transparent.
    const size_t first = 4 * (8 + 18 * 8), second = first + 4;
    auto outline = PreviewSelectionOverlay::Outline(mask, 9, 9);
    PreviewSelectionOverlay::Compose(18, 18, outline, 9, 9, output);
    expect_true((output.size() == 18 * 18 * 4));
    size_t covered = 0;
    for (size_t i = 0; i < output.size(); i += 4) {
      covered += output[i + 3] == 255;
    }
    expect_true((covered == 160));
    expect_true((output[3] == 0));
    expect_true((output[first + 3] == 0 && output[second + 3] == 0));
    std::fill(mask.begin(), mask.end(), 0);
    outline = PreviewSelectionOverlay::Outline(mask, 9, 9);
    PreviewSelectionOverlay::Compose(18, 18, outline, 9, 9, output);
    expect_true((std::all_of(output.begin(), output.end(), [](uint8_t v) {
      return v == 0;
    })));
  }

  test_that(
      "cached outline fills diagonal boundaries and leaves non-border pixels transparent") {
    using Overlay = PreviewSelectionOverlay;
    const uint32_t width = 13;
    std::vector<uint8_t> mask(width * width, 0), output;
    for (uint32_t y = 0; y < width; ++y) {
      for (uint32_t x = 0; x < width; ++x) {
        // A slanted diamond plus a hole exercises diagonal and concave edges.
        mask[x + width * y] = std::abs(int(x) - 6) + std::abs(int(y) - 6) <= 6;
      }
    }
    mask[6 + width * 6] = 0;
    const auto outline = Overlay::Outline(mask, width, width);
    expect_true((outline[5 + width * 5] == 1)); // Diagonal neighbor of the hole.
    expect_true((outline[6 + width * 6] == 0));
    expect_true((std::count(outline.begin(), outline.end(), 2) > 0));
    Overlay::Compose(width, width, outline, width, width, output);
    bool transparent_interior = true, boundary_opaque = true;
    for (size_t i = 0; i < mask.size(); ++i) {
      if (outline[i]) {
        const uint8_t value = outline[i] == 1 ? 10 : 245;
        boundary_opaque = boundary_opaque && output[4 * i] == value &&
                          output[4 * i + 1] == value && output[4 * i + 2] == value &&
                          output[4 * i + 3] == 255;
      } else {
        transparent_interior = transparent_interior && output[4 * i + 3] == 0;
      }
    }
    expect_true(transparent_interior);
    expect_true(boundary_opaque);
    // One-pixel components and image corners must remain fully covered.
    mask.assign(width * width, 0);
    mask[0] = mask[width + 1] = mask.back() = 1;
    const auto thin = Overlay::Outline(mask, width, width);
    expect_true((thin[0] == 1 && thin[width + 1] == 1 && thin.back() == 1));
  }

  test_that("top-level groups and instance placements pick and rebuild independently") {
    auto input = EditorFixture(R"(
      native_editor_scene(process_scene(add_object(
        group_objects(add_object(sphere(x=-3),sphere(x=-1)), translate=c(0,0,0)),
        create_instances(sphere(material=diffuse(color='red')), x=c(2,5))
      ))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    expect_true((editor.roots.size() == 4));
    expect_true((editor.roots[0].id == editor.roots[1].id &&
                 editor.roots[2].id != editor.roots[3].id));
    PreviewObjectState ui;
    expect_true((editor.Pick(Ray(point3f(-3, 0, -8), vec3f(0, 0, 1), .5), ui)));
    expect_true((ui.id == editor.roots[0].id && ui.materials.size() == 2));
    ui.translation[0] += 1;
    ui.numeric_transform = ui.apply_transform = true;
    expect_true((editor.Apply(ui)));
    expect_true((editor.settings.transforms.size() == 2));
    aabb bound;
    editor.roots[0].object->bounding_box(0, 1, bound);
    expect_true((std::abs(bound.Centroid()[0] + 2) < 1e-5));
    expect_true((editor.Pick(Ray(point3f(2, 0, -8), vec3f(0, 0, 1), .5), ui)));
    expect_true((ui.id == editor.roots[2].id && ui.materials.size() == 1));
    expect_true((PreviewMaterials(editor.roots[2].object.get())[0] ==
                 PreviewMaterials(editor.roots[3].object.get())[0]));
    ui.materials[0].fields[0].values = {0, 1, 0};
    ui.materials[0].changed = true;
    ui.apply_material = true;
    expect_true((editor.Apply(ui)));
    expect_true((PreviewMaterials(editor.roots[2].object.get())[0] !=
                 PreviewMaterials(editor.roots[3].object.get())[0]));
    expect_true(((EditorColor(editor, 2) - point3f(0, 1, 0)).length() < 1e-6));
    expect_true(((EditorColor(editor, 3) - point3f(1, 0, 0)).length() < 1e-6));
    expect_true((editor.ExportEdits().size() == 3));
    ui.revert = true;
    expect_true((editor.Apply(ui)));
    expect_true(((EditorColor(editor, 2) - point3f(1, 0, 0)).length() < 1e-6));
    expect_true((editor.ExportEdits().size() == 2));
    expect_false(editor.Pick(Ray(point3f(100, 0, -8), vec3f(0, 0, 1), .5), ui));
    expect_false(ui.selected);
  }
  test_that("invalid transforms and failed rebuilds leave the live scene intact") {
    auto input = EditorFixture("native_editor_scene(process_scene(sphere())$scene)");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    bool published = false;
    editor.prepare_rebuild = [&](PreviewScene& next) -> std::function<void()> {
      throw std::runtime_error("rebuild failed");
    };
    auto original = editor.roots[0].object;
    ui.scale[0] = 0;
    ui.numeric_transform = ui.apply_transform = true;
    expect_false(editor.Apply(ui));
    expect_true((!ui.error.empty() && editor.settings.transforms.empty()));
    ui.translation[0] = 7;
    ui.numeric_transform = ui.apply_transform = true;
    expect_false(editor.Apply(ui));
    expect_true((ui.error == "rebuild failed" && editor.roots[0].object == original));
    ui.model[0] = std::numeric_limits<float>::quiet_NaN();
    ui.apply_transform = true;
    expect_false(editor.Apply(ui));
    expect_true((editor.roots[0].object == original));
    ui.model[12] = 3;
    ui.transform_pending = ui.cancel_transform = true;
    editor.CancelTransform(ui);
    expect_true((ui.model[12] == 0 && !ui.transform_pending && !ui.cancel_transform));
    expect_false(published);
  }
  test_that("mesh selection keeps its root and light edits rebuild light sampling") {
    auto input = EditorFixture(R"(
      path = tempfile(fileext='.obj')
      writeLines(c('v -1 -1 0','v 1 -1 0','v 1 1 0','v -1 1 0','f 1 2 3','f 1 3 4'),path)
      native_editor_scene(process_scene(add_object(
        obj_model(path, material=diffuse(color='blue'),load_material=FALSE),
        sphere(x=4,material=light(color='red',intensity=2))
      ))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    expect_true((editor.roots.size() == 2));
    PreviewObjectState ui;
    expect_true((editor.Pick(Ray(point3f(0, 0, -10), vec3f(0, 0, 1), .5), ui)));
    expect_true((ui.id == editor.roots[0].id));
    editor.Describe(editor.roots[1].id, ui);
    auto old_sampler = live->lights.volume_scene->light_sampler;
    ui.translation[0] = 6;
    ui.numeric_transform = ui.apply_transform = true;
    auto light_field = std::find_if(ui.materials[0].fields.begin(),
                                    ui.materials[0].fields.end(),
                                    [](const PreviewField& field) {
                                      return field.name == "Light intensity";
                                    });
    const size_t intensity_index = size_t(light_field - ui.materials[0].fields.begin());
    ui.materials[0].fields[intensity_index].values[0] = 5;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true((editor.Apply(ui)));
    expect_true((live->lights.volume_scene->light_sampler != old_sampler));
    aabb bound;
    live->lights.objects[0]->bounding_box(0, 1, bound);
    expect_true((std::abs(bound.Centroid()[0] - 6) < 1e-5));
    expect_true((ui.materials[0].fields[intensity_index].values[0] == 5));
  }
  test_that("Shift-click queues object picking without changing the camera target") {
    RayrenderGui gui;
    gui.can_edit = gui.object.enabled = true;
    rayimgui_viewport_v1 view{sizeof(view)};
    view.width = 100;
    view.height = 200;
    view.mouse_x = 25;
    view.mouse_y = 50;
    rayimgui_input_v1 input{sizeof(input)};
    input.mouse_available = 1;
    input.mouse_clicked = RAYIMGUI_MOUSE_LEFT;
    input.modifiers = RAYIMGUI_SHIFT;
    gui.collect_input(input, view);
    expect_true((gui.object.pick_pending && !gui.pick_pending));
    expect_true((gui.object.pick_u == .75f && gui.object.pick_v == .75f));
    gui.object.pick_pending = false;
    input.modifiers = 0;
    gui.collect_input(input, view);
    expect_false(gui.object.pick_pending);
    expect_true((gui.pick_pending && gui.pick_focus));
    gui.pick_pending = false;
    input.modifiers = RAYIMGUI_SHIFT | RAYIMGUI_ALT;
    gui.collect_input(input, view);
    expect_false(gui.object.pick_pending);
    expect_true((gui.pick_pending && gui.pick_focus));
    gui.pick_pending = false;
    input.modifiers = RAYIMGUI_SHIFT;
    input.mouse_clicked = RAYIMGUI_MOUSE_RIGHT;
    gui.collect_input(input, view);
    expect_false(gui.object.pick_pending);
    expect_true((gui.pick_pending && !gui.pick_focus));
    gui.pick_pending = false;
    input.mouse_clicked = RAYIMGUI_MOUSE_LEFT;
    input.mouse_available = 0;
    gui.collect_input(input, view);
    expect_false(gui.object.pick_pending);
    PreviewObjectState ui;
    ui.translation = {1, 2, 3};
    ui.rotation = {10, 20, 30};
    ui.scale = {-2, 3, 4};
    auto transform = PreviewNumericTransform(ui);
    PreviewMatrix(transform, ui.model);
    PreviewDecompose(ui);
    auto again = PreviewNumericTransform(ui);
    const auto& a = transform.GetMatrix();
    const auto& b = again.GetMatrix();
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        expect_true((std::abs(a.m[i][j] - b.m[i][j]) < 1e-5));
      }
    }
  }
  test_that("gizmo projection agrees with perspective and orthographic film rays") {
    auto input = EditorFixture("native_editor_scene(process_scene(sphere())$scene)");
    PreviewScene editor(input);
    for (int ortho = 0; ortho < 2; ++ortho) {
      Transform object, world;
      std::unique_ptr<RayCamera> cam;
      if (ortho) {
        cam = std::make_unique<ortho_camera>(
            point3f(4, 3, -10), point3f(0), vec3f(0, 1, 0), 8, 6, 0, 1, 1);
      } else {
        cam = std::make_unique<camera>(point3f(4, 3, -10),
                                       point3f(0),
                                       vec3f(0, 1, 0),
                                       60,
                                       4.f / 3,
                                       0,
                                       10,
                                       0,
                                       1,
                                       1);
      }
#ifdef HAS_OIDN
      PreviewDisplay display(80,
                             60,
                             false,
                             true,
                             false,
                             10,
                             cam.get(),
                             &object,
                             &world,
                             nullptr,
                             nullptr,
                             nullptr,
                             false,
                             false);
#else
      PreviewDisplay display(
          80, 60, false, true, false, 10, cam.get(), &object, &world, false);
#endif
      RayrenderGui gui;
      gui.width = 80;
      gui.height = 60;
      display.AttachNativeGui(&gui, true, true);
      display.scene_editor = &editor;
      display.UpdateNativeObjectCamera();
      expect_true(gui.object.projection_valid);
      auto multiply = [](const std::array<float, 16>& m,
                         const std::array<double, 4>& v) {
        std::array<double, 4> r{};
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            r[i] += m[i + 4 * j] * v[j];
          }
        }
        return r;
      };
      const float x = .2, y = .7;
      auto ray = cam->get_ray(1 - x, 1 - y, point3f(0), .5);
      auto p = ray(1);
      auto clip = multiply(gui.object.projection,
                           multiply(gui.object.view, {p[0], p[1], p[2], 1}));
      expect_true((std::abs((clip[0] / clip[3] + 1) / 2 - x) < 1e-5));
      expect_true((std::abs((1 - clip[1] / clip[3]) / 2 - y) < 1e-5));
    }
  }
  test_that("material panel defaults remain valid for each supported material class") {
    auto input = EditorFixture(R"(
      mats=list(diffuse(),diffuse(sigma=30),metal(),dielectric(),microfacet(),
        microfacet(transmission=TRUE),glossy(),light(),light(spotlight_focus=c(0,0,0)),hair())
      scene=sphere(material=mats[[1]])
      for(i in 2:length(mats))scene=add_object(scene,sphere(x=i*3,material=mats[[i]]))
      native_editor_scene(process_scene(scene)$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    for (auto& root : editor.roots) {
      for (auto mat : PreviewMaterials(root.object.get())) {
        auto fields = PreviewMaterialFields(mat);
        PreviewMaterialEdit edit;
        edit.type = mat->GetName();
        for (auto& binding : fields) {
          edit.fields.push_back(binding.field);
        }
        expect_true((!edit.fields.empty()));
        bool success = true;
        try {
          PreviewApplyMaterial(mat, edit);
        } catch (...) {
          success = false;
        }
        expect_true(success);
        expect_true((mat->GetName() == edit.type));
      }
    }
  }

  test_that(
      "selected materials and gizmos draw through the provider's headless backend") {
    Rcpp::Function require = Rcpp::Environment::base_env()["requireNamespace"];
    if (!Rcpp::as<bool>(require("rayimgui", Rcpp::_["quietly"] = true))) {
      return;
    }
    Rcpp::Function acquire =
        Rcpp::Environment::namespace_env("rayimgui")["acquire_api"];
    Rcpp::RObject handle = acquire();
    RayrenderGui gui;
    expect_true(
        (rayimgui_api_from_R_v1(
             handle, sizeof(rayimgui_api_v1), RAYIMGUI_CAP_HEADLESS, &gui.api) == 0));
    if (gui.api->header.abi_minor < 8) {
      return;
    }
    rayimgui_session_desc_v1 desc{
        sizeof(desc), RAYIMGUI_HEADLESS, "Object editor test", 18, 1100, 760};
    expect_true((gui.api->open(&desc, &gui.session, &gui.error) == 0));
    rayimgui_callback_v1 callback{sizeof(callback), 1, RayrenderGui::draw, &gui};
    expect_true((gui.api->register_callback(
                     gui.session, &callback, &gui.callback, &gui.error) == 0));
    gui.width = 80;
    gui.height = 60;
    gui.pixels.assign(80 * 60 * 4, 128);
    gui.publish();
    auto input = EditorFixture(
        "native_editor_scene(process_scene(sphere(material=glossy()))$scene)");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    Transform object, world;
    EditorCountingCamera cam(
        point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 4.f / 3, 0, 10, 0, 1, 1);
#ifdef HAS_OIDN
    PreviewDisplay display(80,
                           60,
                           false,
                           true,
                           false,
                           10,
                           &cam,
                           &object,
                           &world,
                           nullptr,
                           nullptr,
                           nullptr,
                           false,
                           false);
#else
    PreviewDisplay display(
        80, 60, false, true, false, 10, &cam, &object, &world, false);
#endif
    display.AttachNativeGui(&gui, true, true);
    display.scene_editor = &editor;
    display.UpdateNativeObjectCamera();
    editor.Describe(editor.roots[0].id, gui.object);
    const auto image_pixels = gui.pixels;
    expect_true(display.UpdateNativeSelectionMask());
    const size_t first_trace = cam.visibility_rays;
    expect_true((first_trace == 80 * 60));
    expect_false(display.UpdateNativeSelectionMask());
    expect_true((cam.visibility_rays == first_trace));
    expect_true((gui.selection_id == gui.object.id));
    expect_true(
        (std::count(gui.selection_mask.begin(), gui.selection_mask.end(), 1) > 0));
    expect_true((gui.pixels == image_pixels));
    expect_true((gui.selection_pixels.size() == 80 * 60 * 4));
    const auto original_mask = gui.selection_mask;
    const auto original_outline = gui.selection_outline;
    const auto previous_overlay = gui.selection_pixels;
    for (size_t i = 0; i < gui.pixels.size(); i += 4) {
      gui.pixels[i] = gui.pixels[i + 1] = 255;
      gui.pixels[i + 2] = 0;
    }
    const auto yellow_image = gui.pixels;
    gui.publish();
    expect_true((gui.selection_pixels == previous_overlay));
    expect_true((gui.selection_mask == original_mask));
    expect_true((gui.pixels == yellow_image));
    expect_false(display.UpdateNativeSelectionMask());
    expect_true((cam.visibility_rays == first_trace));
    expect_true((gui.selection_outline == original_outline));
    // Fast mode restarts accumulation, but camera geometry and selection did
    // not change. Restoring a hidden overlay should reuse its cached coverage.
    gui.fast_preview = !display.write_fast_output;
    gui.fast_pending = true;
    expect_true(display.ApplyNativeControls(nullptr));
    gui.selection_visible = false;
    expect_true(display.UpdateNativeSelectionMask());
    expect_true(gui.selection_visible);
    expect_true((cam.visibility_rays == first_trace));
    expect_true((gui.selection_mask == original_mask &&
                 gui.selection_outline == original_outline));
    cam.update_position(vec3f(2, 0, 0), false);
    expect_true(display.UpdateNativeSelectionMask());
    expect_true((gui.selection_mask != original_mask));
    expect_true((cam.visibility_rays == 2 * first_trace));
    expect_true((gui.selection_outline != original_outline));
    // The image sky panel must draw even without native atmosphere controls.
    gui.has_sky_model = gui.has_sun = gui.has_location = true;
    gui.set_datetime("2026-06-21 16:00:00");
    // Exercise text, choice, boolean and numeric material widgets through the provider.
    EditorField(gui.object, "Texture mode").values[0] = 5;
    EditorField(gui.object, "Use roughness map").values[0] = 1;
    EditorField(gui.object, "Use alpha map").values[0] = 1;
    EditorField(gui.object, "Use bump map").values[0] = 1;
    rayimgui_step_v1 step{sizeof(step)};
    for (int i = 0; i < 3; ++i) {
      gui.object.operation = i;
      expect_true((gui.api->step(gui.session, &step, &gui.error) == 0));
    }
    editor.Describe(0, gui.object);
    expect_true(display.UpdateNativeSelectionMask());
    expect_true((gui.selection_id == 0 && !gui.selection_visible));
    expect_true((cam.visibility_rays == 2 * first_trace));
    gui.close();
  }
  test_that(
      "rebuilding edited instances retains later alpha textures and volume boundaries") {
    auto input = EditorFixture(R"(
      native_editor_scene(process_scene(add_object(
        create_instances(sphere(),x=c(-4,-2)),
        add_object(sphere(x=2,material=diffuse(alpha_texture=matrix(1,2,2))),
                   set_medium(sphere(x=5),homogeneous_medium(sigma_s=.1)))
      ))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    ui.materials[0].fields[0].values = {0, 1, 0};
    ui.materials[0].changed = ui.apply_material = true;
    expect_true((editor.Apply(ui)));
    expect_true((editor.Pick(Ray(point3f(2, 0, -8), vec3f(0, 0, 1), .5), ui)));
    expect_true((ui.id == editor.roots[2].id));
    auto volume = live->lights.volume_scene;
    expect_true((volume->BoundaryCount() == 1));
    editor.Describe(editor.roots[3].id, ui);
    ui.translation[0] = 7;
    ui.numeric_transform = ui.apply_transform = true;
    expect_true((editor.Apply(ui)));
    aabb bound;
    live->lights.volume_scene->boundary_bvh->bounding_box(0, 1, bound);
    expect_true((std::abs(bound.Centroid()[0] - 7) < 1e-5));
    auto next_sampler = live->lights.volume_scene->light_sampler;
    volume->ReplacePreparedScene(*live->lights.volume_scene);
    expect_true(
        (volume->light_sampler == next_sampler && volume->BoundaryCount() == 1));
    volume->boundary_bvh->bounding_box(0, 1, bound);
    expect_true((std::abs(bound.Centroid()[0] - 7) < 1e-5));
  }
  test_that(
      "hierarchy selection can edit a group child without changing its siblings") {
    auto input = EditorFixture(R"(
      inner = group_objects(rbind(sphere(x=-3), sphere(x=-1)))
      native_editor_scene(process_scene(group_objects(rbind(inner, cube(x=2))))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    auto tree = editor.Hierarchy();
    expect_true((tree.size() == 5));
    expect_true((tree[0].parent == 0 && tree[1].parent == tree[0].id));
    expect_true((tree[2].parent == tree[1].id && tree[3].parent == tree[1].id));
    expect_true((tree[4].parent == tree[0].id));
    PreviewObjectState ui;
    editor.Describe(tree[2].id, ui);
    expect_true((ui.selected && ui.materials.size() == 1));
    ui.translation[0] += 1;
    ui.numeric_transform = ui.apply_transform = true;
    expect_true(editor.Apply(ui));
    expect_true((editor.settings.transforms.size() == 1));
    editor.Describe(tree[1].id, ui);
    expect_true((ui.materials.size() == 2));
    ui.translation[1] += 2;
    ui.numeric_transform = ui.apply_transform = true;
    expect_true(editor.Apply(ui));
    expect_true((editor.settings.transforms.size() == 2));
    aabb bounds;
    editor.roots[2].object->bounding_box(0, 1, bounds);
    expect_true((std::abs(bounds.Centroid()[1]) < 1e-5));
    editor.Describe(tree[0].id, ui);
    expect_true((ui.materials.size() == 3));
    editor.Describe(UINT64_C(123), ui);
    expect_false(ui.selected);
  }

  test_that(
      "Shift-click reaches selection through an idle gizmo but respects active drags") {
    struct InputState {
      uint32_t modifiers = RAYIMGUI_SHIFT;
      unsigned gizmos = 0;
      bool captured = false, active = false;
    } state;
    const auto owner = uint64_t(reinterpret_cast<uintptr_t>(&state));
    rayimgui_api_v1 api{};
    api.widget = [](uint64_t,
                    rayimgui_widget_v1*,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      item->flags = RAYIMGUI_VISIBLE;
      return 0;
    };
    api.viewport =
        [](uint64_t, rayimgui_viewport_v1* view, rayimgui_error_v1*) -> int32_t {
      view->width = view->height = 100;
      view->mouse_x = view->mouse_y = 50;
      return 0;
    };
    api.input =
        [](uint64_t owner, rayimgui_input_v1* input, rayimgui_error_v1*) -> int32_t {
      const auto& state = *reinterpret_cast<InputState*>(uintptr_t(owner));
      input->modifiers = state.modifiers;
      input->mouse_available = !state.captured;
      input->mouse_clicked = state.captured ? 0 : RAYIMGUI_MOUSE_LEFT;
      return 0;
    };
    api.gizmo = [](uint64_t owner,
                   const rayimgui_gizmo_v1*,
                   rayimgui_item_v1* item,
                   rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<InputState*>(uintptr_t(owner));
      ++state.gizmos;
      state.captured = true;
      item->flags = state.active ? RAYIMGUI_ACTIVE : RAYIMGUI_HOVERED;
      return 0;
    };
    RayrenderGui gui;
    gui.can_edit = gui.object.enabled = gui.object.selected =
        gui.object.projection_valid = true;
    gui.texture = 1;
    expect_true((gui.draw_impl(&api, owner, &gui.error) == 0));
    expect_true((state.gizmos == 0 && gui.object.pick_pending && !gui.pick_pending));
    gui.object.pick_pending = false;
    gui.object.transform_active = state.active = true;
    state.captured = false;
    expect_true((gui.draw_impl(&api, owner, &gui.error) == 0));
    expect_true(
        (state.gizmos == 1 && !gui.object.pick_pending && gui.object.transform_active));
    state.modifiers = 0;
    state.captured = state.active = gui.object.transform_active = false;
    expect_true((gui.draw_impl(&api, owner, &gui.error) == 0));
    expect_true((state.gizmos == 2 && !gui.object.pick_pending));
  }

  test_that("instance-child masks include foreground siblings in visibility") {
    auto input = EditorFixture(R"(
      source = rbind(sphere(material=diffuse(color='red')),sphere(z=-3))
      native_editor_scene(process_scene(create_instances(source))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    const auto parent = editor.roots[0].id;
    const auto rear = EditorChild(editor, parent);
    const auto front = EditorChild(editor, parent, 1);
    auto ray = [](float, float, Ray& out) {
      out = Ray(point3f(0, 0, -8), vec3f(0, 0, 1), .5);
      return true;
    };
    expect_true((editor.SelectionMask(rear, 1, 1, ray) == std::vector<uint8_t>{0}));
    PreviewObjectState ui;
    editor.Describe(front, ui);
    ui.translation[0] += 4;
    ui.numeric_transform = ui.apply_transform = true;
    expect_true(editor.Apply(ui));
    expect_true((editor.SelectionMask(rear, 1, 1, ray) == std::vector<uint8_t>{1}));
  }

  test_that("repeated picking drills through groups and stays at the chosen depth") {
    auto input = EditorFixture(R"(
      inner = group_objects(rbind(sphere(x=-3),sphere(x=-1)))
      native_editor_scene(process_scene(rbind(
        group_objects(rbind(inner,sphere(x=2))),sphere(x=6)))$scene)
    )");
    PreviewScene editor(input);
    EditorBuild live(input, editor);
    const auto group = EditorChild(editor, 0);
    const auto inner = EditorChild(editor, group);
    const auto first = EditorChild(editor, inner);
    const auto second = EditorChild(editor, inner, 1);
    PreviewObjectState ui;
    const Ray first_ray(point3f(-3, 0, -8), vec3f(0, 0, 1), .5);
    expect_true(editor.Pick(first_ray, ui));
    expect_true((ui.id == group));
    expect_true(editor.Pick(first_ray, ui));
    expect_true((ui.id == inner));
    expect_true(editor.Pick(first_ray, ui));
    expect_true((ui.id == first));
    ui.transform_pending = true;
    expect_true(editor.Pick(first_ray, ui));
    expect_true((ui.id == first && ui.transform_pending));
    expect_true(editor.Pick(Ray(point3f(-1, 0, -8), vec3f(0, 0, 1), .5), ui));
    expect_true((ui.id == second));
    expect_true(editor.Pick(Ray(point3f(6, 0, -8), vec3f(0, 0, 1), .5), ui));
    expect_true((ui.id == editor.roots[3].id));
    expect_false(editor.Pick(first_ray, ui, []() {
      return true;
    }));
    expect_true((ui.id == editor.roots[3].id));
    expect_false(editor.Pick(Ray(point3f(100, 0, -8), vec3f(0, 0, 1), .5), ui));
    expect_false(ui.selected);
  }

  test_that("instance child selection edits only that copy in world coordinates") {
    auto input = EditorFixture(R"(
      source = rbind(sphere(x=-1,radius=.4,material=diffuse(color='red')),
                     cube(x=1,width=.8,material=diffuse(color='blue')))
      native_editor_scene(process_scene(create_instances(source,x=c(-4,4),
        angle_z=c(90,0),scale_x=2,scale_y=2,scale_z=2))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    const auto parent = editor.roots[0].id;
    const auto child = EditorChild(editor, parent);
    const auto tree = editor.Hierarchy();
    PreviewObjectState ui;
    const Ray ray(point3f(-4, -2, -8), vec3f(0, 0, 1), .5);
    expect_true(editor.Pick(ray, ui));
    expect_true((ui.id == parent));
    auto source_material =
        PreviewMaterials(editor.roots[0].contents->roots[0].object.get())[0];
    const auto slots = PreviewMaterials(editor.roots[0].object.get());
    const size_t slot =
        std::find(slots.begin(), slots.end(), source_material) - slots.begin();
    ui.materials[slot].fields[0].values = {1, 1, 0};
    ui.materials[slot].fields[0].changed = ui.materials[slot].changed =
        ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true(editor.Pick(ray, ui));
    expect_true((ui.id == child));
    expect_true((std::abs(ui.translation[0] + 4) < 1e-5 &&
                 std::abs(ui.translation[1] + 2) < 1e-5));
    ui.materials[0].fields[0].values = {0, 1, 0};
    ui.materials[0].fields[0].changed = ui.materials[0].changed = ui.apply_material =
        true;
    ui.translation[0] += 1;
    ui.translation[1] += 6; // Cross the sibling so the child BVH changes order.
    ui.numeric_transform = ui.apply_transform = true;
    expect_true(editor.Apply(ui));
    expect_true((ui.id == child && std::abs(ui.translation[0] + 3) < 1e-5));
    expect_true(
        ((EditorColor(*editor.roots[0].contents, 0) - point3f(0, 1, 0)).length() <
         1e-6));
    expect_true(
        ((EditorColor(*editor.roots[1].contents, 0) - point3f(1, 0, 0)).length() <
         1e-6));
    expect_true(
        ((EditorColor(*editor.roots[0].contents, 1) - point3f(0, 0, 1)).length() <
         1e-6));
    const auto placed =
        editor.roots[0].frame * editor.roots[0].contents->roots[0].frame;
    expect_true(((placed(point3f(0)) - point3f(-3, 4, 0)).length() < 1e-5));
    auto new_tree = editor.Hierarchy();
    expect_true((tree.size() == new_tree.size()));
    bool stable = true;
    for (size_t i = 0; i < tree.size(); ++i) {
      stable = stable && tree[i].id == new_tree[i].id &&
               tree[i].parent == new_tree[i].parent;
    }
    expect_true(stable);
    const auto rays = [](float u, float, Ray& out) {
      const point3f points[] = {
          point3f(-3, 4, -8), point3f(-4, 2, -8), point3f(2, 0, -8), point3f(0, 8, -8)};
      out = Ray(points[std::min(3, int(u * 4))], vec3f(0, 0, 1), .5);
      return true;
    };
    expect_true(
        (editor.SelectionMask(child, 4, 1, rays) == std::vector<uint8_t>{1, 0, 0, 0}));
    expect_true(
        (editor.SelectionMask(parent, 4, 1, rays) == std::vector<uint8_t>{1, 1, 0, 0}));
    Rcpp::List edits = editor.ExportEdits();
    Rcpp::List parent_edit = edits[0];
    Rcpp::List children = parent_edit["children"];
    expect_true((edits.size() == 1 && children.size() == 1));
    // Rebuild from only the original R scene and the exported value records.
    PreviewScene restored(input);
    restored.ImportEdits(input, edits);
    EditorBuild replay(input, restored);
    expect_true(
        ((EditorColor(*restored.roots[0].contents, 0) - point3f(0, 1, 0)).length() <
         1e-6));
    expect_true(
        ((EditorColor(*restored.roots[1].contents, 0) - point3f(1, 0, 0)).length() <
         1e-6));
    const auto replayed =
        restored.roots[0].frame * restored.roots[0].contents->roots[0].frame;
    expect_true(((replayed(point3f(0)) - point3f(-3, 4, 0)).length() < 1e-5));
    // A subsequent parent edit also updates that child's narrower override.
    editor.Describe(parent, ui);
    ui.materials[slot].fields[0].values = {.2, .3, .4};
    ui.materials[slot].fields[0].changed = ui.materials[slot].changed =
        ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true(
        ((EditorColor(*editor.roots[0].contents, 0) - point3f(.2, .3, .4)).length() <
         1e-6));
    editor.Describe(child, ui);
    ui.revert = true;
    expect_true(editor.Apply(ui));
    expect_true((std::abs(ui.translation[0] + 4) < 1e-5));
    expect_true(
        ((EditorColor(*editor.roots[0].contents, 0) - point3f(.2, .3, .4)).length() <
         1e-6));
    expect_true(
        ((EditorColor(*editor.roots[1].contents, 0) - point3f(1, 0, 0)).length() <
         1e-6));
  }

  test_that("nested instances retain stable child identity and isolate descendants") {
    auto input = EditorFixture(R"(
      source = create_instances(sphere(radius=.3,material=diffuse(color='red')),x=c(-1,1))
      native_editor_scene(process_scene(create_instances(source,x=c(-4,4)))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    const auto outer = editor.roots[0].id;
    const auto collection = EditorChild(editor, outer);
    const auto placement = EditorChild(editor, collection);
    const auto sphere = EditorChild(editor, placement);
    PreviewObjectState ui;
    const Ray ray(point3f(-5, 0, -8), vec3f(0, 0, 1), .5);
    for (auto expected : {outer, collection, placement, sphere}) {
      expect_true(editor.Pick(ray, ui));
      expect_true((ui.id == expected));
    }
    ui.materials[0].fields[0].values = {0, 1, 0};
    ui.materials[0].fields[0].changed = ui.materials[0].changed = ui.apply_material =
        true;
    expect_true(editor.Apply(ui));
    expect_true((ui.id == sphere));
    auto& first = *editor.roots[0].contents;
    auto& second = *editor.roots[1].contents;
    expect_true((
        (EditorColor(*first.roots[0].contents, 0) - point3f(0, 1, 0)).length() < 1e-6));
    expect_true((
        (EditorColor(*first.roots[1].contents, 0) - point3f(1, 0, 0)).length() < 1e-6));
    expect_true(
        ((EditorColor(*second.roots[0].contents, 0) - point3f(1, 0, 0)).length() <
         1e-6));
    const auto rays = [](float u, float, Ray& out) {
      out = Ray(point3f(-5 + 2 * int(u * 2), 0, -8), vec3f(0, 0, 1), .5);
      return true;
    };
    expect_true(
        (editor.SelectionMask(sphere, 2, 1, rays) == std::vector<uint8_t>{1, 0}));
    const auto previous = editor.roots[0].object;
    editor.prepare_rebuild = [](PreviewScene&) -> std::function<void()> {
      throw std::runtime_error("deliberate nested failure");
    };
    ui.translation[0] += 1;
    ui.numeric_transform = ui.apply_transform = true;
    expect_false(editor.Apply(ui));
    expect_true((editor.roots[0].object == previous && ui.id == sphere));
  }

  test_that("instance hierarchy exposes stable placement selection") {
    auto input = EditorFixture(
        "native_editor_scene(process_scene(create_instances(sphere(), x=c(2,5)))$scene)");
    PreviewScene editor(input);
    EditorBuild live(input, editor);
    auto tree = editor.Hierarchy();
    expect_true((tree.size() == 5 && tree[0].parent == 0));
    const auto second = EditorChild(editor, tree[0].id, 1);
    expect_true((EditorChild(editor, second) != second));
    PreviewObjectState ui;
    editor.Describe(second, ui);
    expect_true((ui.id == editor.roots[1].id && ui.materials.size() == 1));
    editor.Describe(tree[0].id, ui);
    expect_true((ui.materials.size() == 1 && ui.materials[0].targets.size() == 2));
  }

  test_that(
      "instance parent material edits preserve other slots and untouched properties") {
    auto input = EditorFixture(R"(
      source = rbind(sphere(x=-.5, radius=.3, material=metal(color='red', fuzz=.1)),
                     sphere(x=.5, radius=.3, material=diffuse(color='blue')))
      native_editor_scene(process_scene(create_instances(source, x=c(-3,3)))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    const auto tree = editor.Hierarchy();
    PreviewObjectState ui;
    editor.Describe(EditorChild(editor, tree[0].id, 1), ui);
    size_t slot = 0, fuzz = 0;
    for (size_t i = 0; i < ui.materials.size(); ++i) {
      for (size_t j = 0; j < ui.materials[i].fields.size(); ++j) {
        if (ui.materials[i].fields[j].name == "Fuzz") {
          slot = i;
          fuzz = j;
        }
      }
    }
    ui.materials[slot].fields[fuzz].values[0] = .7;
    ui.materials[slot].fields[fuzz].changed = true;
    ui.materials[slot].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((editor.settings.materials.size() == 1));
    editor.Describe(tree[0].id, ui);
    expect_true((ui.materials.size() == 2));
    expect_true((ui.materials[slot].targets.size() == 2));
    expect_true(ui.materials[slot].fields[fuzz].mixed);
    const size_t source_slot = ui.materials[slot].slot;
    ui.materials[slot].fields[0].values = {0, 1, 0};
    ui.materials[slot].fields[0].changed = true;
    ui.materials[slot].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((editor.settings.materials.size() == 2));
    for (size_t i = 0; i < editor.roots.size(); ++i) {
      auto materials = PreviewMaterials(editor.roots[i].object.get());
      auto fields = PreviewMaterialFields(materials[source_slot]);
      expect_true((std::abs(fields[fuzz].field.values[0] - (i == 0 ? .1 : .7)) < 1e-6));
      expect_true((fields[0].field.values[0] == 0 && fields[0].field.values[1] == 1));
      auto other = PreviewMaterialFields(materials[1 - source_slot]);
      expect_true((other[0].field.values[2] == 1 && other[0].field.values[1] == 0));
    }
    expect_true((editor.ExportEdits().size() == 2));
    // A failed batch must leave every live placement unchanged.
    editor.Describe(tree[0].id, ui);
    ui.materials[slot].fields[0].values = {1, 0, 0};
    ui.materials[slot].fields[0].changed = true;
    ui.materials[slot].changed = ui.apply_material = true;
    editor.prepare_rebuild = [](PreviewScene&) -> std::function<void()> {
      throw std::runtime_error("deliberate batch failure");
    };
    expect_false(editor.Apply(ui));
    for (const auto& root : editor.roots) {
      auto fields =
          PreviewMaterialFields(PreviewMaterials(root.object.get())[source_slot]);
      expect_true((fields[0].field.values[1] == 1));
    }
  }

  test_that(
      "selection masks respect instance membership, occlusion and rebuilt geometry") {
    auto input = EditorFixture(R"(
      native_editor_scene(process_scene(rbind(
        create_instances(sphere(), x=c(-3,3)),
        sphere(x=3,z=-3)))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    const auto tree = editor.Hierarchy();
    auto rays = [](float u, float v, Ray& ray) {
      ray = Ray(point3f((u - .5f) * 9, 0, -8), vec3f(0, 0, 1), .5);
      return true;
    };
    auto mask = editor.SelectionMask(tree[0].id, 3, 1, rays);
    expect_true((mask == std::vector<uint8_t>{1, 0, 0}));
    expect_true((editor.SelectionMask(editor.roots[1].id, 3, 1, rays) ==
                 std::vector<uint8_t>{0, 0, 0}));
    expect_true((editor.SelectionMask(editor.roots[2].id, 3, 1, rays) ==
                 std::vector<uint8_t>{0, 0, 1}));
    expect_true(editor
                    .SelectionMask(tree[0].id,
                                   3,
                                   1,
                                   rays,
                                   []() {
                                     return true;
                                   })
                    .empty());
    PreviewObjectState ui;
    editor.Describe(editor.roots[2].id, ui);
    ui.translation[0] += 10;
    ui.numeric_transform = ui.apply_transform = true;
    const auto revision = editor.Revision();
    expect_true(editor.Apply(ui));
    expect_true((editor.Revision() != revision));
    expect_true((editor.SelectionMask(tree[0].id, 3, 1, rays) ==
                 std::vector<uint8_t>{1, 0, 1}));
    expect_true((editor.SelectionMask(editor.roots[1].id, 3, 1, rays) ==
                 std::vector<uint8_t>{0, 0, 1}));
    expect_true((editor.SelectionMask(0, 3, 1, rays) == std::vector<uint8_t>{0, 0, 0}));
  }
  test_that("each material exposes and commits its shading inputs") {
    auto input = EditorFixture(R"(
      mats=list(diffuse(),diffuse(sigma=30),metal(),dielectric(),microfacet(roughness=.2),
        microfacet(transmission=TRUE,roughness=.2,eta=1.5),glossy(),light(),
        light(spotlight_focus=c(0,0,0)),hair())
      scene=sphere(material=mats[[1]])
      for(i in 2:length(mats))scene=add_object(scene,sphere(x=i*3,material=mats[[i]]))
      native_editor_scene(process_scene(scene)$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged] {
        live = staged;
      };
    };
    const std::vector<std::pair<std::string, double>> edits = {
        {"Sigma (degrees)", 40},
        {"Sigma (degrees)", 60},
        {"Eta RGB", 2},
        {"Priority", 3},
        {"Microfacet distribution", 1},
        {"Index of refraction", 1.8},
        {"Gloss X/Y", .4},
        {"Invisible", 1},
        {"Spotlight width (degrees)", 45},
        {"Scale angle (degrees)", 8}};
    for (size_t row = 0; row < edits.size(); ++row) {
      PreviewObjectState ui;
      editor.Describe(editor.roots[row].id, ui);
      auto& field = EditorField(ui, edits[row].first);
      field.values[0] = edits[row].second;
      field.changed = ui.materials[0].changed = ui.apply_material = true;
      expect_true(editor.Apply(ui));
      expect_true((std::abs(EditorField(ui, edits[row].first).values[0] -
                            edits[row].second) < .001));
      expect_true(ui.error.empty());
    }
    auto diffuse =
        dynamic_cast<lambertian*>(PreviewMaterials(editor.roots[0].object.get())[0]);
    expect_true((diffuse && diffuse->preview_rough_model));
    expect_true((diffuse->preview_rough_model->A < 1));
    auto glass =
        dynamic_cast<dielectric*>(PreviewMaterials(editor.roots[3].object.get())[0]);
    expect_true((glass && glass->priority == 3));
    auto emitter =
        dynamic_cast<diffuse_light*>(PreviewMaterials(editor.roots[7].object.get())[0]);
    expect_true((emitter && emitter->invisible));
    auto fiber = dynamic_cast<hair*>(PreviewMaterials(editor.roots[9].object.get())[0]);
    expect_true((fiber && std::abs(fiber->alpha - 8) < 1e-6));
    hit_record sample{};
    sample.dpdu = vec3f(1, 0, 0);
    sample.dpdv = vec3f(0, 1, 0);
    sample.normal = normal3f(0, 0, 1);
    sample.v = .4;
    hair reference(fiber->sigma_a, fiber->eta, fiber->beta_m, fiber->beta_n, 8);
    Ray incoming(point3f(0, 0, 1), vec3f(.2, .3, -1), .5);
    const auto evaluated = fiber->f(incoming, sample, vec3f(.3, .4, 1));
    const auto expected = reference.f(incoming, sample, vec3f(.3, .4, 1));
    for (int channel = 0; channel < 3; ++channel) {
      expect_true(std::isfinite(evaluated[channel]));
      expect_true((std::abs(evaluated[channel] - expected[channel]) < 1e-5));
    }
    PreviewObjectState ui;
    editor.Describe(editor.roots[9].id, ui);
    EditorField(ui, "Hair color mode").values[0] = 2;
    EditorField(ui, "Pigment").values[0] = .3;
    EditorField(ui, "Red pigment").values[0] = .7;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    fiber = dynamic_cast<hair*>(PreviewMaterials(editor.roots[9].object.get())[0]);
    expect_true((std::abs(fiber->sigma_a[0] - (.3 * .419 + .7 * .187)) < 1e-6));
  }

  test_that(
      "color textures retain recipes and invalid joint edits preserve the scene") {
    auto input = EditorFixture(R"(
      native_editor_scene(process_scene(create_instances(sphere(material=diffuse(color='red')),x=c(-2,2)))$scene)
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged] {
        live = staged;
      };
    };
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    EditorField(ui, "Texture mode").values[0] = 1;
    EditorField(ui, "Secondary color").values = {0, 0, 1};
    EditorField(ui, "Checker period").values[0] = .7;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((EditorField(ui, "Texture mode").values[0] == 1));
    auto first =
        dynamic_cast<lambertian*>(PreviewMaterials(editor.roots[0].object.get())[0]);
    auto other =
        dynamic_cast<lambertian*>(PreviewMaterials(editor.roots[1].object.get())[0]);
    auto blue = first->albedo->value(.5, .5, point3f(.1));
    auto red = other->albedo->value(.5, .5, point3f(.1));
    expect_true((blue[2] == 1 && red[0] == 1));
    auto before = live;
    EditorField(ui, "Texture mode").values[0] = 4;
    EditorField(ui, "Gradient start XYZ").values = {2, 2, 2};
    EditorField(ui, "Gradient end XYZ").values = {2, 2, 2};
    ui.materials[0].changed = ui.apply_material = true;
    expect_false(editor.Apply(ui));
    expect_true((live == before));
    expect_true((ui.error.find("different points") != std::string::npos));
    EditorField(ui, "Gradient end XYZ").values = {3, 2, 2};
    ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((EditorField(ui, "Texture mode").values[0] == 4));
    EditorField(ui, "Texture mode").values[0] = 2;
    EditorField(ui, "Noise phase (degrees)").values[0] = 120;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((EditorField(ui, "Noise phase (degrees)").values[0] == 120));
  }

  test_that(
      "initial roughness files affect shading and retain mappings through edits") {
    auto input = EditorFixture(R"(
      asset = tempfile(fileext='.png')
      png::writePNG(matrix(rep(c(0,1), each=2), nrow=2), asset)
      scene = native_editor_scene(process_scene(rbind(
        sphere(material=microfacet(roughness=.18, eta=2, kappa=3,
          roughness_texture=asset, roughness_range=c(.03,.45))),
        sphere(x=3, material=microfacet(roughness=.18, eta=2, kappa=3,
          roughness_texture=asset, roughness_range=c(.2,.7), roughness_flip=TRUE))
      ))$scene)
      attr(scene,'test_asset') = asset
      scene
    )");
    const std::string path = Rcpp::as<std::string>(input.attr("test_asset"));
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged]() {
        live = staged;
      };
    };
    // A second material may reuse the bytes, but it must not remap that cache.
    expect_true((live->roughness_textures.size() == 2));
    expect_true((live->roughness_textures[0] == live->roughness_textures[1]));
    expect_true(
        (live->roughness_textures[0][0] == 0 && live->roughness_textures[0][3] == 255));
    auto reflection = [&](size_t row, Float u) {
      auto* material = PreviewMaterials(editor.roots[row].object.get())[0];
      hit_record rec{};
      rec.normal = normal3f(0, 0, 1);
      rec.dpdu = vec3f(1, 0, 0);
      rec.dpdv = vec3f(0, 1, 0);
      rec.u = u;
      rec.v = .5;
      return material->f(
          Ray(point3f(0, 0, 1), vec3f(0, 0, -1)), rec, vec3f(0, 0, 1))[0];
    };
    expect_true((reflection(0, .25) > 2 * reflection(0, .75)));
    expect_true((reflection(1, .75) > 2 * reflection(1, .25)));
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    expect_true((EditorField(ui, "Roughness map file").text == path));
    expect_true((EditorField(ui, "Roughness map range").values[0] == Approx(.03)));
    expect_true((EditorField(ui, "Roughness map range").values[1] == Approx(.45)));
    expect_false(EditorField(ui, "Flip roughness map").values[0]);
    EditorField(ui, "Color").values = {.2, .4, .6};
    EditorField(ui, "Color").changed = true;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true(ui.error.empty());
    expect_true((reflection(0, .25) > 2 * reflection(0, .75)));
    expect_true((EditorField(ui, "Roughness map file").text == path));
    expect_true((EditorField(ui, "Roughness map range").values[1] == Approx(.45)));
    editor.Describe(editor.roots[1].id, ui);
    expect_true((EditorField(ui, "Roughness map range").values[0] == Approx(.2)));
    expect_true((EditorField(ui, "Roughness map range").values[1] == Approx(.7)));
    expect_true(EditorField(ui, "Flip roughness map").values[0]);
    expect_true((reflection(1, .75) > 2 * reflection(1, .25)));
  }

  test_that("texture files apply atomically and stay isolated to one instance") {
    auto input = EditorFixture(R"(
      asset=tempfile(fileext='.png')
      image=array(0,dim=c(8,8,3)); image[,,3]=1
      png::writePNG(image,asset)
      scene=native_editor_scene(process_scene(create_instances(
        sphere(material=microfacet(roughness=.2,eta=2)), x=c(-2,2)))$scene)
      attr(scene,'test_asset')=asset
      scene
    )");
    const std::string path = Rcpp::as<std::string>(input.attr("test_asset"));
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged] {
        live = staged;
      };
    };
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    EditorField(ui, "Texture mode").values[0] = 5;
    EditorField(ui, "Color texture file").text = path;
    EditorField(ui, "Use roughness map").values[0] = 1;
    EditorField(ui, "Roughness map file").text = path;
    EditorField(ui, "Roughness map range").values = {.1, .4, 0};
    EditorField(ui, "Use bump map").values[0] = 1;
    EditorField(ui, "Bump map file").text = path;
    EditorField(ui, "Bump intensity").values[0] = .2;
    EditorField(ui, "Use alpha map").values[0] = 1;
    EditorField(ui, "Alpha map file").text = path;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true(ui.error.empty());
    expect_true((EditorField(ui, "Color texture file").text == path));
    expect_true((EditorField(ui, "Bump map file").text == path));
    expect_true((EditorField(ui, "Alpha map file").text == path));
    expect_true((EditorField(ui, "Roughness map file").text == path));
    auto material = PreviewMaterials(editor.roots[0].object.get())[0];
    expect_true((material->preview_maps && material->preview_maps->bump_enabled));
    auto sibling = PreviewMaterials(editor.roots[1].object.get())[0];
    expect_false(bool(sibling->preview_maps));
    hit_record rec{};
    rec.u = rec.v = .5;
    rec.p = point3f(0);
    const auto color = material->get_albedo(rec);
    expect_true((color[2] > .99 && color[0] < .01));
    const auto original = sibling->get_albedo(rec);
    expect_true((original[0] > .99));
    auto saved = editor.ExportEdits();
    Rcpp::List entry = saved[0];
    Rcpp::List materials = entry["materials"], edit = materials[0],
               values = edit["values"];
    expect_true((Rcpp::as<std::string>(values["Color texture file"]) == path));
    PreviewScene restored(input);
    restored.ImportEdits(input, saved);
    EditorBuild replay(input, restored);
    auto restored_material = PreviewMaterials(restored.roots[0].object.get())[0];
    expect_true(((restored_material->get_albedo(rec) - color).length() < 1e-6));
    expect_true(restored_material->preview_maps->bump_enabled);
    restored.Describe(restored.roots[0].id, ui);
    expect_true((EditorField(ui, "Roughness map range").values[0] == .1));
    auto before = live;
    EditorField(ui, "Color texture file").text = path + ".missing";
    EditorField(ui, "Color texture file").changed = true;
    ui.materials[0].changed = ui.apply_material = true;
    expect_false(editor.Apply(ui));
    expect_true((live == before));
    expect_true((ui.error.find("does not exist") != std::string::npos));
    editor.Describe(ui.id, ui);
    EditorField(ui, "Use bump map").values[0] = 0;
    EditorField(ui, "Use bump map").changed = true;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_false(EditorField(ui, "Use bump map").values[0]);
    EditorField(ui, "Use bump map").values[0] = 1;
    EditorField(ui, "Use bump map").changed = true;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((EditorField(ui, "Bump map file").text == path));
  }

  test_that("mesh alpha map changes refresh cached opaque shadow classification") {
    auto input = EditorFixture(R"(
      path=tempfile(fileext='.obj')
      writeLines(c('v -1 -1 0','v 1 -1 0','v 1 1 0','v -1 1 0',
        'vt 0 0','vt 1 0','vt 1 1','vt 0 1','f 1/1 2/2 3/3','f 1/1 3/3 4/4'),path)
      asset=tempfile(fileext='.png'); png::writePNG(matrix(0,8,8),asset)
      scene=native_editor_scene(process_scene(obj_model(path,load_material=FALSE))$scene)
      attr(scene,'test_asset')=asset
      scene
    )");
    PreviewScene editor(input);
    auto live = std::make_shared<EditorBuild>(input, editor);
    editor.prepare_rebuild = [&](PreviewScene& next) {
      auto staged = std::make_shared<EditorBuild>(input, next);
      return [&, staged] {
        live = staged;
      };
    };
    const auto before = live->root->ShadowType();
    expect_true((before != OpaqueShadowType::Unsupported));
    PreviewObjectState ui;
    editor.Describe(editor.roots[0].id, ui);
    EditorField(ui, "Use alpha map").values[0] = 1;
    EditorField(ui, "Alpha map file").text =
        Rcpp::as<std::string>(input.attr("test_asset"));
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((live->root->ShadowType() == OpaqueShadowType::Unsupported));
    EditorField(ui, "Use alpha map").values[0] = 0;
    EditorField(ui, "Use alpha map").changed = true;
    ui.materials[0].changed = ui.apply_material = true;
    expect_true(editor.Apply(ui));
    expect_true((live->root->ShadowType() == before));
  }
}
#endif
