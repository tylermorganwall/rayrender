#include "animation_rcpp.h"
#include "math/quaternion.h"

std::shared_ptr<const animation::CameraTrack> read_camera_track(const Rcpp::DataFrame &input,
                                                              double fps, double start) {
  if (!std::isfinite(fps) || fps <= 0 || !std::isfinite(start))
    Rcpp::stop("Camera frame rate and start time must be finite; fps must be positive.");
  const char *columns[] = {"x", "y", "z", "dx", "dy", "dz", "aperture", "fov",
                           "focal", "orthox", "orthoy", "upx", "upy", "upz"};
  std::vector<animation::CameraKey> keys(input.nrows());
  for (size_t i = 0; i < keys.size(); ++i) keys[i].time = start + double(i)/fps;
  if (input.containsElementNamed("time")) {
    auto times = Rcpp::as<std::vector<double>>(input["time"]);
    if (times.size() != keys.size()) Rcpp::stop("Camera time count does not match its rows.");
    for (size_t i = 0; i < keys.size(); ++i) keys[i].time = times[i];
  }
  for (size_t column = 0; column < 14; ++column) {
    if (!input.containsElementNamed(columns[column])) Rcpp::stop("Camera track is missing a required column.");
    auto values = Rcpp::as<std::vector<double>>(input[columns[column]]);
    if (values.size() != keys.size()) Rcpp::stop("Camera column length does not match its rows.");
    for (size_t i = 0; i < keys.size(); ++i) keys[i].values[column] = values[i];
  }
  return std::make_shared<const animation::CameraTrack>(std::move(keys));
}

namespace {
void check_schema(const Rcpp::List &input) {
  if (!input.containsElementNamed("version") || Rcpp::as<int>(input["version"]) != animation::schema_version)
    Rcpp::stop("Unsupported animation schema version.");
}
std::vector<animation::Value> read_values(SEXP data, size_t count, int dimension) {
  if (!Rf_isMatrix(data) || !Rf_isNumeric(data)) Rcpp::stop("Track values and tangents must be numeric matrices.");
  Rcpp::NumericMatrix matrix = Rcpp::as<Rcpp::NumericMatrix>(data);
  if (matrix.nrow() != count || matrix.ncol() != dimension) Rcpp::stop("Track value dimensions do not match the property and key count.");
  std::vector<animation::Value> result(count);
  for (size_t i = 0; i < count; ++i) for (int k = 0; k < dimension; ++k) result[i][k] = matrix(i, k);
  return result;
}
}

std::shared_ptr<const animation::Track> read_keyframe_track(const Rcpp::List &input) {
  check_schema(input);
  std::string property = Rcpp::as<std::string>(input["property"]), method = Rcpp::as<std::string>(input["interpolation"]);
  animation::Property p;
  int dimension;
  if (property == "position") { p = animation::Property::Position; dimension = 3; }
  else if (property == "rotation") { p = animation::Property::Rotation; dimension = 4; }
  else if (property == "scale") { p = animation::Property::Scale; dimension = 3; }
  else if (property == "scalar") { p = animation::Property::Scalar; dimension = 1; }
  else Rcpp::stop("Unknown animation property.");
  animation::Interpolation interpolation;
  if (method == "step") interpolation = animation::Interpolation::Step;
  else if (method == "linear") interpolation = animation::Interpolation::Linear;
  else if (method == "cubic") interpolation = animation::Interpolation::Cubic;
  else if (method == "smooth") interpolation = animation::Interpolation::Smooth;
  else Rcpp::stop("Unknown animation interpolation.");
  if (!Rf_isNumeric(input["time"])) Rcpp::stop("Track times must be numeric.");
  auto times = Rcpp::as<std::vector<double>>(input["time"]);
  auto values = read_values(input["value"], times.size(), dimension);
  std::vector<animation::Value> incoming, outgoing;
  if (!Rf_isNull(input["incoming"])) incoming = read_values(input["incoming"], times.size(), dimension);
  if (!Rf_isNull(input["outgoing"])) outgoing = read_values(input["outgoing"], times.size(), dimension);
  return std::make_shared<const animation::Track>(p, interpolation, std::move(times), std::move(values), std::move(incoming), std::move(outgoing));
}

std::shared_ptr<const animation::TransformTrack> read_transform_track(const Rcpp::List &input) {
  check_schema(input);
  auto result = std::make_shared<animation::TransformTrack>();
  for (const char *name : {"position", "rotation", "scale"}) {
    if (Rf_isNull(input[name])) continue;
    auto track = read_keyframe_track(Rcpp::as<Rcpp::List>(input[name]));
    if (std::string(name) == "position" && track->property == animation::Property::Position) result->position = track;
    else if (std::string(name) == "rotation" && track->property == animation::Property::Rotation) result->rotation = track;
    else if (std::string(name) == "scale" && track->property == animation::Property::Scale) result->scale = track;
    else Rcpp::stop("Transform channel target does not match its track property.");
  }
  auto pivot = Rcpp::as<std::vector<double>>(input["pivot"]);
  if (pivot.size() != 3) Rcpp::stop("Pivot must contain three finite numbers.");
  for (int k = 0; k < 3; ++k) {
    if (!std::isfinite(pivot[k]) || std::abs(pivot[k]) > 1e12)
      Rcpp::stop("Pivot must contain three finite numbers with magnitude at most 1e12.");
    result->pivot[k] = pivot[k];
  }
  if (input.containsElementNamed("player") && !Rf_isNull(input["player"])) {
    Rcpp::List player = input["player"];
    check_schema(player);
    if (!Rf_isNull(player["clip"])) Rcpp::stop("A bare transform track has no named clips; use clip = NULL.");
    if (Rcpp::as<std::string>(player["root_motion"]) != "preserve")
      Rcpp::stop("Root-motion removal requires an explicit hierarchy and locomotion root.");
    std::string mode = Rcpp::as<std::string>(player["loop"]);
    animation::Loop loop;
    if (mode == "hold") loop = animation::Loop::Hold;
    else if (mode == "repeat") loop = animation::Loop::Repeat;
    else if (mode == "pingpong") loop = animation::Loop::PingPong;
    else Rcpp::stop("Unknown clip loop mode.");
    double begin = INFINITY, end = -INFINITY;
    for (const auto &channel : {result->position, result->rotation, result->scale}) {
      if (!channel) continue;
      begin = std::min(begin, channel->times.front());
      end = std::max(end, channel->times.back());
    }
    if (!std::isfinite(begin)) begin = end = 0;
    auto scalar = [&](const char *name) {
      SEXP value = player[name];
      if (!Rf_isNumeric(value) || Rf_xlength(value) != 1) Rcpp::stop("Clip clock fields must be finite scalars.");
      return Rcpp::as<double>(value);
    };
    result->clock = std::make_shared<const animation::ClipClock>(begin, end,
      scalar("start_time"), scalar("offset"), scalar("speed"), loop);
  }
  return result;
}

// Internal authoring validation and a headless oracle for the render evaluator.
// [[Rcpp::export]]
Rcpp::NumericMatrix animation_evaluate_track(Rcpp::List track, Rcpp::NumericVector times) {
  auto compiled = read_keyframe_track(track);
  Rcpp::NumericMatrix result(times.size(), compiled->Dimensions());
  for (R_xlen_t i = 0; i < times.size(); ++i) {
    auto v = compiled->At(times[i]);
    for (int k = 0; k < compiled->Dimensions(); ++k) result(i, k) = v[k];
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List animation_evaluate_transform(Rcpp::List track, Rcpp::NumericVector times) {
  auto compiled = read_transform_track(track);
  Rcpp::List result(times.size());
  for (R_xlen_t i = 0; i < times.size(); ++i) {
    if (!std::isfinite(times[i])) Rcpp::stop("Evaluation time must be finite.");
    auto matrix = compiled->At(times[i]).GetMatrix();
    Rcpp::NumericMatrix m(4, 4);
    for (int row = 0; row < 4; ++row) for (int col = 0; col < 4; ++col) m(row, col) = matrix.m[row][col];
    result[i] = m;
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix animation_euler_quaternions(Rcpp::NumericMatrix angles) {
  if (angles.ncol() != 3) Rcpp::stop("Euler angles require three columns.");
  Rcpp::NumericMatrix result(angles.nrow(), 4);
  for (int i = 0; i < angles.nrow(); ++i) {
    for (int k = 0; k < 3; ++k) if (!std::isfinite(angles(i, k))) Rcpp::stop("Euler angles must be finite.");
    Quaternion q(RotateZ(angles(i, 2)) * RotateY(angles(i, 1)) * RotateX(angles(i, 0)));
    q = Normalize(q);
    for (int k = 0; k < 3; ++k) result(i, k) = q.v[k];
    result(i, 3) = q.w;
  }
  return result;
}
