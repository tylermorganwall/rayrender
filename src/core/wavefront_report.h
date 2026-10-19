#ifndef RAYRENDER_WAVEFRONT_REPORT_H
#define RAYRENDER_WAVEFRONT_REPORT_H

#include <Rcpp.h>
#include "wavefront.h"

inline Rcpp::List WavefrontReportToR(const WavefrontReport &report) {
  using Rcpp::_;
  return Rcpp::List::create(
    _["used"] = report.used, _["fallback"] = report.fallback,
    _["scene_reused"] = report.scene_reused,
    _["acceleration_builds"] = double(report.acceleration_builds),
    _["acceleration_refits"] = double(report.acceleration_refits),
    _["geometry_bytes_updated"] = double(report.geometry_bytes_updated),
    _["triangles"] = double(report.triangles),
    _["tessellated_primitives"] = double(report.tessellated_primitives),
    _["analytic_primitives"] = double(report.analytic_primitives),
    _["unique_meshes"] = double(report.unique_meshes),
    _["mesh_instances"] = double(report.mesh_instances),
    _["scattering_events"] = double(report.scattering_events),
    _["discarded_paths"] = double(report.discarded_paths),
    _["csg_failed_queries"] = double(report.csg_failed_queries),
    _["upload_seconds"] = report.upload_seconds,
    _["sample_seconds"] = report.sample_seconds,
    _["atmosphere_seconds"] = report.atmosphere_seconds,
    _["completed_samples"] = double(report.completed_samples));
}
#endif
