/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; see LICENSE.protocol.
 * Optional R adapter, copied alongside the pure C header by optional consumers.
 * Keep the SEXP rooted for the complete native session; do not cache the table.
 */
#ifndef RAYIMGUI_R_H
#define RAYIMGUI_R_H
#include <Rinternals.h>
#include "rayimgui_api.h"
static inline int32_t rayimgui_api_from_R_v1(SEXP handle, uint32_t required_size,
                                           uint64_t caps, const rayimgui_api_v1** out) {
  const void* ptr;
  if (!out) {
    return RAYIMGUI_INVALID;
  }
  *out = NULL;
  if (TYPEOF(handle) != EXTPTRSXP ||
      R_ExternalPtrTag(handle) != Rf_install("rayimgui.api.v1") ||
      !Rf_inherits(handle, "rayimgui_api_v1")) {
    return RAYIMGUI_ABI;
  }
  ptr = R_ExternalPtrAddr(handle);
  if (rayimgui_validate_api_v1(ptr, sizeof(rayimgui_api_header_v1), required_size, caps)) {
    return RAYIMGUI_ABI;
  }
  {
    const rayimgui_api_v1* table = (const rayimgui_api_v1*)ptr;
#define RAYIMGUI_CHECK_MEMBER(member)                                                    \
  if (required_size >= RAYIMGUI_MEMBER_SIZE(member) && !table->member)                   \
  return RAYIMGUI_ABI
    RAYIMGUI_CHECK_MEMBER(open);
    RAYIMGUI_CHECK_MEMBER(close);
    RAYIMGUI_CHECK_MEMBER(register_callback);
    RAYIMGUI_CHECK_MEMBER(unregister_callback);
    RAYIMGUI_CHECK_MEMBER(step);
    RAYIMGUI_CHECK_MEMBER(request_close);
    RAYIMGUI_CHECK_MEMBER(texture_create);
    RAYIMGUI_CHECK_MEMBER(texture_update);
    RAYIMGUI_CHECK_MEMBER(texture_destroy);
    RAYIMGUI_CHECK_MEMBER(widget);
    RAYIMGUI_CHECK_MEMBER(gizmo);
    RAYIMGUI_CHECK_MEMBER(viewport);
    RAYIMGUI_CHECK_MEMBER(event_push);
    RAYIMGUI_CHECK_MEMBER(event_poll);
    RAYIMGUI_CHECK_MEMBER(input);
#undef RAYIMGUI_CHECK_MEMBER
    *out = table;
  }
  return RAYIMGUI_OK;
}
#endif
