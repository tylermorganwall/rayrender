/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; see LICENSE.protocol.
 * Optional R adapter, copied alongside the pure C header by optional consumers.
 * Keep the SEXP rooted for the complete native session; do not cache the table.
 */
#ifndef RIMGUI_R_H
#define RIMGUI_R_H
#include <Rinternals.h>
#include "rimgui_api.h"
static inline int32_t rimgui_api_from_R_v1(SEXP handle, uint32_t required_size,
    uint64_t caps, const rimgui_api_v1 **out) {
    const void *ptr;
    if (!out) return RIMGUI_INVALID;
    *out = NULL;
    if (TYPEOF(handle) != EXTPTRSXP ||
        R_ExternalPtrTag(handle) != Rf_install("rimgui.api.v1") ||
        !Rf_inherits(handle, "rimgui_api_v1")) return RIMGUI_ABI;
    ptr = R_ExternalPtrAddr(handle);
    if (rimgui_validate_api_v1(ptr, sizeof(rimgui_api_header_v1), required_size, caps))
        return RIMGUI_ABI;
    {
        const rimgui_api_v1 *table = (const rimgui_api_v1 *)ptr;
#define RIMGUI_CHECK_MEMBER(member) \
        if (required_size >= RIMGUI_MEMBER_SIZE(member) && !table->member) return RIMGUI_ABI
        RIMGUI_CHECK_MEMBER(open); RIMGUI_CHECK_MEMBER(close);
        RIMGUI_CHECK_MEMBER(register_callback); RIMGUI_CHECK_MEMBER(unregister_callback);
        RIMGUI_CHECK_MEMBER(step); RIMGUI_CHECK_MEMBER(request_close);
        RIMGUI_CHECK_MEMBER(texture_create); RIMGUI_CHECK_MEMBER(texture_update);
        RIMGUI_CHECK_MEMBER(texture_destroy); RIMGUI_CHECK_MEMBER(widget);
        RIMGUI_CHECK_MEMBER(gizmo); RIMGUI_CHECK_MEMBER(viewport);
        RIMGUI_CHECK_MEMBER(event_push); RIMGUI_CHECK_MEMBER(event_poll);
        RIMGUI_CHECK_MEMBER(input);
#undef RIMGUI_CHECK_MEMBER
        *out = table;
    }
    return RIMGUI_OK;
}
#endif
