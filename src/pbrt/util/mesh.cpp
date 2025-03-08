// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/util/mesh.h>

#include <pbrt/util/buffercache.h>
#include <pbrt/util/check.h>
#include <pbrt/util/error.h>
#include <pbrt/util/log.h>
#include <pbrt/util/print.h>
#include <pbrt/util/stats.h>
#include <pbrt/util/transform.h>

#include <rply/rply.h>

namespace pbrt {

STAT_RATIO("Geometry/Triangles per mesh", nTris, nTriMeshes);
STAT_MEMORY_COUNTER("Memory/Triangles", triangleBytes);

// TriangleMesh Method Definitions
TriangleMesh::TriangleMesh(const Transform &renderFromObject, bool reverseOrientation,
                           std::vector<int> indices, std::vector<Point3f> p,
                           std::vector<Vector3f> s, std::vector<Normal3f> n,
                           std::vector<Point2f> uv, std::vector<int> faceIndices,
                           Allocator alloc)
    : nTriangles(indices.size() / 3), nVertices(p.size()) {
    CHECK_EQ((indices.size() % 3), 0);
    ++nTriMeshes;
    nTris += nTriangles;
    triangleBytes += sizeof(*this);
    // Initialize mesh _vertexIndices_
    vertexIndices = intBufferCache->LookupOrAdd(indices, alloc);

    // Transform mesh vertices to rendering space and initialize mesh _p_
    for (Point3f &pt : p)
        pt = renderFromObject(pt);
    this->p = point3BufferCache->LookupOrAdd(p, alloc);

    // Remainder of _TriangleMesh_ constructor
    this->reverseOrientation = reverseOrientation;
    this->transformSwapsHandedness = renderFromObject.SwapsHandedness();

    if (!uv.empty()) {
        CHECK_EQ(nVertices, uv.size());
        this->uv = point2BufferCache->LookupOrAdd(uv, alloc);
    }
    if (!n.empty()) {
        CHECK_EQ(nVertices, n.size());
        for (Normal3f &nn : n) {
            nn = renderFromObject(nn);
            if (reverseOrientation)
                nn = -nn;
        }
        this->n = normal3BufferCache->LookupOrAdd(n, alloc);
    }
    if (!s.empty()) {
        CHECK_EQ(nVertices, s.size());
        for (Vector3f &ss : s)
            ss = renderFromObject(ss);
        this->s = vector3BufferCache->LookupOrAdd(s, alloc);
    }

    if (!faceIndices.empty()) {
        CHECK_EQ(nTriangles, faceIndices.size());
        this->faceIndices = intBufferCache->LookupOrAdd(faceIndices, alloc);
    }

    // Make sure that we don't have too much stuff to be using integers to
    // index into things.
    CHECK_LE(p.size(), std::numeric_limits<int>::max());
    // We could be clever and check indices.size() / 3 if we were careful
    // to promote to a 64-bit int before multiplying by 3 when we look up
    // in the indices array...
    CHECK_LE(indices.size(), std::numeric_limits<int>::max());
}

std::string TriangleMesh::ToString() const {
    std::string np = "(nullptr)";
    return StringPrintf(
        "[ TriangleMesh reverseOrientation: %s transformSwapsHandedness: %s "
        "nTriangles: %d nVertices: %d vertexIndices: %s p: %s n: %s "
        "s: %s uv: %s faceIndices: %s ]",
        reverseOrientation, transformSwapsHandedness, nTriangles, nVertices,
        vertexIndices ? StringPrintf("%s", pstd::MakeSpan(vertexIndices, 3 * nTriangles))
                      : np,
        p ? StringPrintf("%s", pstd::MakeSpan(p, nVertices)) : np,
        n ? StringPrintf("%s", pstd::MakeSpan(n, nVertices)) : np,
        s ? StringPrintf("%s", pstd::MakeSpan(s, nVertices)) : np,
        uv ? StringPrintf("%s", pstd::MakeSpan(uv, nVertices)) : np,
        faceIndices ? StringPrintf("%s", pstd::MakeSpan(faceIndices, nTriangles)) : np);
}

static void PlyErrorCallback(p_ply, const char *message) {
    Error("PLY writing error: %s", message);
}

bool TriangleMesh::WritePLY(std::string filename) const {
    if (s)
        Warning(R"(%s: PLY mesh will be missing tangent vectors "S".)", filename);

    return pbrt::WritePLY(
        filename, pstd::span<const int>(vertexIndices, 3 * nTriangles),
        pstd::span<const int>(), pstd::span<const Point3f>(p, nVertices),
        pstd::span<const Normal3f>(n, n ? nVertices : 0),
        pstd::span<const Point2f>(uv, uv ? nVertices : 0),
        pstd::span<const int>(faceIndices, faceIndices ? nTriangles : 0));
}

bool WritePLY(std::string filename, pstd::span<const int> triIndices,
              pstd::span<const int> quadIndices, pstd::span<const Point3f> p,
              pstd::span<const Normal3f> n, pstd::span<const Point2f> uv,
              pstd::span<const int> faceIndices) {
    p_ply plyFile =
        ply_create(filename.c_str(), PLY_DEFAULT, PlyErrorCallback, 0, nullptr);
    if (!plyFile)
        return false;

    int nVertices = p.size();
    ply_add_element(plyFile, "vertex", nVertices);
    ply_add_scalar_property(plyFile, "x", PLY_FLOAT);
    ply_add_scalar_property(plyFile, "y", PLY_FLOAT);
    ply_add_scalar_property(plyFile, "z", PLY_FLOAT);
    if (!n.empty()) {
        ply_add_scalar_property(plyFile, "nx", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "ny", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "nz", PLY_FLOAT);
    }
    if (!uv.empty()) {
        ply_add_scalar_property(plyFile, "u", PLY_FLOAT);
        ply_add_scalar_property(plyFile, "v", PLY_FLOAT);
    }

    int nTriangles = triIndices.size() / 3, nQuads = quadIndices.size() / 4;
    ply_add_element(plyFile, "face", nTriangles + nQuads);
    ply_add_list_property(plyFile, "vertex_indices", PLY_UINT8, PLY_INT);
    if (!faceIndices.empty())
        ply_add_scalar_property(plyFile, "face_indices", PLY_INT);

    ply_write_header(plyFile);

    for (int i = 0; i < nVertices; ++i) {
        ply_write(plyFile, p[i].x);
        ply_write(plyFile, p[i].y);
        ply_write(plyFile, p[i].z);
        if (!n.empty()) {
            ply_write(plyFile, n[i].x);
            ply_write(plyFile, n[i].y);
            ply_write(plyFile, n[i].z);
        }
        if (!uv.empty()) {
            ply_write(plyFile, uv[i].x);
            ply_write(plyFile, uv[i].y);
        }
    }

    for (int i = 0; i < nTriangles; ++i) {
        ply_write(plyFile, 3);
        ply_write(plyFile, triIndices[3 * i]);
        ply_write(plyFile, triIndices[3 * i + 1]);
        ply_write(plyFile, triIndices[3 * i + 2]);
        if (!faceIndices.empty())
            ply_write(plyFile, faceIndices[i]);
    }
    for (int i = 0; i < nQuads; ++i) {
        ply_write(plyFile, 4);
        ply_write(plyFile, quadIndices[4 * i]);
        ply_write(plyFile, quadIndices[4 * i + 1]);
        ply_write(plyFile, quadIndices[4 * i + 2]);
        ply_write(plyFile, quadIndices[4 * i + 3]);
        if (!faceIndices.empty())
            ply_write(plyFile, faceIndices[i]);
    }

    ply_close(plyFile);
    return true;
}

STAT_RATIO("Geometry/Bilinear patches per mesh", nBlps, nBilinearMeshes);
STAT_MEMORY_COUNTER("Memory/Bilinear patches", blpBytes);

BilinearPatchMesh::BilinearPatchMesh(const Transform &renderFromObject,
                                     bool reverseOrientation, std::vector<int> indices,
                                     std::vector<Point3f> P, std::vector<Normal3f> N,
                                     std::vector<Point2f> UV, std::vector<int> fIndices,
                                     PiecewiseConstant2D *imageDist, Allocator alloc)
    : reverseOrientation(reverseOrientation),
      transformSwapsHandedness(renderFromObject.SwapsHandedness()),
      nPatches(indices.size() / 4),
      nVertices(P.size()),
      imageDistribution(imageDist) {
    CHECK_EQ((indices.size() % 4), 0);
    ++nBilinearMeshes;
    nBlps += nPatches;

    // Make sure that we don't have too much stuff to be using integers to
    // index into things.
    CHECK_LE(P.size(), std::numeric_limits<int>::max());
    CHECK_LE(indices.size(), std::numeric_limits<int>::max());

    vertexIndices = intBufferCache->LookupOrAdd(indices, alloc);

    blpBytes += sizeof(*this);

    // Transform mesh vertices to rendering space
    for (Point3f &p : P)
        p = renderFromObject(p);
    p = point3BufferCache->LookupOrAdd(P, alloc);

    // Copy _UV_ and _N_ vertex data, if present
    if (!UV.empty()) {
        CHECK_EQ(nVertices, UV.size());
        uv = point2BufferCache->LookupOrAdd(UV, alloc);
    }
    if (!N.empty()) {
        CHECK_EQ(nVertices, N.size());
        for (Normal3f &n : N) {
            n = renderFromObject(n);
            if (reverseOrientation)
                n = -n;
        }
        n = normal3BufferCache->LookupOrAdd(N, alloc);
    }

    if (!fIndices.empty()) {
        CHECK_EQ(nPatches, fIndices.size());
        faceIndices = intBufferCache->LookupOrAdd(fIndices, alloc);
    }
}

std::string BilinearPatchMesh::ToString() const {
    std::string np = "(nullptr)";
    return StringPrintf(
        "[ BilinearMatchMesh reverseOrientation: %s transformSwapsHandedness: "
        "%s "
        "nPatches: %d nVertices: %d vertexIndices: %s p: %s n: %s "
        "uv: %s faceIndices: %s ]",
        reverseOrientation, transformSwapsHandedness, nPatches, nVertices,
        vertexIndices ? StringPrintf("%s", pstd::MakeSpan(vertexIndices, 4 * nPatches))
                      : np,
        p ? StringPrintf("%s", pstd::MakeSpan(p, nVertices)) : np,
        n ? StringPrintf("%s", pstd::MakeSpan(n, nVertices)) : np,
        uv ? StringPrintf("%s", pstd::MakeSpan(uv, nVertices)) : np,
        faceIndices ? StringPrintf("%s", pstd::MakeSpan(faceIndices, nPatches)) : np);
}

struct FaceCallbackContext {
    int face[4];
    std::vector<int> triIndices, quadIndices;
};

void rply_message_callback(p_ply ply, const char *message) {
    Warning("rply: %s", message);
}

/* Callback to handle vertex data from RPly */
int rply_vertex_callback(p_ply_argument argument) {
    Float *buffer;
    long index, flags;

    ply_get_argument_user_data(argument, (void **)&buffer, &flags);
    ply_get_argument_element(argument, nullptr, &index);

    int stride = (flags & 0x0F0) >> 4;
    int offset = flags & 0x00F;

    buffer[index * stride + offset] = (float)ply_get_argument_value(argument);

    return 1;
}

/* Callback to handle face data from RPly */
int rply_face_callback(p_ply_argument argument) {
    FaceCallbackContext *context;
    long flags;
    ply_get_argument_user_data(argument, (void **)&context, &flags);

    long length, value_index;
    ply_get_argument_property(argument, nullptr, &length, &value_index);

    if (length != 3 && length != 4) {
        Warning("plymesh: Ignoring face with %i vertices (only triangles and quads "
                "are supported!)",
                (int)length);
        return 1;
    } else if (value_index < 0) {
        return 1;
    }

    if (value_index >= 0)
        context->face[value_index] = (int)ply_get_argument_value(argument);

    if (value_index == length - 1) {
        if (length == 3)
            for (int i = 0; i < 3; ++i)
                context->triIndices.push_back(context->face[i]);
        else {
            CHECK_EQ(length, 4);

            // Note: modify order since we're specifying it as a blp...
            context->quadIndices.push_back(context->face[0]);
            context->quadIndices.push_back(context->face[1]);
            context->quadIndices.push_back(context->face[3]);
            context->quadIndices.push_back(context->face[2]);
        }
    }

    return 1;
}

int rply_faceindex_callback(p_ply_argument argument) {
    std::vector<int> *faceIndices;
    long flags;
    ply_get_argument_user_data(argument, (void **)&faceIndices, &flags);

    faceIndices->push_back((int)ply_get_argument_value(argument));

    return 1;
}

TriQuadMesh TriQuadMesh::ReadPLY(const std::string &filename) {
    TriQuadMesh mesh;

    p_ply ply = ply_open(filename.c_str(), rply_message_callback, 0, nullptr);
    if (!ply)
        ErrorExit("Couldn't open PLY file \"%s\"", filename);

    if (ply_read_header(ply) == 0)
        ErrorExit("Unable to read the header of PLY file \"%s\"", filename);

    p_ply_element element = nullptr;
    size_t vertexCount = 0, faceCount = 0;

    /* Inspect the structure of the PLY file */
    while ((element = ply_get_next_element(ply, element)) != nullptr) {
        const char *name;
        long nInstances;

        ply_get_element_info(element, &name, &nInstances);
        if (strcmp(name, "vertex") == 0)
            vertexCount = nInstances;
        else if (strcmp(name, "face") == 0)
            faceCount = nInstances;
    }

    if (vertexCount == 0 || faceCount == 0)
        ErrorExit("%s: PLY file is invalid! No face/vertex elements found!", filename);

    mesh.p.resize(vertexCount);
    if (ply_set_read_cb(ply, "vertex", "x", rply_vertex_callback, mesh.p.data(), 0x30) ==
            0 ||
        ply_set_read_cb(ply, "vertex", "y", rply_vertex_callback, mesh.p.data(), 0x31) ==
            0 ||
        ply_set_read_cb(ply, "vertex", "z", rply_vertex_callback, mesh.p.data(), 0x32) ==
            0) {
        ErrorExit("%s: Vertex coordinate property not found!", filename);
    }

    mesh.n.resize(vertexCount);
    if (ply_set_read_cb(ply, "vertex", "nx", rply_vertex_callback, mesh.n.data(), 0x30) ==
            0 ||
        ply_set_read_cb(ply, "vertex", "ny", rply_vertex_callback, mesh.n.data(), 0x31) ==
            0 ||
        ply_set_read_cb(ply, "vertex", "nz", rply_vertex_callback, mesh.n.data(), 0x32) ==
            0)
        mesh.n.resize(0);

    /* There seem to be lots of different conventions regarding UV coordinate
     * names */
    mesh.uv.resize(vertexCount);
    if (((ply_set_read_cb(ply, "vertex", "u", rply_vertex_callback, mesh.uv.data(),
                          0x20) != 0) &&
         (ply_set_read_cb(ply, "vertex", "v", rply_vertex_callback, mesh.uv.data(),
                          0x21) != 0)) ||
        ((ply_set_read_cb(ply, "vertex", "s", rply_vertex_callback, mesh.uv.data(),
                          0x20) != 0) &&
         (ply_set_read_cb(ply, "vertex", "t", rply_vertex_callback, mesh.uv.data(),
                          0x21) != 0)) ||
        ((ply_set_read_cb(ply, "vertex", "texture_u", rply_vertex_callback,
                          mesh.uv.data(), 0x20) != 0) &&
         (ply_set_read_cb(ply, "vertex", "texture_v", rply_vertex_callback,
                          mesh.uv.data(), 0x21) != 0)) ||
        ((ply_set_read_cb(ply, "vertex", "texture_s", rply_vertex_callback,
                          mesh.uv.data(), 0x20) != 0) &&
         (ply_set_read_cb(ply, "vertex", "texture_t", rply_vertex_callback,
                          mesh.uv.data(), 0x21) != 0)))
        ;
    else
        mesh.uv.resize(0);

    FaceCallbackContext context;
    context.triIndices.reserve(faceCount * 3);
    context.quadIndices.reserve(faceCount * 4);
    if (ply_set_read_cb(ply, "face", "vertex_indices", rply_face_callback, &context, 0) ==
        0)
        ErrorExit("%s: vertex indices not found in PLY file", filename);

    if (ply_set_read_cb(ply, "face", "face_indices", rply_faceindex_callback,
                        &mesh.faceIndices, 0) != 0)
        mesh.faceIndices.reserve(faceCount);

    if (ply_read(ply) == 0)
        ErrorExit("%s: unable to read the contents of PLY file", filename);

    mesh.triIndices = std::move(context.triIndices);
    mesh.quadIndices = std::move(context.quadIndices);

    ply_close(ply);

    for (int idx : mesh.triIndices)
        if (idx < 0 || idx >= mesh.p.size())
            ErrorExit("plymesh: Vertex index %i is out of bounds! "
                      "Valid range is [0..%i)",
                      idx, int(mesh.p.size()));
    for (int idx : mesh.quadIndices)
        if (idx < 0 || idx >= mesh.p.size())
            ErrorExit("plymesh: Vertex index %i is out of bounds! "
                      "Valid range is [0..%i)",
                      idx, int(mesh.p.size()));

    return mesh;
}

void TriQuadMesh::ConvertToOnlyTriangles() {
    if (quadIndices.empty())
        return;

    triIndices.reserve(triIndices.size() + 3 * quadIndices.size() / 2);

    for (size_t i = 0; i < quadIndices.size(); i += 4) {
        triIndices.push_back(quadIndices[i]);  // 0, 1, 2 of original
        triIndices.push_back(quadIndices[i + 1]);
        triIndices.push_back(quadIndices[i + 3]);

        triIndices.push_back(quadIndices[i]);  // 0, 2, 3 of original
        triIndices.push_back(quadIndices[i + 3]);
        triIndices.push_back(quadIndices[i + 2]);
    }

    quadIndices.clear();
}

void TriQuadMesh::ComputeNormals() {
    n.resize(p.size());
    for (size_t i = 0; i < n.size(); ++i)
        n[i] = Normal3f(0, 0, 0);

    for (size_t i = 0; i < triIndices.size(); i += 3) {
        int v[3] = {triIndices[i], triIndices[i + 1], triIndices[i + 2]};
        Vector3f v10 = p[v[1]] - p[v[0]];
        Vector3f v21 = p[v[2]] - p[v[1]];

        Normal3f vn(Cross(v10, v21));
        if (LengthSquared(vn) > 0) {
            vn = Normalize(vn);
            n[v[0]] += vn;
            n[v[1]] += vn;
            n[v[2]] += vn;
        }
    }
    CHECK_EQ(0, quadIndices.size());  // TODO: handle this...

    for (size_t i = 0; i < n.size(); ++i)
        if (LengthSquared(n[i]) > 0)
            n[i] = Normalize(n[i]);
}

std::string TriQuadMesh::ToString() const {
    return StringPrintf("[ TriQuadMesh p: %s n: %s uv: %s faceIndices: %s triIndices: %s "
                        "quadIndices: %s ]",
                        p, n, uv, faceIndices, triIndices, quadIndices);
}

}  // namespace pbrt
