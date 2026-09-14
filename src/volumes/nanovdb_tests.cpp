#ifdef NOT_CRAN
#include <Rcpp.h>
#include <nanovdb/GridHandle.h>
#include <nanovdb/tools/GridChecksum.h>
#include <algorithm>
#include <atomic>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <testthat.h>
#include <thread>

namespace {
using Handles = std::vector<nanovdb::GridHandle<>>;

Handles read_uncompressed(const char *filename) {
  return nanovdb::io::readUncompressedGrids<nanovdb::GridHandle<>, std::vector>(filename);
}

// Supply enough bytes for the upstream reader's initial GridData probe, then
// let it rewind and diagnose the file header without touching grid storage.
struct HeaderStream : std::stringstream {
  explicit HeaderStream(const nanovdb::io::FileHeader &header) {
    std::string bytes(sizeof(nanovdb::GridData), '\0');
    std::memcpy(&bytes[0], &header, sizeof(header));
    str(bytes);
  }
  void skip(int64_t offset) { seekg(offset, std::ios_base::cur); }
};

Handles read_header(const nanovdb::io::FileHeader &header) {
  HeaderStream stream(header);
  return nanovdb::io::readUncompressedGrids<nanovdb::GridHandle<>, HeaderStream,
                                         std::vector>(stream);
}
} // namespace

context("NanoVDB host error handling") {
  test_that("range callbacks stay on the caller and propagate exceptions") {
    const auto caller = std::this_thread::get_id();
    std::vector<int> visits(1024, 0);
    std::atomic<bool> same_thread(true);
    nanovdb::util::forEach(0, visits.size(), 1, [&](const nanovdb::util::Range1D &range) {
      if (std::this_thread::get_id() != caller)
        same_thread = false;
      for (auto i = range.begin(); i != range.end(); ++i)
        ++visits[i];
    });
    expect_true(same_thread.load());
    expect_true(std::all_of(visits.begin(), visits.end(), [](int n) { return n == 1; }));

    auto fail = [](const nanovdb::util::Range1D &) {
      throw std::runtime_error("NanoVDB callback failure");
    };
    expect_error_as(nanovdb::util::forEach(0, 1024, 1, fail), std::runtime_error);
    nanovdb::util::forEach(0, 0, 1, fail); // An empty range must not call fail.
  }

  test_that("failed file opens unwind without closing a null FILE pointer") {
    Rcpp::Function tempfile = Rcpp::Environment::base_env()["tempfile"];
    const std::string filename = Rcpp::as<std::string>(tempfile()) + "/missing.nvdb";
    expect_error_as(read_uncompressed(filename.c_str()), std::runtime_error);
    expect_error_as(nanovdb::io::writeUncompressedGrids(filename.c_str(), Handles{}),
                    std::runtime_error);
  }

  test_that("invalid file headers throw instead of terminating the process") {
    nanovdb::io::FileHeader header{0, nanovdb::Version(), 1, nanovdb::io::Codec::NONE};
    expect_error_as(read_header(header), std::invalid_argument);
    header.magic = NANOVDB_MAGIC_FILE;
    header.version = nanovdb::Version(NANOVDB_MAJOR_VERSION_NUMBER + 1, 0, 0);
    expect_error_as(read_header(header), std::invalid_argument);
    header.version = nanovdb::Version();
    header.codec = nanovdb::io::Codec::ZIP;
    expect_error_as(read_header(header), std::invalid_argument);
  }

  test_that("full checksum validation still detects changes in grid data") {
    Rcpp::Function test_path = Rcpp::Environment::namespace_env("testthat")["test_path"];
    const std::string filename = Rcpp::as<std::string>(test_path("fixtures", "volumes", "tiles.nvdb"));
    auto handles = read_uncompressed(filename.c_str());
    expect_true(handles.size() == 1);
    auto &handle = handles.at(0);
    auto *grid = reinterpret_cast<nanovdb::GridData *>(handle.data());
    auto *bytes = static_cast<uint8_t *>(handle.data());
    // Validate the existing upstream checksum before explicitly enabling full
    // checking (the format also permits files without a checksum).
    expect_true(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));
    nanovdb::tools::updateChecksum(grid, nanovdb::CheckMode::Full);
    expect_true(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));
    bytes[grid->mGridSize - 1] ^= 1;
    expect_false(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));
    bytes[grid->mGridSize - 1] ^= 1;
    expect_true(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));

    // Exercise the legacy checksum algorithm with empty node levels as well.
    // It builds a NodeManager instead of hashing the payload in 4 KB blocks.
    grid->mVersion = nanovdb::Version(32, 6, 0);
    nanovdb::tools::updateChecksum(grid, nanovdb::CheckMode::Full);
    expect_true(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));
    bytes[grid->mGridSize - 1] ^= 1;
    expect_false(nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full));
  }
}
#endif
