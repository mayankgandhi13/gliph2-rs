## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The local check also reports two notes that are artefacts of the checking
machine rather than the package: "unable to verify current time" (no network
access) and HTML validation messages about `<main>` from an older `tidy`.

## Rust

This package contains Rust source code and follows the CRAN policy in
"Using Rust in CRAN packages":

* All Cargo dependencies are vendored in `src/rust/vendor.tar.xz` (18 crates,
  produced by `cargo vendor` and unmodified, so their checksums verify). The
  build unpacks them and compiles with `--offline`, so nothing is downloaded
  at install time.
* `cargo build` is invoked with `-j 2`.
* `configure` and `configure.win` check for `cargo` and `rustc`, look in
  `~/.cargo/bin` as well as on the `PATH`, verify the minimum supported Rust
  version declared in `SystemRequirements`, and report both versions in the
  installation log before compilation starts.
* The vendored crates are licensed MIT or Apache-2.0; their authors and
  licences are listed in `inst/AUTHORS`.

## Test environments

* macOS 15.1 (aarch64), R 4.4.2, rustc 1.95.0
* Rocky Linux 9 (x86_64), R 4.4, rustc 1.98.1
