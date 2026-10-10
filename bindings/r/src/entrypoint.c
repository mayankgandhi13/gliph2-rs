// Forwards routine registration from C to Rust, so the linker keeps the
// static library's symbols.
void R_init_gliph2rs_extendr(void *dll);

void R_init_gliph2rs(void *dll) {
    R_init_gliph2rs_extendr(dll);
}
