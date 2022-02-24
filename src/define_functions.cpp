// [[Rcpp::export]]
bool has_gui_capability() {
#ifdef RAY_HAS_X11
  return(true);
#else
  return(false);
#endif
}