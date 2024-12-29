// [[Rcpp::export]]
bool has_gui_capability() {
#ifdef RAY_HAS_X11
  return(true);
#else
#ifdef RAY_WINDOWS
  return(true);
#else
  return(false);
#endif
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_OIDN() {
#ifdef HAS_OIDN
    return(true);
#else
    return(false);
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_NEON() {
#ifdef HAS_NEON
    return(true);
#else
    return(false);
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_SSE() {
#ifdef HAS_SSE
    return(true);
#else
    return(false);
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_SSE2() {
#ifdef HAS_SSE2
    return(true);
#else
    return(false);
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_SSE3() {
#ifdef HAS_SSE3
    return(true);
#else
    return(false);
#endif
}

// [[Rcpp::export]]
bool cppdef_HAS_SSE41() {
#ifdef HAS_SSE41
    return(true);
#else
    return(false);
#endif
}