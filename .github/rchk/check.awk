# Compare function, diagnostic and source line, ignoring the installation prefix.
# The first input is the baseline for the installed Rcpp version (or /dev/null).
/^Function / { fn = $0; next }
/^  \[[A-Z][A-Z]\] / {
    # R-hub treats these two analysis limitations as informational.
    if (/unsupported form of unprotect|has address taken, results will be incomplete/) next

    diagnostic = $0
    sub(/ \/[^ ]*\/Rcpp\/include\/Rcpp\//, " Rcpp/", diagnostic)
    key = fn SUBSEP diagnostic
    if (FILENAME == ARGV[1]) known[key] = 1
    else if (!(key in known)) {
        print fn "\n" $0
        failed = 1
    }
}
/^Analyzed [1-9][0-9]* functions, traversed [1-9][0-9]* states\.$/ {
    if (FILENAME == ARGV[2]) complete = 1
}
END {
    if (!complete) print "Missing or incomplete rchk report"
    exit (failed || !complete)
}
