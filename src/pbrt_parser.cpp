// PBRT lexical/statement grammar, using PEGTL from the piton R package.
// Numeric arrays become numeric vectors during parsing, avoiding the character
// token copies that dominate memory for large PBRT volume grids.
#include <R_ext/Utils.h>
#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <pegtl.hpp>
#include <sstream>
namespace pbrt_peg {
using namespace tao::pegtl;
struct comment : seq<one<'#'>, until<eolf>> {};
struct spacing : star<sor<space, comment>> {};
struct quote_start : one<'"'> {};
struct quote_end : one<'"'> {};
struct quoted
    : seq<quote_start, must<star<sor<seq<one<'\\'>, any>, not_one<'"', '\n', '\r'>>>, quote_end>> {
};
struct number : seq<opt<one<'+', '-'>>,
                    sor<seq<plus<digit>, opt<one<'.'>, star<digit>>>, seq<one<'.'>, plus<digit>>>,
                    opt<one<'e', 'E'>, opt<one<'+', '-'>>, plus<digit>>,
                    at<sor<space, one<'[', ']', '"', '#'>, eof>>> {};
struct bare : plus<not_one<' ', '\t', '\n', '\r', '[', ']', '"', '#'>> {};
struct array_start : one<'['> {};
struct array_end : one<']'> {};
struct array
    : seq<array_start, spacing, must<star<seq<sor<quoted, number, bare>, spacing>>, array_end>> {};
struct reserved
    : sor<TAO_PEGTL_STRING("All"), TAO_PEGTL_STRING("StartTime"), TAO_PEGTL_STRING("EndTime"),
          TAO_PEGTL_STRING("Inf"), TAO_PEGTL_STRING("NaN")> {};
struct directive : seq<not_at<seq<reserved, not_at<alpha>>>, range<'A', 'Z'>, plus<alpha>,
                       at<sor<space, one<'[', ']', '"', '#'>, eof>>> {};
// A statement begins with a capitalized directive. Its operand arity and typed
// parameter contracts are validated by R against PBRT's directive table.
struct operand : sor<quoted, array, seq<not_at<directive>, sor<number, bare>>> {};
struct statement : seq<directive, spacing, star<seq<operand, spacing>>> {};
struct document : must<spacing, star<statement>, eof> {};
struct State {
  std::vector<std::string> tokens;
  std::vector<int> lines;
  std::vector<Rcpp::RObject> arrays;
  std::vector<double> numbers;
  std::vector<std::string> strings;
  bool in_array = false, numeric = true;
  std::string filename, directive_name;
  int directive_line = 1, quote_line = 1;
  size_t values_seen = 0;
  void token(const std::string &s, int line) {
    tokens.push_back(s);
    lines.push_back(line);
  }
  void text(const std::string &s, int line) {
    if (!in_array) {
      token(s, line);
      return;
    }
    if (numeric) {
      strings.reserve(numbers.size() + 1);
      for (double n : numbers) {
        std::ostringstream o;
        o << std::setprecision(17) << n;
        strings.push_back(o.str());
      }
      numbers.clear();
      numeric = false;
    }
    strings.push_back(s);
  }
};
template <typename Rule> struct action : nothing<Rule> {};
template <> struct action<directive> {
  template <typename Input> static void apply(const Input &in, State &s) {
    s.directive_name = in.string();
    s.directive_line = in.position().line;
    s.token(in.string(), in.position().line);
  }
};
template <> struct action<quote_start> {
  template <typename Input> static void apply(const Input &in, State &s) {
    s.quote_line = in.position().line;
  }
};
template <> struct action<quoted> {
  template <typename Input> static void apply(const Input &in, State &s) {
    s.text(in.string(), in.position().line);
  }
};
template <> struct action<bare> {
  template <typename Input> static void apply(const Input &in, State &s) {
    const std::string text = in.string();
    if (text.compare(0, 7, "@array:") == 0)
      Rcpp::stop(s.filename + ":" + std::to_string(in.position().line) + ": Invalid bare token.");
    s.text(text, in.position().line);
  }
};
template <> struct action<number> {
  template <typename Input> static void apply(const Input &in, State &s) {
    if ((++s.values_seen & 65535) == 0)
      Rcpp::checkUserInterrupt();
    // Transform operands remain textual until R's operand validation; preserve
    // their spelling and full precision rather than double -> text rounding.
    if (s.in_array && s.numeric && s.directive_name != "Transform" &&
        s.directive_name != "ConcatTransform")
      // Match R's numeric conversion exactly, including extreme exponents.
      s.numbers.push_back(R_strtod(in.string().c_str(), nullptr));
    else
      s.text(in.string(), in.position().line);
  }
};
template <> struct action<array_start> {
  template <typename Input> static void apply(const Input &in, State &s) {
    s.in_array = true;
    s.numeric = true;
    s.numbers.clear();
    s.strings.clear();
    s.token("[", in.position().line);
  }
};
template <> struct action<array_end> {
  template <typename Input> static void apply(const Input &in, State &s) {
    s.arrays.push_back(s.numeric && !s.numbers.empty() ? Rcpp::wrap(s.numbers)
                                                       : Rcpp::wrap(s.strings));
    s.token("@array:" + std::to_string(s.arrays.size()), in.position().line);
    s.token("]", in.position().line);
    s.in_array = false;
  }
};
template <typename Rule> struct control : normal<Rule> {};
template <> struct control<eof> : normal<eof> {
  template <typename Input> static void raise(const Input &in, State &s) {
    // A gzip source is parsed from a temporary file; diagnostics must still
    // identify the user's original source path and line.
    Rcpp::stop(s.filename + ":" + std::to_string(in.position().line) +
               ": expected a PBRT directive.");
  }
};
template <> struct control<array_end> : normal<array_end> {
  template <typename Input> static void raise(const Input &, State &s) {
    Rcpp::stop(s.filename + ":" + std::to_string(s.directive_line) + ": " + s.directive_name +
               ": Unclosed array.");
  }
};
template <> struct control<quote_end> : normal<quote_end> {
  template <typename Input> static void raise(const Input &, State &s) {
    Rcpp::stop(s.filename + ":" + std::to_string(s.quote_line) + ": malformed quoted string.");
  }
};
} // namespace pbrt_peg
// [[Rcpp::export]]
Rcpp::List pbrt_lex_cpp(std::string filename, std::string source = "") {
  tao::pegtl::file_input<> input(filename);
  pbrt_peg::State state;
  state.filename = source.empty() ? filename : source;
  if (input.empty())
    Rcpp::stop("Empty PBRT file: " + state.filename);
  tao::pegtl::parse<pbrt_peg::document, pbrt_peg::action, pbrt_peg::control>(input, state);
  Rcpp::List arrays(state.arrays.size());
  for (size_t i = 0; i < state.arrays.size(); ++i)
    arrays[i] = state.arrays[i];
  return Rcpp::List::create(Rcpp::_["tokens"] = state.tokens, Rcpp::_["line"] = state.lines,
                            Rcpp::_["arrays"] = arrays);
}
