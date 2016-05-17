#include <Rcpp.h>

#include <string>
#include <tuple>
#include <vector>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/Sequence.h>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using namespace PacBio::Consensus;  // NOLINT
using namespace Rcpp;

Read MkRead(const std::string& seq, const SNR& snr, const std::string& mdl, const std::vector<uint8_t>& pw)
{
  std::vector<uint8_t> ipd(0, seq.length());
  return Read("NA", seq, ipd, pw, snr, mdl);
}


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {

  auto mdl = "P6-C4";
  const string longTpl =
    "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCAT"
    "AACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGCTTTTTGGCC"
    "TCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATT"
    "TTCGGTGCGAGTATCCGTACCATTCAGAACT";
  const string longRead =
    "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCAT"
    "AACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGCTTTTTGGCC"
    "TCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATT"
    "TTCGGTGCGAGTATCCGTACCATTCAGAACT";
  const SNR snr(10, 7, 5, 11);

  const IntegratorConfig cfg(std::numeric_limits<double>::quiet_NaN());
  std::vector<uint8_t> pw(longTpl.length(), 2);
  MonoMolecularIntegrator ai(longTpl, cfg, snr, mdl);
  ai.AddRead(MappedRead(MkRead(longRead, snr, mdl, pw), StrandEnum::FORWARD, 0,
                        longTpl.length(), true, true));
  return ai.LL();
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
