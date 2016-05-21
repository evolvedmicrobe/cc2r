#include <Rcpp.h>

#include <string>
#include <tuple>
#include <vector>

#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Mutation.h>
#include <pacbio/consensus/Sequence.h>
#include <pacbio/consensus/ModelFactory.h>
#include <pacbio/consensus/ModelConfig.h>

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


//' Get the score for a read/template pair under a given model at a particular
//' SNR value.
//'
//' @param ref A string that is the reference window we are aligning to.
//' @param read A list with elements name/read.
//'
//' @return Returns a list with the aligned read, ref and score
//' @export
// [[Rcpp::export]]
DataFrame getTransitionParameters(std::string tpl, std::string mdl, NumericVector snrs) {
  if (snrs.length() != 4) {
    stop("SNR must have 4 entries!");
  }
  const SNR snr(snrs.at(0), snrs.at(1), snrs.at(2), snrs.at(3));
  auto model = ModelFactory::Create(mdl, snr);
  std::vector<TemplatePosition> tparams =  model->Populate(tpl);
  size_t length = tpl.size();
  NumericVector m(length), s(length), b(length), d(length);
  CharacterVector bp(length);
  for(int i=0; i < length; i++) {
    TemplatePosition pos = tparams.at(i);
    m[i] = pos.Match;
    b[i] = pos.Branch;
    s[i] = pos.Stick;
    d[i] = pos.Deletion;
    bp[i] = pos.Base;
  }

  return DataFrame::create(_["match"] = m,
                           _["branch"] = b,
                           _["stick"] = s,
                           _["delete"] = d,
                           _["bp"] = bp);
}




//' Get the score for a read/template pair under a given model at a particular
//' SNR value.
//'
//' @param ref A string that is the reference window we are aligning to.
//' @param read A list with elements name/read.
//'
//' @return Returns a list with the aligned read, ref and score
//' @export
// [[Rcpp::export]]
DataFrame getScore(std::string read, std::string tpl, std::string mdl, std::vector<int>& pws, NumericVector snrs) {
  if (read.size() != pws.size()) {
    stop("Read and Pulse Width had Different Sizes");
  }
  if (snrs.length() != 4) {
    stop("SNR must have 4 entries!");
  }

  const SNR snr(snrs.at(0), snrs.at(1), snrs.at(2), snrs.at(3));

  const IntegratorConfig cfg(std::numeric_limits<double>::quiet_NaN());
  std::vector<uint8_t> pw(pws.begin(), pws.end());

  MonoMolecularIntegrator ai(tpl, cfg, snr, mdl);
  ai.AddRead(MappedRead(MkRead(read, snr, mdl, pw), StrandEnum::FORWARD, 0,
                        tpl.length(), true, true));

  NumericVector ll(1);
  NumericVector mean(1);
  NumericVector var(1);
  auto params = ai.NormalParameters();;
  if (params.size() == 1) {
    auto norms = params.at(0);
    ll[0] = ai.LL();
    mean[0] = norms.first;
    var[0] = norms.second;
  } else {
    ll[0] = NA_REAL;
    mean[0] = NA_REAL;
    var[0] = NA_REAL;
  }
  // Rcout << norms.first << " " << norms.second << endl;
  // Rcout << ai.ZScores()[0] << endl;
  return DataFrame::create(_["ll"] = ll,  _["mean"] = mean, _["var"] = var);
}


