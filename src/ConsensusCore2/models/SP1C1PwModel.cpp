
#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <pacbio/consensus/ModelConfig.h>
#include <pacbio/consensus/Read.h>
#include <pacbio/consensus/ModelFactory.h>
#include "../Recursor.h"

namespace PacBio {
namespace Consensus {
namespace {

constexpr double kCounterWeight = 15.0;
const size_t OUTCOME_NUMBER = 12;
const size_t CONTEXT_NUMBER = 16;
class SP1C1PwModel : public ModelConfig
{
    REGISTER_MODEL(SP1C1PwModel);

public:
    static std::set<std::string> Names() { return {"S/P1-C1", "S/P2-C2/prospective-compatible"}; }
    SP1C1PwModel(const SNR& snr);
    std::unique_ptr<AbstractRecursor> CreateRecursor(std::unique_ptr<AbstractTemplate>&& tpl,
                                                     const MappedRead& mr, double scoreDiff) const;
    std::vector<TemplatePosition> Populate(const std::string& tpl) const;
    double ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;
    double ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const;



private:
    SNR snr_;
    double ctxTrans_[CONTEXT_NUMBER][4];
    double cachedEmissionExpectations_[CONTEXT_NUMBER][3][2];
    double ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const;
};

REGISTER_MODEL_IMPL(SP1C1PwModel);

// TODO(lhepler) comments regarding the CRTP
class SP1C1PwRecursor : public Recursor<SP1C1PwRecursor>
{
public:
    SP1C1PwRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                    double scoreDiff);
    static inline std::vector<uint8_t> EncodeRead(const MappedRead& read);
    static inline double EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr);
    virtual double UndoCounterWeights(size_t nEmissions) const;
};

constexpr double emissionPmf[3][CONTEXT_NUMBER][OUTCOME_NUMBER] = {
    {// matchPmf
        {   0.0490685786,  0.000808342716,  0.000198524496,  0.000148732926,    0.0655943844,   0.00107162279,  0.000983837213,  3.94796477e-05,     0.877135173,    0.0040625605,  0.000708838315,  0.000121097575},
        {  0.00387883926,     0.011507013,  9.13948572e-05,  5.16430582e-05,    0.0153488408,    0.0371426471,     0.000224268,  0.000255072982,    0.0203800638,     0.910442401,  0.000321465913,  0.000280323168},
        { 0.000432233922,  9.44785334e-05,    0.0173299044,   0.00256972779,   0.00197528288,  4.44643626e-05,    0.0724760208,   0.00323612799,  0.000574385807,  0.000255955416,     0.884915124,     0.016011158},
        { 0.000454426932,  0.000112200233,   0.00836788548,    0.0159362152,  9.25335646e-05,  0.000103163774,    0.0349616999,    0.0460292171,  0.000363212547,  0.000433181139,    0.0773685313,     0.815715149},
        {   0.0205388073,  0.000836497056,  0.000588824539,  0.000187100857,    0.0586792371,   0.00096878204,   0.00115400144,  3.87177042e-05,     0.913238232,   0.00211728628,   0.00114695021,   0.00044362853},
        {   0.0438189114,    0.0109384165,  0.000150015065,  9.35020377e-05,    0.0129245291,    0.0318527699,  0.000185897665,  0.000132825435,    0.0168541502,     0.882889389,  5.40348332e-05,  3.33949266e-05},
        { 0.000360861879,  3.84117097e-05,   0.00720852457,   0.00103747404,  0.000535742001,  9.15036164e-05,    0.0333113081,   0.00170295802,  0.000526151133,  0.000498106276,     0.942324413,     0.012307329},
        { 0.000265238546,  0.000176139842,    0.0104070994,    0.0164876561,  8.76448456e-05,  8.81823058e-05,    0.0322729483,    0.0475522835,  0.000346258242,  0.000683284248,    0.0775451031,     0.814004031},
        {  0.00987645414,  0.000418107961,  0.000196775866,  0.000126944678,    0.0326467357,  0.000493721702,  0.000656948821,  2.34422213e-05,      0.95082053,   0.00395960317,  0.000406271368,   0.00030159854},
        {  0.00208206263,   0.00590612695,  0.000198296332,  3.18473899e-05,   0.00736826574,    0.0195638107,  0.000248222378,  0.000153461089,    0.0124171336,     0.951406855,  0.000333588933,  0.000239582648},
        {   0.0374083923,  7.51712737e-05,   0.00665374143,   0.00110622602,  0.000113449438,  3.09671401e-05,    0.0262360916,   0.00135132587,  0.000524666323,  6.74751739e-05,     0.915035331,    0.0113242266},
        { 0.000216715854,  3.58854699e-05,     0.004419479,   0.00921522491,  0.000248910156,  9.00216206e-05,    0.0181401014,    0.0273831797,  0.000267395037,  0.000742503313,     0.051868398,     0.887294927},
        {   0.0192362915,   0.00087686256,   0.00133296974,  0.000166816026,    0.0562744523,  0.000748684538,   0.00155309126,  4.37453256e-05,     0.912324398,    0.0055886634,   0.00144806003,  0.000312863856},
        {  0.00183069899,   0.00629984804,  0.000164784294,  0.000110642093,   0.00702967674,    0.0195211718,  0.000243273161,  0.000154381448,    0.0117671081,     0.952492955,  8.33559045e-05,  0.000228624592},
        { 0.000150538355,  3.66282657e-05,     0.008512658,   0.00114292242,  0.000454485373,  6.03123555e-05,    0.0341609161,    0.0012982219,  0.000972247894,   0.00015899966,     0.940440965,    0.0125498232},
        {   0.0290674502,  0.000699127513,   0.00697586658,    0.0074766871,   0.00324403268,   0.00101436563,    0.0171499416,    0.0228744096,    0.0283644831,    0.0262597707,    0.0761732956,     0.780648571}},


    {// branchPmf
        {    0.342625724,  0.000416360813,  0.000416360813,  0.000416360813,     0.200310897,  0.000416360813,  0.000416360813,  0.000416360813,     0.451234328,  0.000416360813,  0.000416360813,  0.000416360813},
        { 0.000956676671,      0.10706572,  0.000956676671,  0.000956676671,  0.000956676671,    0.0858499434,  0.000956676671,  0.000956676671,  0.000956676671,     0.793690863,  0.000956676671,  0.000956676671},
        {  0.00027270246,   0.00027270246,     0.273881665,   0.00027270246,   0.00027270246,   0.00027270246,     0.163397174,   0.00027270246,   0.00027270246,   0.00027270246,     0.558903327,   0.00027270246},
        { 0.000681543674,  0.000681543674,  0.000681543674,     0.148804662,  0.000681543674,  0.000681543674,  0.000681543674,     0.106754268,  0.000681543674,  0.000681543674,  0.000681543674,     0.734899459},
        {    0.168053431,  0.000118965157,  0.000118965157,  0.000118965157,      0.13198635,  0.000118965157,  0.000118965157,  0.000118965157,     0.698294707,  0.000118965157,  0.000118965157,  0.000118965157},
        { 0.000205108667,    0.0531328883,  0.000173569826,  0.000173569826,  0.000173786628,    0.0457788373,  0.000173569826,  0.000173569826,  0.000186168486,     0.898613943,  0.000173569826,  0.000173569826},
        { 0.000146043562,  0.000146043562,     0.218268977,  0.000146043562,  0.000146043562,  0.000146043562,     0.140286078,  0.000146043562,  0.000146043562,  0.000146043562,     0.639400335,  0.000146043562},
        { 0.000688738719,  0.000688738719,  0.000688738719,       0.1041046,  0.000688738719,  0.000688738719,  0.000688738719,      0.11315224,  0.000688738719,  0.000688738719,  0.000688738719,     0.773100818},
        {    0.124999283,  0.000148121105,  0.000148121105,  0.000148121105,    0.0910614333,  0.000148121105,  0.000148121105,  0.000148121105,     0.781865588,  0.000148121105,  0.000148121105,  0.000148121105},
        {  0.00045773541,    0.0859072817,   0.00045773541,   0.00045773541,   0.00045773541,    0.0820955628,   0.00045773541,   0.00045773541,   0.00045773541,      0.82558886,   0.00045773541,   0.00045773541},
        { 0.000421318581,   0.00042126564,     0.396171182,   0.00042126564,  0.000518297203,   0.00042126564,     0.196822139,   0.00042126564,  0.000426833657,   0.00042126564,     0.401006307,   0.00042126564},
        {  0.00051094738,   0.00051094738,   0.00051094738,    0.0725177578,   0.00051094738,   0.00051094738,   0.00051094738,    0.0619382495,   0.00051094738,   0.00051094738,   0.00051094738,     0.858390729},
        {    0.131136418,  0.000164812265,  0.000164812265,  0.000164812265,     0.120277362,  0.000164812265,  0.000164812265,  0.000164812265,     0.746278849,  0.000164812265,  0.000164812265,  0.000164812265},
        { 0.000329798356,    0.0510585908,  0.000329798356,  0.000329798356,  0.000329798356,    0.0611592526,  0.000329798356,  0.000329798356,  0.000329798356,      0.88316498,  0.000329798356,  0.000329798356},
        { 0.000180865401,  0.000180865401,     0.205582698,  0.000180865401,  0.000180865401,  0.000180865401,     0.123360032,  0.000180865401,  0.000180865401,  0.000180865401,     0.668525155,  0.000180865401},
        { 0.000186658164,   0.00012839928,   0.00012839928,    0.0559315702,  0.000128444271,   0.00012839928,   0.00012839928,    0.0770335504,   0.00030571667,   0.00012839928,   0.00012839928,     0.865001668}},

    {// stickPmf
        { 0.000277685074,    0.0417252045,     0.364922788,    0.0500229872,  0.000277685074,    0.0189918334,     0.136206612,    0.0200429607,  0.000277685074,     0.107998074,     0.170819489,    0.0870485712},
        {    0.133578351,  0.000188315822,     0.308078563,    0.0400178107,     0.061982011,  0.000188315822,     0.140354129,    0.0203735931,     0.107984177,  0.000188315822,     0.107575551,    0.0785492868},
        {    0.336284312,     0.051820026,  0.000587615385,     0.134309535,    0.0629375938,    0.0167340057,  0.000587615385,    0.0738947979,    0.0605851001,    0.0793139492,  0.000587615385,     0.179419757},
        {    0.240026428,    0.0441930049,     0.241465266,   0.00024036228,    0.0820393901,    0.0193107959,     0.150745472,   0.00024036228,    0.0540448043,    0.0524441393,       0.1138078,   0.00024036228},
        { 0.000180797891,    0.0320889304,     0.372050698,    0.0519061217,  0.000180797891,   0.00760625602,     0.142299495,    0.0201454343,  0.000180797891,     0.175549117,     0.132371406,    0.0645361585},
        {    0.114758774,  9.51379064e-05,      0.28133416,    0.0403378649,    0.0750935043,  8.42767139e-05,     0.173947953,    0.0166338642,     0.131259847,  8.61046238e-05,      0.12028389,    0.0456681878},
        {    0.291000124,    0.0407605367,  0.000198597625,    0.0597961964,     0.135790286,    0.0272185085,  0.000198597625,    0.0244154796,     0.146095649,     0.148845059,  0.000198597625,      0.12448938},
        {    0.213420218,    0.0298572908,     0.186745555,  0.000267027493,     0.093752031,    0.0152431981,     0.115286566,  0.000267027493,     0.086018103,     0.141042664,     0.116498155,  0.000267027493},
        {  0.00034456068,    0.0352939245,     0.398266288,    0.0542827844,   0.00034456068,    0.0124175456,     0.141878989,    0.0170548733,   0.00034456068,     0.122289029,    0.0932034947,     0.122556586},
        {    0.114541438,  8.86227971e-05,      0.35179944,     0.043330166,    0.0611624765,  8.86227971e-05,     0.226336937,    0.0178851686,    0.0436941493,  8.86227971e-05,    0.0834659753,    0.0570752663},
        {    0.318416407,    0.0579399913,   0.00143058768,    0.0793018387,      0.13473448,    0.0242078947,  0.000480640884,    0.0204456085,     0.138556017,     0.075744249,  0.000456896993,     0.146078987},
        {    0.289121416,    0.0450738183,     0.271758023,  0.000388554739,     0.102061905,    0.0111931178,     0.122103464,  0.000388554739,    0.0916006041,    0.0365694207,    0.0274097923,  0.000388554739},
        { 0.000347420644,    0.0360243203,     0.312236676,    0.0494711108,  0.000347420644,    0.0136543489,     0.133997713,    0.0198445274,  0.000347420644,    0.0851807391,     0.210300962,     0.136510237},
        {    0.145899101,  0.000166862191,     0.242177432,    0.0387027921,    0.0533086833,  0.000166862191,     0.110673618,    0.0175713693,    0.0677301103,  0.000166862191,     0.124938792,     0.197663205},
        {     0.35446875,    0.0448309155,  0.000364014885,    0.0777168986,     0.141455244,    0.0195149694,  0.000364014885,     0.024359256,     0.118317246,    0.0552441261,  0.000364014885,     0.161180475},
        {    0.162167278,    0.0268508415,     0.209924244,  0.000193471148,    0.0739113275,    0.0102206136,     0.172347554,  0.000163585965,    0.0633846456,     0.051012026,     0.228860686,  0.000167429533}}};

constexpr double transProbs[16][3][4] = {
    // Fit for context:  AA

    {
        { -4.58216653258377, 0.687233863051845, -0.154083259755313, 0.0114321184985665  },
        { -0.566307049216117, -1.17721138761644, 0.183688871696491, -0.0102230250471649  },
        { 1.55964949069417, -1.49408937589455, 0.192409907621455, -0.0107363759180426  }
    },
    // Fit for context:  AC
    {
        { -9.12062708424519, 2.69046791086471, -0.478689729927781, 0.0284461989238249  },
        { -4.52166209104686, 1.3046852217299, -0.293023836249338, 0.0219143007296014  },
        { -5.66904537637777, 1.71920288092646, -0.405844858353435, 0.0281955889895175  }
    },
    // Fit for context:  AG
    {
        { -6.74566659882406, 2.07641443499033, -0.347033216992796, 0.0185629529478924  },
        { -2.37768657840025, -0.370906490876627, 0.0553441488513224, -0.00505838466319259  },
        { 1.0268062000401, -0.992101482052789, 0.0232410961826395, 0.00498641246117365  }
    },
    // Fit for context:  AT
    {
        { 0.805566650997283, -2.37447906504741, 0.403885965682443, -0.0235391014707725  },
        { -4.39071710204575, 0.62852368286356, -0.0925318542625281, 0.00473923421245102  },
        { -1.74378290827118, -0.588623249302254, 0.0407041596645099, 0.000782287562388532  }
    },
    // Fit for context:  CA
    {
        { -4.78597212823989, 0.880172219607289, -0.0763870449781868, 0.000226154498806874  },
        { -1.81182754179299, -0.780917422749851, 0.183174657962876, -0.0127220800289381  },
        { 0.696491929495958, -1.06156216334184, 0.0860837518209075, -0.00249776503379324  }
    },
    // Fit for context:  CC
    {
        { 13.1313671779978, -8.18494160805398, 1.40628613179, -0.0793210826767796  },
        { 4.32431038347223, -3.18803372082616, 0.532761343838732, -0.0281471451467225  },
        { 2.1469627367621, -3.09095113181557, 0.628481119888392, -0.0420018038375466  }
    },
    // Fit for context:  CG
    {
        { -1.88963050600106, -0.772602331471341, 0.212056991260888, -0.0167011968658315  },
        { -7.94622958292638, 2.76237097279153, -0.492857837051516, 0.0287801869740674  },
        { 0.731637618821747, -1.38514974486228, 0.118725630338243, -0.00123153166482531  }
    },
    // Fit for context:  CT
    {
        { 2.41969457412114, -3.16833322038735, 0.577279266997955, -0.0370366438782084  },
        { -13.822631635269, 5.88533305065379, -1.03018477573883, 0.0592515975151891  },
        { -0.875545738208091, -0.724639160525869, 0.0234867193924595, 0.00470715553387549  }
    },
    // Fit for context:  GA
    {
        { 17.6186929207558, -11.3441944513817, 2.12226205764468, -0.130541424460828  },
        { -0.21892757000539, -1.62557811545721, 0.300175180569162, -0.0186217436024484  },
        { -5.35466156815154, 2.16708096630271, -0.549398837732772, 0.0383644669980167  }
    },
    // Fit for context:  GC
    {
        { -4.00540299833362, 0.129394070344217, -0.00215612840246416, -0.00248970697890298  },
        { 0.0354053301640659, -1.80685511508197, 0.400864938033335, -0.0262694990357512  },
        { 1.46328423271853, -2.46601779445692, 0.36351693567861, -0.0174310397532786  }
    },
    // Fit for context:  GG
    {
        { -7.37089873903281, 2.12287596486111, -0.36120094766727, 0.0197531430986257  },
        { 2.95095337934729, -3.19420253691392, 0.548805416066258, -0.0326191904660653  },
        { -1.93719686965216, -0.00260678228586515, -0.035271165135989, 0.00213311742302812  }
    },
    // Fit for context:  GT
    {
        { -20.9089083023533, 9.49214996072597, -1.66413304931095, 0.0937389668515694  },
        { -3.51723305331226, -0.241565020692246, 0.111473889995382, -0.0103914893152198  },
        { -5.96713960173248, 1.49157732856665, -0.291788891790897, 0.0179501023886058  }
    },
    // Fit for context:  TA
    {
        { -9.65959242262948, 4.02257794513185, -0.700108169437621, 0.03935010943199  },
        { -8.6043505311291, 2.26933168858594, -0.259692578874624, 0.00660330430930349  },
        { 6.25352569939231, -3.85441465393031, 0.532652512502431, -0.0256805818192053  }
    },
    // Fit for context:  TC
    {
        { 0.878783463347493, -2.17490864187508, 0.385324390474283, -0.0221380876241637  },
        { -7.30005117647813, 2.11778474057213, -0.289593226559052, 0.0119514469416285  },
        { -5.0160011902027, 0.387297563790036, -0.0504636293118311, 0.00311089059364794  }
    },
    // Fit for context:  TG
    {
        { -7.60555673192205, 2.18668042183426, -0.299901640988149, 0.0121413300238935  },
        { 1.80667391148329, -3.17590248876971, 0.646433064745285, -0.0437772973369292  },
        { 13.7747092945745, -8.71434608350849, 1.48853785404397, -0.0859038022870709  }
    },
    // Fit for context:  TT
    {
        { -8.70473590657829, 3.23086795205131, -0.548184540910583, 0.0308453661016667  },
        { 3.23182615565406, -3.0800614694476, 0.533215728487221, -0.0303679753377579  },
        { 4.57212258109004, -4.88550946707369, 1.0133460021107, -0.0693366128266997  }
    }
  };


inline double CalculateExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t row, const bool secondMoment)  {
    double expectedLL = 0;
    for(size_t i = 0; i < OUTCOME_NUMBER; i++) {
        double curProb = emissionPmf[index][row][i];
        double lgCurProb = std::log(curProb);
        if(!secondMoment) {
            expectedLL +=  curProb * lgCurProb;
        } else {
            expectedLL += curProb * pow(lgCurProb, 2.0);
        }
    }
    return expectedLL;
}

double SP1C1PwModel::ExpectedLogLikelihoodOfOutcomeRow(const int index, const uint8_t prev, const uint8_t curr, const bool secondMoment) const {
    const auto row = (prev << 2) | curr;
    const auto moment = secondMoment ? 1 : 0;
    return cachedEmissionExpectations_[row][index][moment];
}

SP1C1PwModel::SP1C1PwModel(const SNR& snr) : snr_(snr)
{
    // Generate cached transistion probabilities
    const double snr1 = snr_.A, snr2 = snr1 * snr1, snr3 = snr2 * snr1;
    for (int ctx = 0; ctx < CONTEXT_NUMBER; ++ctx) {
        double sum = 1.0;
        ctxTrans_[ctx][0] = 1.0;
        for (size_t j = 0; j < 3; ++j) {
            double xb = transProbs[ctx][j][0] + snr1 * transProbs[ctx][j][1] +
                        snr2 * transProbs[ctx][j][2] + snr3 * transProbs[ctx][j][3];
            xb = std::exp(xb);
            ctxTrans_[ctx][j + 1] = xb;
            sum += xb;
        }
        for (size_t j = 0; j < 4; ++j)
            ctxTrans_[ctx][j] /= sum;
    }

    // Generate cached emission expectations
    // TODO: These are identical for all instances, either we should enrich the model or avoid doing this in a context dependent way
    for(int ctx = 0; ctx < CONTEXT_NUMBER; ctx++) {
        for (int index = 0; index < 3; index++) {
            cachedEmissionExpectations_[ctx][index][0] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, false);
            cachedEmissionExpectations_[ctx][index][1] = CalculateExpectedLogLikelihoodOfOutcomeRow(index, ctx, true);
        }
    }
}

std::unique_ptr<AbstractRecursor> SP1C1PwModel::CreateRecursor(
    std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr, double scoreDiff) const
{
    return std::unique_ptr<AbstractRecursor>(
        new SP1C1PwRecursor(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff));
}

std::vector<TemplatePosition> SP1C1PwModel::Populate(const std::string& tpl) const
{
    std::vector<TemplatePosition> result;

    if (tpl.empty()) return result;

    result.reserve(tpl.size());

    // Calculate probabilities in all 16 Contexts
    uint8_t prev = detail::TranslationTable[static_cast<uint8_t>(tpl[0])];
    if (prev > 3) throw std::invalid_argument("invalid character in template!");

    for (size_t i = 1; i < tpl.size(); ++i) {
        const uint8_t curr = detail::TranslationTable[static_cast<uint8_t>(tpl[i])];
        if (curr > 3) throw std::invalid_argument("invalid character in template!");
        const auto row = (prev << 2) | curr;
        const auto params = ctxTrans_[row];
        result.emplace_back(TemplatePosition{
            tpl[i - 1], prev,
            params[0],  // match
            params[1],  // branch
            params[2],  // stick
            params[3]   // deletion
        });
        prev = curr;
    }
    result.emplace_back(TemplatePosition{tpl.back(), prev, 1.0, 0.0, 0.0, 0.0});

    return result;
}



double SP1C1PwModel::ExpectedLogLikelihoodForMatchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::MATCH), prev, curr, secondMoment);
}
double SP1C1PwModel::ExpectedLogLikelihoodForStickEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::STICK), prev, curr, secondMoment);
}
double SP1C1PwModel::ExpectedLogLikelihoodForBranchEmission(uint8_t prev, uint8_t curr, bool secondMoment) const {
    return ExpectedLogLikelihoodOfOutcomeRow(static_cast<uint8_t>(MoveType::BRANCH), prev, curr, secondMoment);
}

SP1C1PwRecursor::SP1C1PwRecursor(std::unique_ptr<AbstractTemplate>&& tpl, const MappedRead& mr,
                                 double scoreDiff)
    : Recursor<SP1C1PwRecursor>(std::forward<std::unique_ptr<AbstractTemplate>>(tpl), mr, scoreDiff)
{
}

std::vector<uint8_t> SP1C1PwRecursor::EncodeRead(const MappedRead& read)
{
    std::vector<uint8_t> result;

    result.reserve(read.Length());

    for (size_t i = 0; i < read.Length(); ++i) {
        if (read.PulseWidth[i] < 1U) throw std::runtime_error("invalid PulseWidth in read!" + std::to_string(read.PulseWidth[i]));
        const uint8_t pw = std::min(2, read.PulseWidth[i] - 1);
        const uint8_t bp = detail::TranslationTable[static_cast<uint8_t>(read.Seq[i])];
        if (bp > 3) throw std::invalid_argument("invalid character in read!");
        const uint8_t em = (pw << 2) | bp;
        if (em > 11) throw std::runtime_error("read encoding error!");
        result.emplace_back(em);
    }

    return result;
}

double SP1C1PwRecursor::EmissionPr(MoveType move, uint8_t emission, uint8_t prev, uint8_t curr)
{
    assert(move != MoveType::DELETION);
    const auto row = (prev << 2) | curr;
    return emissionPmf[static_cast<uint8_t>(move)][row][emission] * kCounterWeight;
}

double SP1C1PwRecursor::UndoCounterWeights(const size_t nEmissions) const
{
    return -std::log(kCounterWeight) * nEmissions;
}
}  // namespace anonymous
}  // namespace Consensus
}  // namespace PacBio
