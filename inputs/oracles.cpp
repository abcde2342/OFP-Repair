#include <mpreal.h>
#include <gsl/gsl_cdf.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using mpfr::mpreal;
using boost::multiprecision::cpp_dec_float_50;

const mpfr_prec_t highprec=512;

mpreal oracle_erfcinv(mpreal q) {
    if (q <= 0) {
        return std::numeric_limits<double>::infinity();
    }
    if (q >= 2.0) {
        return -std::numeric_limits<double>::infinity();
    }

    mpreal lo = 0.0;
    mpreal hi = 4.0;

    while (mpfr::erfc(hi) > q) {
        hi *= 2.0;
    }

    for (int i = 0; i < 96; ++i) {
        mpreal mid = (lo + hi) / 2.0;
        if (mpfr::erfc(mid) > q) {
            lo = mid;
        }
        else {
            hi = mid;
        }
    }

    return (lo + hi) / 2.0;
}

mpreal oracle_cos_x2(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = 0.5 + x*x*(-1.0/24.0 + x*x*(1.0/720.0));
    }
    else {
        y = (1-cos(x))/(x*x);
    }
    return y;
}

mpreal oracle_exp_2(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(1.0 + x*x*(1.0/12.0 + x*x*(1.0/360.0)));
    }
    else {
        y = ((exp(x)-2)+exp(-x))/x;
    }
    return y;
}

mpreal oracle_cos_sin(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(0.5 + x*x*(1.0/24.0 + x*x*(1.0/240.0)));
    }
    else {
        y = (1-cos(x))/sin(x);
    }
    return y;
}

mpreal oracle_sin_sin(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal eps = 1e-6;
    y = sin(x+eps) - sin(x);
    return y;
}

mpreal oracle_tan_tan(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal eps = 1e-6;
    y = tan(x+eps) - tan(x);
    return y;
}

mpreal oracle_cos_cos(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal eps = 1e-6;
    y = cos(x+eps) - cos(x);
    return y;
}

mpreal oracle_exp_exp(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    y = 2 * sinh(x);
    return y;
}

mpreal oracle_exp_1(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(1.0 + x*(0.5 + x*(1.0/6.0 + x*(1.0/24.0))));
    }
    else {
        y = exp(x)-1;
    }
    return y;
}

mpreal oracle_x_tan(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(1.0/3.0 + x*x*(1.0/45.0 + x*x*(2.0/945.0)));
    }
    else {
        y = 1/x - 1/tan(x);
    }
    return y;
}

mpreal oracle_log_log(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-4) {
        y = -1.0 + x*(-1.0 + x*(-0.5 + x*(-5.0/12.0 + x*(-7.0/24.0 + x*(-191.0/720.0)))));
    }
    else {
        y = log(1-x) / log(1+x);
    }
    return y;
}

mpreal oracle_log_x(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(-2.0 + x*x*(-2.0/3.0 + x*x*(-2.0/5.0)));
    }
    else {
        y = log((1-x)/(1+x));
    }
    return y;
}

mpreal oracle_sqrt_exp(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        mpreal sqrt2 = sqrt((mpreal)2);
        y = sqrt2 + x*(sqrt2/4.0 + x*(3.0*sqrt2/32.0 + x*(7.0*sqrt2/384.0)));
    }
    else {
        y = sqrt((exp(2*x)-1) / (exp(x)-1));
    }
    return y;
}

mpreal oracle_sin_tan(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = -0.5 + x*x*(9.0/40.0 + x*x*(-27.0/2800.0));
    }
    else {
        y = (x-sin(x))/(x-tan(x));
    }
    return y;
}

mpreal oracle_exp_x(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = 1.0 + x*(0.5 + x*(1.0/6.0 + x*(1.0/24.0 + x*(1.0/120.0))));
    }
    else {
        y = ((exp(x) - 1) / x);
    }
    return y;
}

mpreal oracle_x_x2(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = 1.0 + x*(-1.0 + x*(1.0 + x*(-1.0 + x*(1.0 + x*(-1.0 + x)))));
    }
    else {
        y = 1 / (x + 1);
    }
    return y;
}

mpreal oracle_exp_BI(double dx) {
    // https://www.wolframalpha.com/input/?i=series+%28exp%28besselI%281%2Cx%29%29-1.0%29%2Fx+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    int num = 11;
    mpreal c[11] = {
        0.5,
        ((mpreal)1.0) / 8.0,
        ((mpreal)1.0) / 12.0,
        ((mpreal)13.0) / 384.0,
        ((mpreal)41.0) / 3840.0,
        ((mpreal)211.0) / 46080.0,
        ((mpreal)109.0) / 71680.0,
        ((mpreal)5209.0) / mpreal("10321920"),
        ((mpreal)33139.0) / mpreal("185794560"),
        ((mpreal)203641.0) / mpreal("3715891200"),
        ((mpreal)725819.0) / mpreal("40874803200"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_bJ_sin(double dx) {
    // https://www.wolframalpha.com/input/?i=series+%281-besselj%280%2Cx%29%29%2Fsin%28x%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = x*(0.25 + x*x*(5.0/192 + x*x*(31.0/11520.0)));
    }
    else {
        y = (1 - mpfr::besselj0(x)) / sin(x);
    }
    return y;
}

mpreal oracle_di_tan(double dx) {
    // https://www.wolframalpha.com/input/?i=series+1%2Fdilog%28x%29+-+1%2Ftan%28x%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        y = -0.25 + x*(41.0/144.0 + x*(-13.0/576.0 + x*(4609.0/518400.0 + x*(-6151.0/691200.0))));
    }
    else {
        y = (1 / mpfr::li2(x)) - (1 / tan(x));
    }
    return y;
}

mpreal oracle_log_erf(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    y = log1p(-mpfr::erf(x)) / log1p(x);
    return y;
}

mpreal oracle_acos_fd(double dx) {
    // https://www.wolframalpha.com/input/?i=series+acos%28x%29*acos%28x%29+%2B+3*polylog%282%2C+-exp%28x%29%29+at+x%3D0
    // dilog == li2 == polylog(2, x)
    // fermi_dirac_1 == -polylog(2, -exp(x))
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    if (fabs(x) < 1e-6) {
        mpreal c1 = -mpfr::const_pi()-mpfr::log(8);
        mpreal c2 = 0.25;
        mpreal c3 = -1.0/8.0 - mpfr::const_pi()/6.0;
        y = (c1 + x*(c2 + x*c3));
    }
    else {
        y = (acos(x)*acos(x) + 3*mpfr::li2(-exp(x))) / x;
    }
    return y;
}

mpreal oracle_sci3545_2(double dx) {
    mpreal::set_default_prec(highprec);
    boost::math::normal_distribution<cpp_dec_float_50> norm;
    cpp_dec_float_50 q(dx);
    cpp_dec_float_50 y = boost::math::quantile(boost::math::complement(norm, q / 2));
    return mpreal(y.str(60));
}

mpreal oracle_sci4034_1(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    y = 1-exp(-0.5*x*x);
    return y;
}

mpreal oracle_sci4034_2(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    y = sqrt(-2*log(1-x));
    return y;
}

mpreal oracle_sci3547(double dx) {
    mpreal::set_default_prec(highprec);
    cpp_dec_float_50 y(dx);
    cpp_dec_float_50 val = boost::math::erfc_inv(y);
    return mpreal(val.str(60));
}

mpreal oracle_sci3545_1(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx;
    mpreal y = mpfr::erfc(sqrt(0.5 / x));
    return y;
}

mpreal oracle_ei(double dx) {
    // https://www.wolframalpha.com/input?i=series+sin%28erf%5E%28-1%29%28x%29%29%2F%28cos%28x%29-e%5Ex%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal pi = mpfr::const_pi();
    mpreal spi = sqrt(pi);
    int num = 10;
    mpreal c[10] = {
        - spi/2.0,
        spi/2.0,
        -spi*(20.0+pi)/48.0,
        spi*(16.0+pi)/48.0,
        -spi*(2992.0+200.0*pi+27.0*pi*pi)/11520.0,
        spi*(2320.0+160.0*pi+27.0*pi*pi)/11520.0,
        -spi*(301760.0+20944.0*pi+3789.0*pi*pi+651.0*pi*pi*pi)/mpreal("1935360"),
        spi*(233472.0+16240.0*pi+3024.0*pi*pi+651.0*pi*pi*pi)/mpreal("1935360"),
        -spi*(mpreal("86697216")+mpreal("6035200")*pi+mpreal("1130976")*pi*pi
            +mpreal("260400")*pi*pi*pi+47925.0*pi*pi*pi*pi)/mpreal("928972800"),
        spi*(mpreal("13413632")+mpreal("933888")*pi+mpreal("175392")*pi*pi
            +41664.0*pi*pi*pi+9585.0*pi*pi*pi*pi)/mpreal("928972800"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_Q1_W(double dx) {
    // https://www.wolframalpha.com/input?i=series+%281%2BLegendreQ%5B1%2Cx%5D%29%2F%28W%28x%29%5E2%29+at+x+%3D+0
    // https://math.stackexchange.com/questions/1700919/how-to-derive-the-lambert-w-function-series-expansion
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    int num = 12;
    mpreal c[12] = {
        1.0,
        2.0,
        ((mpreal)1.0) / 3.0,
        1.0,
        -((mpreal)7.0) / 15.0,
        ((mpreal)67.0) / 36.0,
        -((mpreal)307.0) / 105.0,
        ((mpreal)2521.0) / 360.0,
        -((mpreal)14039.0) / 945.0,
        mpreal("31188749") / mpreal("907200"),
        -mpreal("493817") / mpreal("6237"),
        mpreal("1017580243") / mpreal("5443200"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

// mpreal oracle27(double dx) {
//     // https://www.wolframalpha.com/input?i=series+x%5E2*trigamma%28x%29%2Fsqrt%283%29+at+x+%3D+0
//     // https://mathworld.wolfram.com/AiryFunctions.html
//     mpreal::set_default_prec(highprec);
//     mpreal x = dx, y;
//     mpreal pi = mpfr::const_pi();
//     mpreal c3 = 1.0/sqrt(mpreal("3"));
//     int num = 10;
//     mpreal c[10] = {
//         c3,
//         0,
//         pow(pi, 2) * c3 / 6.0,
//         -2.0 * zeta(mpreal("3")) * c3,
//         pow(pi, 4) * c3 / 30.0,
//         -24.0 * zeta(mpreal("5")) * c3 / 6.0,
//         pow(pi, 6) * c3 / 189.0,
//         -720.0 * zeta(mpreal("7")) * c3 / 120.0,
//         pow(pi, 8) * c3 / 1350.0,
//         -40320 * zeta(mpreal("9")) * c3 / 5040.0,
//     };
//     y = 0;
//     for (int i = num-1; i >= 0; i--) {
//         y = y*x + c[i];
//     }
//     return y;
// }

mpreal oracle_bj_tan(double dx) {
    // https://www.wolframalpha.com/input?i=series+%281-sin%28x%29%2Fx%29%2F%28x*tan%28x%29%29+at+x+%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    int num = 11;
    mpreal c[11] = {
        mpreal("1") / 6.0,
        0,
        - mpreal("23") / 360.0,
        0,
        - mpreal("11") / 15120.0,
        0,
        - 143.0 / mpreal("604800"),
        0,
        - 361.0 / mpreal("17107200"),
        0,
        - mpreal("1416533") / mpreal("65383718400"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_Si_tan(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28Si%28x%29-tan%28x%29%29%2F%28x%5E3%29+at+x+%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    int num = 12;
    mpreal c[12] = {
        - mpreal("7") / 18.0,
        0,
        - mpreal("79") / 600.0,
        0,
        - mpreal("127") / mpreal("2352"),
        0,
        - mpreal("71423") / mpreal("3265920"),
        0,
        - mpreal("555959") / mpreal("62726400"),
        0,
        - mpreal("3589967") / mpreal("999398400"),
        0,
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_by_psi(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28cos%28x%29%2Fx%29%5E2+-+trigamma%28x%29+at+x+%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    static mpreal pi = mpfr::const_pi();
    static mpreal pi2 = pi*pi;
    static mpreal pi4 = pi2*pi2;
    static mpreal pi6 = pi4*pi2;
    static mpreal pi8 = pi6*pi2;
    static mpreal pi10 = pi8*pi2;
    static mpreal pi12 = pi10*pi2;
    int num = 11;
    static mpreal c[11] = {
        -1.0-pi2/6.0,
        2.0 * zeta(mpreal("3")),
        mpreal("1")/3.0-pi4/30.0,
        4.0 * zeta(mpreal("5")), // 24.0 * zeta(mpreal("5")) / 6.0
        (-42.0 - 5*pi6) / mpreal("945"),
        6.0 * zeta(mpreal("7")), // 720.0 * zeta(mpreal("7")) / 120.0
        (30 - 7*pi8) / mpreal("9450"),
        8.0 * zeta(mpreal("9")),
        (-22-15*pi10) / mpreal("155925"),
        10.0 * zeta(mpreal("11")),
        (2730 - 7601*pi12) / mpreal("638512875"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_fdm_log(double dx) {
    // https://www.wolframalpha.com/input?i=series+%282*e%5Ex%2F%281%2Be%5Ex%29-1%29%2Fln%281%2Bx%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    int num = 10;
    mpreal c[10] = {
        mpreal("1") / 2.0,
        mpreal("1") / 4.0,
        -mpreal("1") / 12.0,
        0,
        -mpreal("1") / 180.0,
        mpreal("7") / 720.0,
        -mpreal("823") / mpreal("120960"),
        mpreal("1117") / mpreal("241920"),
        -mpreal("3319") / mpreal("806400"),
        mpreal("10319") / mpreal("2903040"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_eQ_sqrt(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28erfc%28-x%2Fsqrt%282%29%29-%281%2Bx%29%5E%281%2F2%29%29%2Fx+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal pi = mpfr::const_pi();
    mpreal s2pi = sqrt(2*pi);
    int num = 10;
    mpreal c[10] = {
        - 2.0 / s2pi - mpreal("1") / 2.0,
        mpreal("1") / 8.0,
        - mpreal("1") / 16.0 + mpreal("1") / s2pi / 3.0,
        mpreal("5") / 128.0,
        - mpreal("1") / 20.0 / s2pi - mpreal("7") / 256.0,
        mpreal("21") / 1024.0,
        - mpreal("33") / 2048 + mpreal("1") / s2pi / 168.0,
        mpreal("429") / mpreal("32768"),
        - mpreal("1") / 1728.0 / s2pi - mpreal("715") / mpreal("65536"),
        mpreal("2431") / mpreal("262144"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}


// For oracle lambert
namespace {
mpreal lambert_hornor_eval(mpreal x) {
    static const mpreal coefs[] = {
        mpreal("-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00", 2048),
        mpreal("2.33164398159712420336353606216840087638023629918758842300809644777601004941265734950262999179547775424833502631207336244836642977791061698455539548600555307416996882138541742694334174087392312213363683987122700179387442559506512788634356927720937116053557860128366386599733164485065199170625213160832004835937544209214383692602651217262355792696583043564578955225266808603204078284291653000434342979943549056663121370749460830774420660582011918347651141120384933042095664492394514140021279071978771163128207781183132515286382897991170054223044715148811055913935631930685228833658607045692545630494092658006995786287848e+00", 2048),
        mpreal("-1.81218788563936349024019164756844166517149806246663971664464508514938442023569839638092145235011095161831092795466870661454494239775269571526688953017373042049208821908529566050882255325383546350067934104922279195380143605943328992277833949650764044538817653344565160790249489696162473835938516332804636780184122573750755425638866716802995588437353173782474213380621913941829583136482046465147287342779522454601701007243830918140750159318961670463579784719029997997863124296993732528775792615339919541849078547694999948638423432147217993012911202212168591293233097673721395949321965919554690833962874486749213137894414e+00", 2048),
        mpreal("1.93663111449235975536327745766838306382688574831514085485065868872121244237192188993413191870651779067641135602483219996872773343984420346505017610766600323217330040245395116560470356951516806792745605297072913735096298141242687170013842782992659650681925952581473461007533688529937027584289896577827568429063983706201044845897984259131670644534615221107244864704234245223248512468061306334024495828344546381022918064196017739238372257443943159196591190160928360958373820628336786140368579946576975265415663124077201323998773335240782226643125203013848616754424344900453390376065333059756626237177694788222039856362684e+00", 2048),
        mpreal("-2.35355120188161451682154356151648397012410051506466248100552960272941239430666284315204027302520175685386211803270808328390398487992292411187089529760086921719040951489717527710520612763514577281276309830811006172563986297953314888634759917985800916229304645323735423662097429399149177362432665216461599107562532179167481693542983090492129132370910279934295273961357150014269641624385903304112182430363151870950838469859031730220182801990521346200619466953582271591323234983773900455470378248800623669154819847012683145358247759026756370970650548229630690940422557105693696943290775841097392844563903571618634430608796e+00", 2048),
        mpreal("3.06685890105063191289314892270400749848838989930150683701410272579127772496366295229313832650525106780114670786382941818040740537439584900408229711589535326238617854705169763433519873562798156383723336291702490104975311428502671414095150356147233721376397000445683843872349670761107282666250861459463093780624705360134232165898026433652937745709527648348954169124083834197126596892111060730219813663664792625308241628977954864941408436141000024875277427643712614649008038298647454077309628306802079583816184283010914772903123949697496338390776307720315770015394777388258110113581672416068367670739714495711372720767624e+00", 2048),
        mpreal("-4.17533560025817713885498417746037357341265385756187389244702595157055373456397525652787195524130926483334253554688369630747111934902962960944518228276963890546454522197728696555014943640912296371268826353511292621178495848731504686570896075360486062395332468714717673501391834808005708172060343267145200386488513311924003227022939217198321432227423225443240557794931930056017334198597690524979525282480894611239838117299065853274668017531920538873514014433406705165071839766679487032981496263291955523684020704021875358584979774795181722319436946711965146440604148209498249673962737789788878155973580517941929511118303e+00", 2048),
        mpreal("5.85802372987477414881505384611862130415924138010788632736945848782280919514977947533077434042221870411701476126548376249464163680280526653153972768622398457309567747401634463163830092366141157899504606674471794027855319699364392038915268904372798001623948225785868095737154740177735229804140500707406776044624247113776163906079376131621497366149042454170938729022527876693024089871706928139546204612696247658825666638709488665168630399178950476122855248927980367345141767080851288221454461389760946558984719845825703317760797696994521698846629420927358352115422306738237280692974323081314839719646448378214974185872206e+00", 2048),
        mpreal("-8.40103221752397737098416168851388628686484160743087728072095660910373025296289334058647171021732417184285468313138503788316513384571020253037175700290278009353154154816709068235753171889341034421076312914014420703484980598534726034425773530465357340955913304048942084992612978344405594669263267368373439605141647576267460127792804132177289993297243260832935320581414296215967157360558363468045054925889597145476564894818932004618727436395976117439521844028449936910637537944256492716229941192653190045360817681850123699681482211909666741632992658727484357244394612154607995686224853834931124115076096728704446302259583e+00", 2048),
        mpreal("1.22507535013144604237679393600547319876095560250303619518555639177615265162130600075710392451753832164463736074656059001095202683967542641296791780319876560327900683147984241309695657251819192905592975382249231902151417327548650753383668986869841796762052465218339326358009492817678150829868806844849787873194616883006669001296547543961963779758690896481002568395248498226079768689138911453844788015294016240307338930036626817732876323641435157906119168977824822514222571323670297686920863264036885263534184086187202702928145291413786751164904011781135117188200535038542586054993453424591144079570872314885150526582313e+01", 2048),
        mpreal("-1.81006970124724427553771640388630427263597281318579945556499648817756033221685939101164608776582589651357580796955500560438008979827495854386956915769386793136010856264259850406893284476462236507378038306649887117078512281583263829151423022443700012966523599370382117884500090675078054715389907743304549553344822288177997183921877609265589861919292345933312138855438157011171696678609307567799386814769121366037676242115401293872201005880838630904908352560966121816401972942960683138043806407350648770107997626449133898727445792907214448589068225363064154302669103382570785948397994749703377620720656261047592151577357e+01", 2048),
        mpreal("2.70290447990105616503114822804462560520548363497890418088830747911492797291289333459820541189198794299738357927364358056446357737955796830466280745243483523473250636685998022704960564924572931321129342713814161324918741517264247174334455290121713698471128604629017920615113030588165039308794263892307657533638029948238239522834839553117194644861076306239192801608593048460684607125648025708891810439984636665156714262898416850613181139136412044608183032179626057464481245323925077367296978646680958511445619833288141628401808379297419364549753304135035333876695931328838508310480486956897217665511548254854436594692223e+01", 2048),
        mpreal("-4.07154628082606272861347534464736630000701453696083761746349387556322483682431126741279211078343256397950043529318315691230500957432671133346577495868444549723893477665443311310633622591766291830302049246866462307832042307608908846888116971197515000823829657331693588134374019860715153257503019888699625216866407948916073431674015146614579202533033942733159278851995933405426132651402705549220639680248641729628616766801913037109152171255011524513058450207281384886392533964621738908746187450512003697355459843453893966958237625462937563163733346581738145626846313895590723679779643611556264206940282759159995717566086e+01", 2048),
        mpreal("6.17828461870965257412320781663174545149581307721215735418562179522164861587759629108128780480708885569428064826543200620703364107252865701333653095965203722419230582583593896104330317223288106676857143981745573759995916335721937374743355552140207179405121636127591236446017640122963458830283423164092903758959075268110430718263675560498990000844665558627150158439174032235507421157207037199392809413375905958347602739805816706005788610640789699136045060706732104294146306530472504397570540419337536723013597968670141243441745173565129858231733626393034691062380181244774903982958764083363286531379356542787942592505313e+01", 2048),
        mpreal("-9.43366488618669339618453248418100591586049618543683238361589286699249789896178756502693583101753229062354498448184106061411346486892770351680806100607739346485292744413559789283357851470124337945787904302015161041929298661518309203696406831598472400367307376876265091797632726536698918673419150873171510205148343701518970066049676161928685870101769578378666983143142069824338729550178820212572781696516361791798612682155454565846279094640781141281493624107710956759111210497435684664874119605180398341822928777255779943245088503481139780885470089502340863341272137476860274626073032512747413803420266820132457951893037e+01", 2048),
        mpreal("1.44817290387311640026521521190073292596412271946795357505672334647644415592222781642866880985859284035473332259498271743940354822389093851434429555064221513618064776992651543681287493275525323828975781780837888582897532193409532392506308513350773601245611250162900301381510894237469362482976530578867508823608432378454173896278183617001353425546108994604375805824983007487693353875173817549595816883613211519552814239689631548695080838561833899558369182481118839684030753367496155728764418555372101494760991471284213869898670386759598287174594035191708916870301406236316055732211849661666203942289899105821321347518322e+02", 2048),
        mpreal("-2.23349378734845379303554030728851656460915626386153027799902829195548736345380327435531623297053265526838169448253572777347529417050119381179696123207764235258772842758842495132952251711026914173369339333756984882547982624314902413598320083225151597608943727746367059827007576627792484120869688758772860323527035250601219487401779677558976596746523678237899271432334667592376935882269414685796735156832429952438580043235777415389862770575435927025452303308673387928087647013552344077225763497789619553535529711344538614539687105514297545453827008572962378623465760149130998485202744347992719566038548938929397384370744e+02", 2048),
        mpreal("3.45880352864508459383543932858958879008910062234200243883959777217506206203520440357423392990794794293294942924461187089151557984694694676904002561453780848121735439897913461785355553056526723473089960375535721287278804621595838300404082700373245571266491216525196184066441845229403877258844093870803211036157664306772199158577473172305457699413430821197299126606176916106694343030565684839772460061318073695106867553773475714901208490033479020261806257317911218539542322763003022942374866863334580936175580114669719371644485474688468879863511802110963372486571105337453130790447012274714978737704773243046945624072954e+02", 2048),
    };
    int terms = 18;
    mpreal res = 0;
    for (int i = terms-1; i >= 0; i--) {
        res = res * x + coefs[i];
    }
    return res;
}

mpreal lambert_halley_iteration(mpreal x, mpreal w) {
    int max_iters = 64;
    mpreal tol_factor = 1e-50;
    for (int i = 0; i < max_iters; i++) {
        mpreal tol;
        mpreal e = exp(w);
        mpreal p = w + 1.0;
        mpreal t = w*e - x;

        // printf("w, t: %20.16g  %20.16g\n", (double)w, (double)t);

        if (w > 0) {
            t = (t/p)/e; /* Newton iteration */
        }
        else {
            t /= e*p - 0.5*(p + 1.0)*t/p; /* Halley iteration */
        }

        w -= t;

        mpreal temp = (1.0/(fabs(p)*e));
        tol = 10 * tol_factor * max(fabs(w), temp);

        if (fabs(t) < tol) {
            return w;
        }
    }
    // Should not reach here
    return w;
}

mpreal our_lambert(mpreal x) {
    mpreal::set_default_prec(highprec);
    mpreal one_over_e = 1.0/exp((mpreal)1.0);
    mpreal comp_term = 0;
    mpreal q = x + one_over_e - comp_term;

    if (x == 0) {
        return 0.0;
    }
    else if (q <= 0.0) {
        return -1.0;
    }
    else if (q < 1.0e-3) {
        mpreal r = sqrt(q);
        mpreal y = lambert_hornor_eval(r);
        return y;
    }
    else {
        mpreal w;
        if (x < 1.0) {
            mpreal p = sqrt(2.0 * exp((mpreal)1.0) * q);
            w = -1.0 + p*(1.0 + p*(-1.0/3.0 + p*11.0/72.0));
        }
        else {
            w = log(x);
            if (x > 3.0)
                w -= log(w);
        }
        return lambert_halley_iteration(x, w);
    }
}
}// namespace lambert


mpreal oracle_W_var(double dx) {
    // Lambert at 0.
    // mpreal::set_default_prec(highprec);
    // mpreal x = dx, y;
    // int num = 11;
    // mpreal c[11] = {
    //      0.0,
    //      1.0,
    //     -1.0,
    //      1.5,
    //     -mpreal("8") / 3.0,
    //      mpreal("125") / mpreal("24"),
    //     -mpreal("54") / mpreal("5"),
    //      mpreal("16807") / mpreal("720"),
    //     -mpreal("16384") / mpreal("315"),
    //      mpreal("531441") / mpreal("4480"),
    //     -mpreal("156250") / mpreal("567")
    // };
    // y = 0;
    // for (int i = num-1; i >= 0; i--) {
    //     y = y*x + c[i];
    // }

    // Lambert at e.
    // boost::multiprecision::cpp_dec_float_100 x = dx, lx, y;
    // lx = boost::math::lambert_w0(x);

    mpreal x = dx, lx, y;
    lx = our_lambert(x);
    y = 1 / (lx + 1);
    return y; // should be sufficient.
}

mpreal oracle_W_log(double dx) {
    // Lambert at e.
    // boost::multiprecision::cpp_dec_float_100 x = dx, lx, y;
    // lx = boost::math::lambert_w0(x);

    mpreal x = dx, lx, y;
    lx = our_lambert(x);
    y = (lx - 1) / (lx * log(x) - 1);
    return y;
}

mpreal oracle_pow_df(double dx) {
    // https://www.wolframalpha.com/input?i=series+%281%2BF%28x%29%29%5E%281%2Fx%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal e = mpfr::exp(1);
    int num = 10;
    mpreal c[10] = {
        e,
        -e/2.0,
        -5*e/24.0,
        9*e/16.0,
        -2257.0*e/5760.0,
        37.0*e/11520.0,
        24445.0*e/82944.0,
        -mpreal("1989077") * e / mpreal("5806080"),
        mpreal("236533639") * e / mpreal("1393459200"),
        mpreal("22039481") * e / mpreal("309657600"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_chi_ci(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28Chi%28x%29-Ci%28x%29%29%2F%28x%5E2%29+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    y = 1.0/2.0 + x*x*x*x*(1.0/mpreal("2160") + x*x*x*x*1.0/mpreal("18144000"));
    return y;
}

mpreal oracle_fc_bj(double dx) {
    // https://www.wolframalpha.com/input?i=series+1%2FfresnelC%28x%29-cos%28x%29%2Fx+at+x%3D0
    mpreal::set_default_prec(highprec);
    mpreal x = dx, y;
    mpreal pi = mpfr::const_pi();
    int num = 10;
    mpreal c[10] = {
        0,
        mpreal("1") / 2.0,
        0,
        (3.0*pi*pi-5.0) / mpreal("120"),
        0,
        mpreal("1") / 720.0,
        0,
        (203.0*pi*pi*pi*pi-15.0) / mpreal("604800"),
        0,
        mpreal("1") / mpreal("3628800"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_gb_sqrt(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28sqrt%28x%2B1%29-2%29%2F%28x-3%29+at+x%3D3
    mpreal::set_default_prec(highprec);
    mpreal x = dx-3.0, y;
    int num = 10;
    mpreal c[10] = {
        1.0/mpreal("4"),
        -1.0/mpreal("64"),
        1.0/mpreal("512"),
        -5.0/mpreal("16384"),
        7.0/mpreal("131072"),
        -21.0/mpreal("2097152"),
        33.0/mpreal("16777216"),
        -429.0/mpreal("1073741824"),
        715.0/mpreal("8589934592"),
        2431.0/mpreal("137438953472"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_l1_l2(double dx) {
    // https://www.wolframalpha.com/input?i=series+%281%2Blog%28x%29-x%29%2F%281%2F2%28x-1%29%5E2%29+at+x%3D1
    mpreal x = dx-1.0, y;
    int num = 12;
    mpreal c[12] = {
        -1,
        2.0/mpreal("3"),
        -1.0/mpreal("2"),
        2.0/mpreal("5"),
        -1.0/mpreal("3"),
        2.0/mpreal("7"),
        -1.0/mpreal("4"),
        2.0/mpreal("9"),
        -1.0/mpreal("5"),
        2.0/mpreal("11"),
        -1.0/mpreal("6"),
        2.0/mpreal("13"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_hyp_g2(double dx) {
    // https://www.wolframalpha.com/input?i=series+%28sqrt%282%29*sqrt%28x%5E2%2B%281%2Fx%29%5E2%29-2%29+%2F%284*%28x-1%29%5E2%29+at+x%3D1
    mpreal::set_default_prec(highprec);
    mpreal x = dx-1, y;
    mpreal pi = mpfr::const_pi();
    int num = 10;
    mpreal c[10] = {
        1.0/2.0,
        -1.0/mpreal("2"),
        3.0/mpreal("8"),
        -1.0/mpreal("4"),
        1.0/mpreal("4"),
        -3.0/mpreal("8"),
        31.0/mpreal("64"),
        -7.0/mpreal("16"),
        17.0/mpreal("64"),
        -11.0/mpreal("64"),
    };
    y = 0;
    for (int i = num-1; i >= 0; i--) {
        y = y*x + c[i];
    }
    return y;
}

mpreal oracle_sympy_28514_log_ratio(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal b = dx;
    return -log(1 + b);
}

mpreal oracle_sympy_12671_log_cosh(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx;
    mpreal ax = fabs(x);
    return ax - log(2.0) + log(1 + exp(-2.0 * ax));
}

mpreal oracle_statsmodels_1604_logistic(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal eta = dx;
    return 1 / (1 + exp(-eta));
}

mpreal oracle_statsmodels_826_cloglog(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal p = dx;
    return log(-log(1 - p));
}

mpreal oracle_scipy_17194_expm1_kernel(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal g = dx;
    return -(exp(-g) - 1);
}

mpreal oracle_pytorch_39242_log1mexp(double dx) {
    mpreal::set_default_prec(highprec);
    mpreal x = dx;
    return log(1 - exp(x));
}

