#include "..\include\LAMBERBATTIN.h"

void LAMBERTBATTIN(double *ro, double *r, char* dm, double Dtsec, double *vo, double *v) {
    int n = 3;
    double small = 0.000001;
    double mu = 3.986004418e14;   // m3/s2
    double y1 = 0;
    double magr = 0, magro = 0, CosDeltaNu = 0, SinDeltaNu = 0;
    double *rcrossr = NULL, magrcrossr = 0, DNu = 0;
    double RoR = 0, eps = 0, tan2w = 0, rp = 0, L = 0;
    double m = 0, x = 0, xn = 0, chord = 0, s = 0, lim1 = 0;
    int Loops = 0;
    double tempx = 0, Denom = 0, h1 = 0, h2 = 0;
    double b = 0, u = 0, k2 = 0, y = 0, a = 0;
    double arg1 = 0, arg2 = 0, AlpH = 0, BetH = 0, DH = 0;
    double F = 0, GDot = 0, G = 0, Sinv = 0, Cosv = 0;
    double BetE = 0, am = 0, ae = 0, be = 0, tm = 0;
    double AlpE = 0, DE = 0;

    magr = v_norm(r, 3);
    magro = v_norm(ro, 3);
    CosDeltaNu = v_dot(ro, n, r, n) / (magro * magr);
    rcrossr = v_cross(ro, n, r, n);
    magrcrossr = v_norm(rcrossr, n);

    if (strcmp(dm, "pro") == 0)
        SinDeltaNu = magrcrossr / (magro * magr);
    else
        SinDeltaNu = -magrcrossr / (magro * magr);

    DNu = atan2(SinDeltaNu, CosDeltaNu);
    if (DNu < 0.0)
        DNu = 2.0 * M_PI + DNu;

    RoR = magr / magro;
    eps = RoR - 1.0;
    tan2w = 0.25 * eps * eps / (sqrt(RoR) + RoR * (2.0 + sqrt(RoR)));
    rp = sqrt(magro * magr) * (pow(cos(DNu * 0.25), 2) + tan2w);

    if (DNu < M_PI)
        L = (pow(sin(DNu * 0.25), 2) + tan2w) / (pow(sin(DNu * 0.25), 2) + tan2w + cos(DNu * 0.5));
    else
        L = (pow(cos(DNu * 0.25), 2) + tan2w - cos(DNu * 0.5)) / (pow(cos(DNu * 0.25), 2) + tan2w);

    m = mu * Dtsec * Dtsec / (8.0 * rp * rp * rp);
    x = 10.0;
    xn = L;
    chord = sqrt(magro * magro + magr * magr - 2.0 * magro * magr * cos(DNu));
    s = (magro + magr + chord) * 0.5;
    lim1 = sqrt(m / L);

    Loops = 1;
    while (1) {
        x = xn;
        tempx = seebatt(x);
        Denom = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x + tempx * (3.0 + x)));
        h1 = pow((L + x), 2) * (1.0 + 3.0 * x + tempx) * Denom;
        h2 = m * (x - L + tempx) * Denom;

        b = 0.25 * 27.0 * h2 / (pow((1.0 + h1), 3));
        if (b < -1.0)
            xn = 1.0 - 2.0 * L;
        else {
            if (y1 > lim1)
                xn = xn * (lim1 / y1);
            else {
                u = 0.5 * b / (1.0 + sqrt(1.0 + b));
                k2 = seebattk(u);
                y = ((1.0 + h1) / 3.0) * (2.0 + sqrt(1.0 + b) / (1.0 + 2.0 * u * k2 * k2));
                xn = sqrt(pow((1.0 - L) * 0.5, 2) + m / (y * y)) - (1.0 + L) * 0.5;
            }
        }

        Loops = Loops + 1;
        y1 = sqrt(m / ((L + x) * (1.0 + x)));
        if ((fabs(xn - x) < small) && (Loops > 30))
            break;
    }

    a = mu * Dtsec * Dtsec / (16.0 * rp * rp * xn * y * y);
    if (a < -small) {
        arg1 = sqrt(s / (-2.0 * a));
        arg2 = sqrt((s - chord) / (-2.0 * a));
        AlpH = 2.0 * asinh(arg1);
        BetH = 2.0 * asinh(arg2);
        DH = AlpH - BetH;
        F = 1.0 - (a / magro) * (1.0 - cosh(DH));
        GDot = 1.0 - (a / magr) * (1.0 - cosh(DH));
        G = Dtsec - sqrt(-a * a * a / mu) * (sinh(DH) - DH);
    } else if (a > small) {
        arg1 = sqrt(s / (2.0 * a));
        arg2 = sqrt((s - chord) / (2.0 * a));
        Sinv = arg2;
        Cosv = sqrt(1.0 - (magro + magr - chord) / (4.0 * a));
        BetE = 2.0 * acos(Cosv);
        BetE = 2.0 * asin(Sinv);
        if (DNu > M_PI)
            BetE = -BetE;

        Cosv = sqrt(1.0 - s / (2.0 * a));
        Sinv = arg1;
        am = s * 0.5;
        ae = M_PI;
        be = 2.0 * asin(sqrt((s - chord) / s));
        tm = sqrt(am * am * am / mu) * (ae - (be - sin(be)));
        if (Dtsec > tm)
            AlpE = 2.0 * M_PI - 2.0 * asin(Sinv);
        else
            AlpE = 2.0 * asin(Sinv);
        DE = AlpE - BetE;
        F = 1.0 - (a / magro) * (1.0 - cos(DE));
        GDot = 1.0 - (a / magr) * (1.0 - cos(DE));
        G = Dtsec - sqrt(a * a * a / mu) * (DE - sin(DE));
    } else {
        printf("a parabolic orbit");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < 3; i++) {
        vo[i] = (r[i] - F * ro[i]) / G;
        v[i] = (GDot * r[i] - ro[i]) / G;
    }
}
