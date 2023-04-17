
#include "__avx2_basemul.h"

int16_t const_buff[80] = {
    Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q, Q,
    Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv, Qinv,
    Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar, Qbar,
    sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod,
    sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod, sqrt2Rmod,
    sqrt2invRmod, sqrt2inv, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod,
    sqrt2invRmod, sqrt2inv, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod, sqrt2invRmod
};





