
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

void __barrett_int16x16(int16_t *a){

    int16x16_t v = load_int16x16((int16x16_t*)a);

    v = barrett_int16x16(v);

    store_int16x16((int16x16_t*)a, v);

}





