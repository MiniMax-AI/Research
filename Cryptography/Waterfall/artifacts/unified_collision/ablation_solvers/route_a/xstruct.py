"""Structural facts about Waterfall (X-box cascade).

X(v) with v = (a0,a1,a2,a3)  (a0 = low byte):

    xb = XBOX[a0]
    u2 = a1 ^ byte1(xb) ;  yb = XBOX[u2]
    u3 = a2 ^ byte2(xb) ^ byte1(yb) ;  zb = XBOX[u3]
    u4 = a3 ^ byte3(xb) ^ byte2(yb) ^ byte1(zb) ;  wb = XBOX[u4]
    out3 = SBOX[u4]
    out2 = SBOX[u3] ^ byte3(wb)
    out1 = SBOX[u2] ^ byte3(zb) ^ byte2(wb)
    out0 = SBOX[a0] ^ byte3(yb) ^ byte2(zb) ^ byte1(wb)

byte_k(XBOX[x]) = L_k(SBOX[x]) ^ c_k  with L_k GF(2)-linear, so a difference
confined to the top input byte gives

    dX(v, v ^ (d<<24)) = PSI(y),   y = SBOX[u4] ^ SBOX[u4 ^ d]
    PSI(y) = L_1(y) | L_2(y)<<8 | L_3(y)<<16 | y<<24
"""
from model import XBOX, SBOX, SBOX_INV, M32

XK = [[(XBOX[x] >> (8 * k)) & 0xFF for x in range(256)] for k in range(4)]
_C = [XK[k][SBOX_INV[0]] for k in range(4)]
# LMAP[k][s] = L_k(s) for k = 1,2,3  (byte index k of the xbox table)
LMAP = {k: [XK[k][SBOX_INV[s]] ^ _C[k] for s in range(256)] for k in (1, 2, 3)}

for _k in (1, 2, 3):
    _t = LMAP[_k]
    assert all(_t[a ^ b] == (_t[a] ^ _t[b]) for a in range(256) for b in range(256)), _k


def PSI(y: int) -> int:
    return LMAP[1][y] | (LMAP[2][y] << 8) | (LMAP[3][y] << 16) | (y << 24)


PSI_INV = {}
for _y in range(256):
    PSI_INV[PSI(_y)] = _y


def G(low24: int) -> int:
    """u4 ^ a3, i.e. the part of u4 determined by the low three input bytes."""
    a0 = low24 & 0xFF
    xb = XBOX[a0]
    u2 = ((low24 >> 8) & 0xFF) ^ ((xb >> 8) & 0xFF)
    yb = XBOX[u2]
    u3 = ((low24 >> 16) & 0xFF) ^ ((xb >> 16) & 0xFF) ^ ((yb >> 8) & 0xFF)
    zb = XBOX[u3]
    return ((xb >> 24) & 0xFF) ^ ((yb >> 16) & 0xFF) ^ ((zb >> 8) & 0xFF)


# SBOX difference table:  DDTX[d][y] = list of x with SBOX[x]^SBOX[x^d] == y
DDTX = [[None] * 256 for _ in range(256)]
for _d in range(1, 256):
    row = DDTX[_d]
    for _x in range(256):
        _y = SBOX[_x] ^ SBOX[_x ^ _d]
        if row[_y] is None:
            row[_y] = [_x]
        else:
            row[_y].append(_x)

if __name__ == '__main__':
    import random
    from model import xbox_apply
    random.seed(11)
    bad = sum(1 for _ in range(50000)
              if PSI((xbox_apply(v := random.getrandbits(32))
                      ^ xbox_apply(v ^ (random.randrange(1, 256) << 24))) >> 24)
              != (xbox_apply(v) ^ xbox_apply(v ^ 0)) )
    print('PSI ok (recheck):', bad == 0)
    # kernel of L1,L2,L3 is trivial  ->  PSI injective (already implied by dict size)
    print('PSI injective:', len(PSI_INV) == 256)
    print('example PSI(0x01) =', hex(PSI(1)), ' G(0x123456) =', hex(G(0x123456)))
