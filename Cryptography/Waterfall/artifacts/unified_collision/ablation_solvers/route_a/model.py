"""Python model of the Waterfall hash function."""
from __future__ import annotations

SBOX = [
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16,
]

M32 = 0xFFFFFFFF


def _build_xbox():
    xb = [0] * 256
    for i in range(256):
        work1 = SBOX[i]
        v = (work1 ^ i) & 0xFF
        w = (work1 ^ (work1 << 2) ^ (work1 >> 6) ^ (work1 << 3)) & 0xFF
        w = (((w << 7) & 0xFF) | (((w >> 1) ^ 0x3b) & 0xFF)) & 0xFF
        v |= w << 8
        w = (work1 ^ (work1 << 5) ^ (work1 >> 3) ^ (work1 >> 7) ^ 0x95) & 0xFF
        v |= w << 16
        w = (work1 ^ (work1 << 7) ^ (work1 >> 1) ^ (work1 >> 5) ^ 0x6a) & 0xFF
        v |= w << 24
        xb[i] = v
    return xb


XBOX = _build_xbox()
SBOX_INV = [0] * 256
for _i, _v in enumerate(SBOX):
    SBOX_INV[_v] = _i

STREAM1, STREAM2, STREAM3, POOL = 16, 7, 6, 32


def xbox_apply(x: int) -> int:
    for _ in range(4):
        x ^= XBOX[x & 0xFF]
        x = ((x << 24) | (x >> 8)) & M32
    return x


def update_step(st, inp):
    w1 = st['s1'][st['i1']]
    w2 = st['s2'][st['i2']]
    w3 = st['s3'][st['i3']]
    st['i1'] = (st['i1'] + 1) % STREAM1
    st['i2'] = (st['i2'] + 1) % STREAM2
    st['i3'] = (st['i3'] + 1) % STREAM3
    w1 = xbox_apply(w1 ^ st['s1'][st['i1']] ^ inp)
    st['s1'][st['i1']] = w1
    w2 = xbox_apply(w2 ^ st['s2'][st['i2']] ^ inp ^ 1)
    st['s2'][st['i2']] = w2
    w3 = xbox_apply(w3 ^ st['s3'][st['i3']] ^ inp ^ 2)
    st['s3'][st['i3']] = w3
    st['ip'] = (st['ip'] + 1) % POOL
    st['p2'][st['ip']] ^= w2
    st['p3'][st['ip']] ^= w3


def final_step(work_in, pool):
    out = [0] * STREAM1
    for i in range(STREAM1):
        out[i] = ((work_in[i] & 0x14141414)
                  | (work_in[(i + 3) % STREAM1] & 0x28282828)
                  | (work_in[(i + 7) % STREAM1] & 0x42424242)
                  | (work_in[(i + 12) % STREAM1] & 0x81818181))
    for i in range(STREAM1):
        out[i] = xbox_apply(out[i])
    if pool is not None:
        for i in range(STREAM1):
            out[i] ^= pool[i]
    return out


def new_state(hashbitlen=256):
    return {
        'hashbitlen': hashbitlen,
        'buf': bytearray(64),
        'inputbits': 0,
        'inputblocks': 0,
        'inputslabs': 0,
        's1': [0] * STREAM1, 'i1': STREAM1 - 1,
        's2': [hashbitlen] + [0] * (STREAM2 - 1), 'i2': STREAM2 - 1,
        's3': [0] * STREAM3, 'i3': STREAM3 - 1,
        'p2': [0] * POOL, 'p3': [0] * POOL, 'ip': POOL - 1,
    }


def block_words(buf):
    return [buf[4 * j] | (buf[4 * j + 1] << 8) | (buf[4 * j + 2] << 16) | (buf[4 * j + 3] << 24)
            for j in range(16)]


def update(st, data: bytes):
    i = st['inputbits'] // 8
    pos = 0
    n = len(data)
    while pos < n:
        st['inputbits'] += 8
        st['buf'][i] = data[pos]
        i += 1
        pos += 1
        if i == 64:
            for w in block_words(st['buf']):
                update_step(st, w)
            st['inputblocks'] = (st['inputblocks'] + 1) & M32
            if st['inputblocks'] == 0:
                st['inputslabs'] = (st['inputslabs'] + 1) & M32
            update_step(st, st['inputblocks'])
            i = 0
            st['inputbits'] = 0


def final(st):
    if st['inputbits']:
        for i in range((st['inputbits'] + 7) // 8, 64):
            st['buf'][i] = 0
        st['inputblocks'] = (st['inputblocks'] + 1) & M32
        if st['inputblocks'] == 0:
            st['inputslabs'] = (st['inputslabs'] + 1) & M32
        for w in block_words(st['buf']):
            update_step(st, w)
        update_step(st, st['inputblocks'])
    update_step(st, st['inputbits'])
    update_step(st, st['inputslabs'])
    for i in range(16):
        update_step(st, i)

    w = list(st['s1'])
    pp = st['p2']
    for _ in range(2):
        w = final_step(w, pp[:16])
        pp = pp[16:]
    pp = st['p3']
    for _ in range(2):
        w = final_step(w, pp[:16])
        pp = pp[16:]
    for _ in range(4):
        w = final_step(w, None)
    for i in range(STREAM1):
        w[i] ^= st['s1'][i]
    for i in range(STREAM1):
        w[i] ^= st['p2'][i] ^ st['p2'][i + 16]
    for i in range(STREAM1):
        w[i] ^= st['p3'][i] ^ st['p3'][i + 16]
    out = bytearray()
    for i in range(st['hashbitlen'] // 32):
        out += w[i].to_bytes(4, 'little')
    return bytes(out)


def digest(data: bytes, hashbitlen: int = 256) -> bytes:
    st = new_state(hashbitlen)
    update(st, data)
    return final(st)


if __name__ == '__main__':
    expected = {
        b'': '3a9659246e3c9a44ccfcceea81969fab9c99269eac6b01833088c20189fcde33',
        b'abc': '32b8b5871b24d1b48b74fcc1d5bd49765eb37851602f4a116b6467f8f8f2a46a',
    }
    for msg, exp_hex in expected.items():
        got = digest(msg, 256).hex()
        assert got == exp_hex, f'mismatch on {msg!r}: {got} != {exp_hex}'
        print(f'len={len(msg):4d} match=True {got[:32]}')
    print('ALL OK')
