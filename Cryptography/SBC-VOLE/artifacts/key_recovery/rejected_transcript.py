#!/usr/bin/env python3
"""
rejected_transcript.py -- does the key recovery need a VALID signature?

Take one honest signature, corrupt one byte of it, and ask two questions
through the artifact's own code:

  (a) does the artifact's CheckVOLESignature still accept it?
  (b) does the key recovery still return the victim's key?

The attack consumes only Delta*, B, F_x(Delta*) and F_y(Delta*).  Any byte of
the signature that does not feed those four is irrelevant to the attack, and a
signature corrupted there is one the verifier rejects -- so an attacker who
obtains a single *malformed* transcript, one the victim never published and the
verifier would throw away, still recovers the key.

Reports, over every byte offset of the signature:
  usable   = the corrupted signature is REJECTED and the key still recovers
  fragile  = corrupting this byte changes the recovered key
  benign   = the corrupted signature is still accepted (no corruption happened
             in a load-bearing position)

    python3 rejected_transcript.py <driver> [stride]
"""
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recover_key import recover                          # noqa: E402


def main():
    driver = sys.argv[1]
    stride = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    tmp = tempfile.mkdtemp(prefix="rejected_")
    try:
        sig = os.path.join(tmp, "s.bin")
        pk = os.path.join(tmp, "p.bin")
        pub = os.path.join(tmp, "pub.txt")
        msg = "victim message: pay Bob 1 EUR"

        sec = subprocess.run([driver, "keygen", "rejected-transcript-victim"],
                             capture_output=True, text=True).stdout.split()[1:]
        subprocess.run([driver, "sign"] + sec + [msg, sig, pk],
                       capture_output=True, check=True)
        with open(sig, "rb") as f:
            good = bytearray(f.read())
        true = [int(t, 16) for t in sec]

        def attempt(buf):
            csig = os.path.join(tmp, "c.bin")
            with open(csig, "wb") as f:
                f.write(buf)
            if os.path.exists(pub):
                os.remove(pub)
            r = subprocess.run([driver, "verify", msg, csig, pk],
                               capture_output=True, text=True)
            accepted = "FORGE_VERIFY: 1" in r.stdout
            subprocess.run([driver, "extract", msg, csig, pk, pub],
                           capture_output=True, text=True)
            if not os.path.exists(pub):
                return accepted, None
            # a corrupted transcript can leave a large kernel; cap the enumeration
            # so one bad offset cannot dominate the runtime
            words, _, _ = recover(pub, verbose=False, max_cand=8)
            return accepted, words

        acc0, w0 = attempt(good)
        assert acc0 and w0 == true, "baseline failed: %s %s" % (acc0, w0)
        n = len(good)
        usable = fragile = benign = 0
        fragile_pos = []
        usable_pos = []
        for i in range(0, n, stride):
            buf = bytearray(good)
            buf[i] ^= 0x01
            accepted, words = attempt(buf)
            if accepted:
                benign += 1
            elif words == true:
                usable += 1
                usable_pos.append(i)
            else:
                fragile += 1
                if len(fragile_pos) < 24:
                    fragile_pos.append(i)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    tested = usable + fragile + benign
    print("signature length %d bytes, %d byte offsets tested (stride %d)"
          % (n, tested, stride))
    print("  usable  (REJECTED by the artifact's verifier, key still recovered): %d" % usable)
    print("  fragile (corrupting this byte changes the recovered key)          : %d" % fragile)
    print("  benign  (still ACCEPTED -- the byte is not load-bearing)          : %d" % benign)
    if usable_pos:
        print("  usable offsets: %s" % usable_pos)
        print("  -> contiguous run: %s" % (len(usable_pos) > 1 and
              usable_pos[-1] - usable_pos[0] == len(usable_pos) - 1))
    if fragile_pos:
        print("  first fragile offsets: %s" % fragile_pos)
    print()
    if usable:
        print("=> the key recovery does NOT require a valid signature: %d of %d "
              "tested single-byte corruptions are rejected by the artifact's own\n"
              "   verifier and still yield the victim's secret key." % (usable, tested))
    ok = usable > 0
    print("RESULT: %s" % ("SUCCESS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
