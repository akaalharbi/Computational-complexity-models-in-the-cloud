L = [ 0x3015e1cc290000,
     0x30256172970000,
     0x3032f6a9e2bd76,
     0x30372e89580000,
     0x3046d988dafac4,
     0x305151413c0000,
     0x30734577870000,
     0x307db4c7b10000,
     0x30935e63440000,
     0x309f8fac6e0000,
     0x30aeb1753b60a5,
     0x30bb901a92efa4,
     0x30c06036c00000,
     0x30e6dd0aa13df6,
     0x70062bbb43f6a0,
     0x7014de347f18c6,
     0x701e395b1ed2de,
     0x701e9f00594ca7,
     0x7025b027066611,
     0x70269fbdd7ab3f,
     0x7038261a0bd34b,
     0x7062487de82bb7,
     0x706e11745d21e8,
     0x70827dd1c6e1b0,
     0x7083aeeccf0000,
     0x708ebc74050000,
     0x708fc77b539d82,
     0x70b1fe63c90000,
     0x70e7e5fbbe0ea8,
     0x70fd93694bb776,
     0xb0006109ddad28,
     0xb012e542bdec53,
     0xb015107ade4232,
     0xb015e4b3a5aa0e,
     0xb0468b38ce0000,
     0xb05b2bb2fa5668,
     0xb05f0b30f00000,
     0xb09bc5ec8f46e4,
     0xb0fcb7fb005e79,
     0xf008e01c620000,
     0xf00ff282177145,
     0xf0123bcd463bb9,
     0xf014ed288e0000,
     0xf021e6c70d040b,
     0xf02d8d35a20000,
     0xf03f28642b3c98,
     0xf040d036840000,
     0xf0426f5d620000,
     0xf056d1fc9cb70d,
     0xf0851fa1d221fe,
      0xf0f4cdd6eaef88,]



good = 0

with open("all_hashes", "r") as f:
    for line in f:
        try:
            line_int = int(line, 0)
#            print(hex(line_int))
            if line_int in L:
                print(f"line_int={line_int} found in L")
        except:
            pass


