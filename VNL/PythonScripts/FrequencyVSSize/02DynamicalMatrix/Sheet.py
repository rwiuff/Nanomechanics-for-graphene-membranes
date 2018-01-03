# Set up lattice
lattice = Hexagonal(24.612*Angstrom, 20.0*Angstrom)

# Define elements
elements = [Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon,
            Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon, Carbon]

# Define coordinates
fractional_coordinates = [[ 0.0333335,  0.0166665,  0.5      ],
                          [ 0.0666665,  0.0833335,  0.5      ],
                          [ 0.0333335,  0.1166665,  0.5      ],
                          [ 0.0666665,  0.1833335,  0.5      ],
                          [ 0.0333335,  0.2166665,  0.5      ],
                          [ 0.0666665,  0.2833335,  0.5      ],
                          [ 0.0333335,  0.3166665,  0.5      ],
                          [ 0.0666665,  0.3833335,  0.5      ],
                          [ 0.0333335,  0.4166665,  0.5      ],
                          [ 0.0666665,  0.4833335,  0.5      ],
                          [ 0.0333335,  0.5166665,  0.5      ],
                          [ 0.0666665,  0.5833335,  0.5      ],
                          [ 0.0333335,  0.6166665,  0.5      ],
                          [ 0.0666665,  0.6833335,  0.5      ],
                          [ 0.0333335,  0.7166665,  0.5      ],
                          [ 0.0666665,  0.7833335,  0.5      ],
                          [ 0.0333335,  0.8166665,  0.5      ],
                          [ 0.0666665,  0.8833335,  0.5      ],
                          [ 0.0333335,  0.9166665,  0.5      ],
                          [ 0.0666665,  0.9833335,  0.5      ],
                          [ 0.1333335,  0.0166665,  0.5      ],
                          [ 0.1666665,  0.0833335,  0.5      ],
                          [ 0.1333335,  0.1166665,  0.5      ],
                          [ 0.1666665,  0.1833335,  0.5      ],
                          [ 0.1333335,  0.2166665,  0.5      ],
                          [ 0.1666665,  0.2833335,  0.5      ],
                          [ 0.1333335,  0.3166665,  0.5      ],
                          [ 0.1666665,  0.3833335,  0.5      ],
                          [ 0.1333335,  0.4166665,  0.5      ],
                          [ 0.1666665,  0.4833335,  0.5      ],
                          [ 0.1333335,  0.5166665,  0.5      ],
                          [ 0.1666665,  0.5833335,  0.5      ],
                          [ 0.1333335,  0.6166665,  0.5      ],
                          [ 0.1666665,  0.6833335,  0.5      ],
                          [ 0.1333335,  0.7166665,  0.5      ],
                          [ 0.1666665,  0.7833335,  0.5      ],
                          [ 0.1333335,  0.8166665,  0.5      ],
                          [ 0.1666665,  0.8833335,  0.5      ],
                          [ 0.1333335,  0.9166665,  0.5      ],
                          [ 0.1666665,  0.9833335,  0.5      ],
                          [ 0.2333335,  0.0166665,  0.5      ],
                          [ 0.2666665,  0.0833335,  0.5      ],
                          [ 0.2333335,  0.1166665,  0.5      ],
                          [ 0.2666665,  0.1833335,  0.5      ],
                          [ 0.2333335,  0.2166665,  0.5      ],
                          [ 0.2666665,  0.2833335,  0.5      ],
                          [ 0.2333335,  0.3166665,  0.5      ],
                          [ 0.2666665,  0.3833335,  0.5      ],
                          [ 0.2333335,  0.4166665,  0.5      ],
                          [ 0.2666665,  0.4833335,  0.5      ],
                          [ 0.2333335,  0.5166665,  0.5      ],
                          [ 0.2666665,  0.5833335,  0.5      ],
                          [ 0.2333335,  0.6166665,  0.5      ],
                          [ 0.2666665,  0.6833335,  0.5      ],
                          [ 0.2333335,  0.7166665,  0.5      ],
                          [ 0.2666665,  0.7833335,  0.5      ],
                          [ 0.2333335,  0.8166665,  0.5      ],
                          [ 0.2666665,  0.8833335,  0.5      ],
                          [ 0.2333335,  0.9166665,  0.5      ],
                          [ 0.2666665,  0.9833335,  0.5      ],
                          [ 0.3333335,  0.0166665,  0.5      ],
                          [ 0.3666665,  0.0833335,  0.5      ],
                          [ 0.3333335,  0.1166665,  0.5      ],
                          [ 0.3666665,  0.1833335,  0.5      ],
                          [ 0.3333335,  0.2166665,  0.5      ],
                          [ 0.3666665,  0.2833335,  0.5      ],
                          [ 0.3333335,  0.3166665,  0.5      ],
                          [ 0.3666665,  0.3833335,  0.5      ],
                          [ 0.3333335,  0.4166665,  0.5      ],
                          [ 0.3666665,  0.4833335,  0.5      ],
                          [ 0.3333335,  0.5166665,  0.5      ],
                          [ 0.3666665,  0.5833335,  0.5      ],
                          [ 0.3333335,  0.6166665,  0.5      ],
                          [ 0.3666665,  0.6833335,  0.5      ],
                          [ 0.3333335,  0.7166665,  0.5      ],
                          [ 0.3666665,  0.7833335,  0.5      ],
                          [ 0.3333335,  0.8166665,  0.5      ],
                          [ 0.3666665,  0.8833335,  0.5      ],
                          [ 0.3333335,  0.9166665,  0.5      ],
                          [ 0.3666665,  0.9833335,  0.5      ],
                          [ 0.4333335,  0.0166665,  0.5      ],
                          [ 0.4666665,  0.0833335,  0.5      ],
                          [ 0.4333335,  0.1166665,  0.5      ],
                          [ 0.4666665,  0.1833335,  0.5      ],
                          [ 0.4333335,  0.2166665,  0.5      ],
                          [ 0.4666665,  0.2833335,  0.5      ],
                          [ 0.4333335,  0.3166665,  0.5      ],
                          [ 0.4666665,  0.3833335,  0.5      ],
                          [ 0.4333335,  0.4166665,  0.5      ],
                          [ 0.4666665,  0.4833335,  0.5      ],
                          [ 0.4333335,  0.5166665,  0.5      ],
                          [ 0.4666665,  0.5833335,  0.5      ],
                          [ 0.4333335,  0.6166665,  0.5      ],
                          [ 0.4666665,  0.6833335,  0.5      ],
                          [ 0.4333335,  0.7166665,  0.5      ],
                          [ 0.4666665,  0.7833335,  0.5      ],
                          [ 0.4333335,  0.8166665,  0.5      ],
                          [ 0.4666665,  0.8833335,  0.5      ],
                          [ 0.4333335,  0.9166665,  0.5      ],
                          [ 0.4666665,  0.9833335,  0.5      ],
                          [ 0.5333335,  0.0166665,  0.5      ],
                          [ 0.5666665,  0.0833335,  0.5      ],
                          [ 0.5333335,  0.1166665,  0.5      ],
                          [ 0.5666665,  0.1833335,  0.5      ],
                          [ 0.5333335,  0.2166665,  0.5      ],
                          [ 0.5666665,  0.2833335,  0.5      ],
                          [ 0.5333335,  0.3166665,  0.5      ],
                          [ 0.5666665,  0.3833335,  0.5      ],
                          [ 0.5333335,  0.4166665,  0.5      ],
                          [ 0.5666665,  0.4833335,  0.5      ],
                          [ 0.5333335,  0.5166665,  0.5      ],
                          [ 0.5666665,  0.5833335,  0.5      ],
                          [ 0.5333335,  0.6166665,  0.5      ],
                          [ 0.5666665,  0.6833335,  0.5      ],
                          [ 0.5333335,  0.7166665,  0.5      ],
                          [ 0.5666665,  0.7833335,  0.5      ],
                          [ 0.5333335,  0.8166665,  0.5      ],
                          [ 0.5666665,  0.8833335,  0.5      ],
                          [ 0.5333335,  0.9166665,  0.5      ],
                          [ 0.5666665,  0.9833335,  0.5      ],
                          [ 0.6333335,  0.0166665,  0.5      ],
                          [ 0.6666665,  0.0833335,  0.5      ],
                          [ 0.6333335,  0.1166665,  0.5      ],
                          [ 0.6666665,  0.1833335,  0.5      ],
                          [ 0.6333335,  0.2166665,  0.5      ],
                          [ 0.6666665,  0.2833335,  0.5      ],
                          [ 0.6333335,  0.3166665,  0.5      ],
                          [ 0.6666665,  0.3833335,  0.5      ],
                          [ 0.6333335,  0.4166665,  0.5      ],
                          [ 0.6666665,  0.4833335,  0.5      ],
                          [ 0.6333335,  0.5166665,  0.5      ],
                          [ 0.6666665,  0.5833335,  0.5      ],
                          [ 0.6333335,  0.6166665,  0.5      ],
                          [ 0.6666665,  0.6833335,  0.5      ],
                          [ 0.6333335,  0.7166665,  0.5      ],
                          [ 0.6666665,  0.7833335,  0.5      ],
                          [ 0.6333335,  0.8166665,  0.5      ],
                          [ 0.6666665,  0.8833335,  0.5      ],
                          [ 0.6333335,  0.9166665,  0.5      ],
                          [ 0.6666665,  0.9833335,  0.5      ],
                          [ 0.7333335,  0.0166665,  0.5      ],
                          [ 0.7666665,  0.0833335,  0.5      ],
                          [ 0.7333335,  0.1166665,  0.5      ],
                          [ 0.7666665,  0.1833335,  0.5      ],
                          [ 0.7333335,  0.2166665,  0.5      ],
                          [ 0.7666665,  0.2833335,  0.5      ],
                          [ 0.7333335,  0.3166665,  0.5      ],
                          [ 0.7666665,  0.3833335,  0.5      ],
                          [ 0.7333335,  0.4166665,  0.5      ],
                          [ 0.7666665,  0.4833335,  0.5      ],
                          [ 0.7333335,  0.5166665,  0.5      ],
                          [ 0.7666665,  0.5833335,  0.5      ],
                          [ 0.7333335,  0.6166665,  0.5      ],
                          [ 0.7666665,  0.6833335,  0.5      ],
                          [ 0.7333335,  0.7166665,  0.5      ],
                          [ 0.7666665,  0.7833335,  0.5      ],
                          [ 0.7333335,  0.8166665,  0.5      ],
                          [ 0.7666665,  0.8833335,  0.5      ],
                          [ 0.7333335,  0.9166665,  0.5      ],
                          [ 0.7666665,  0.9833335,  0.5      ],
                          [ 0.8333335,  0.0166665,  0.5      ],
                          [ 0.8666665,  0.0833335,  0.5      ],
                          [ 0.8333335,  0.1166665,  0.5      ],
                          [ 0.8666665,  0.1833335,  0.5      ],
                          [ 0.8333335,  0.2166665,  0.5      ],
                          [ 0.8666665,  0.2833335,  0.5      ],
                          [ 0.8333335,  0.3166665,  0.5      ],
                          [ 0.8666665,  0.3833335,  0.5      ],
                          [ 0.8333335,  0.4166665,  0.5      ],
                          [ 0.8666665,  0.4833335,  0.5      ],
                          [ 0.8333335,  0.5166665,  0.5      ],
                          [ 0.8666665,  0.5833335,  0.5      ],
                          [ 0.8333335,  0.6166665,  0.5      ],
                          [ 0.8666665,  0.6833335,  0.5      ],
                          [ 0.8333335,  0.7166665,  0.5      ],
                          [ 0.8666665,  0.7833335,  0.5      ],
                          [ 0.8333335,  0.8166665,  0.5      ],
                          [ 0.8666665,  0.8833335,  0.5      ],
                          [ 0.8333335,  0.9166665,  0.5      ],
                          [ 0.8666665,  0.9833335,  0.5      ],
                          [ 0.9333335,  0.0166665,  0.5      ],
                          [ 0.9666665,  0.0833335,  0.5      ],
                          [ 0.9333335,  0.1166665,  0.5      ],
                          [ 0.9666665,  0.1833335,  0.5      ],
                          [ 0.9333335,  0.2166665,  0.5      ],
                          [ 0.9666665,  0.2833335,  0.5      ],
                          [ 0.9333335,  0.3166665,  0.5      ],
                          [ 0.9666665,  0.3833335,  0.5      ],
                          [ 0.9333335,  0.4166665,  0.5      ],
                          [ 0.9666665,  0.4833335,  0.5      ],
                          [ 0.9333335,  0.5166665,  0.5      ],
                          [ 0.9666665,  0.5833335,  0.5      ],
                          [ 0.9333335,  0.6166665,  0.5      ],
                          [ 0.9666665,  0.6833335,  0.5      ],
                          [ 0.9333335,  0.7166665,  0.5      ],
                          [ 0.9666665,  0.7833335,  0.5      ],
                          [ 0.9333335,  0.8166665,  0.5      ],
                          [ 0.9666665,  0.8833335,  0.5      ],
                          [ 0.9333335,  0.9166665,  0.5      ],
                          [ 0.9666665,  0.9833335,  0.5      ]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# Add tags
bulk_configuration.addTags('Hole',      [ 66,  67,  68,  69,  70,  86,  87,  88,  89,  90,
                                          91,  92, 107, 108, 109, 110, 111, 112, 113, 129,
                                         130, 131, 132, 133])
bulk_configuration.addTags('Substrate', [  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
                                          10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
                                          20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
                                          30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
                                          40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
                                          50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
                                          60,  61,  62,  63,  64,  72,  73,  74,  75,  76,
                                          77,  78,  79,  80,  81,  82,  83,  84,  94,  95,
                                          96,  97,  98,  99, 100, 101, 102, 103, 104, 105,
                                         115, 116, 117, 118, 119, 120, 121, 122, 123, 124,
                                         125, 126, 127, 135, 136, 137, 138, 139, 140, 141,
                                         142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                         152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                                         162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
                                         172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
                                         182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
                                         192, 193, 194, 195, 196, 197, 198, 199])