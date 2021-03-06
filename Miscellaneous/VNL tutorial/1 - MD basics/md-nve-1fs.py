# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

# Set up lattice
vector_a = [21.7224, 0.0, 0.0]*Angstrom
vector_b = [0.0, 21.7224, 0.0]*Angstrom
vector_c = [0.0, 0.0, 21.7224]*Angstrom
lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
elements = [Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
            Silicon]

# Define coordinates
fractional_coordinates = [[ 0.    ,  0.    ,  0.    ],
                          [ 0.    ,  0.    ,  0.25  ],
                          [ 0.    ,  0.    ,  0.5   ],
                          [ 0.    ,  0.    ,  0.75  ],
                          [ 0.    ,  0.25  ,  0.    ],
                          [ 0.    ,  0.25  ,  0.25  ],
                          [ 0.    ,  0.25  ,  0.5   ],
                          [ 0.    ,  0.25  ,  0.75  ],
                          [ 0.    ,  0.5   ,  0.    ],
                          [ 0.    ,  0.5   ,  0.25  ],
                          [ 0.    ,  0.5   ,  0.5   ],
                          [ 0.    ,  0.5   ,  0.75  ],
                          [ 0.    ,  0.75  ,  0.    ],
                          [ 0.    ,  0.75  ,  0.25  ],
                          [ 0.    ,  0.75  ,  0.5   ],
                          [ 0.    ,  0.75  ,  0.75  ],
                          [ 0.25  ,  0.    ,  0.    ],
                          [ 0.25  ,  0.    ,  0.25  ],
                          [ 0.25  ,  0.    ,  0.5   ],
                          [ 0.25  ,  0.    ,  0.75  ],
                          [ 0.25  ,  0.25  ,  0.    ],
                          [ 0.25  ,  0.25  ,  0.25  ],
                          [ 0.25  ,  0.25  ,  0.5   ],
                          [ 0.25  ,  0.25  ,  0.75  ],
                          [ 0.25  ,  0.5   ,  0.    ],
                          [ 0.25  ,  0.5   ,  0.25  ],
                          [ 0.25  ,  0.5   ,  0.5   ],
                          [ 0.25  ,  0.5   ,  0.75  ],
                          [ 0.25  ,  0.75  ,  0.    ],
                          [ 0.25  ,  0.75  ,  0.25  ],
                          [ 0.25  ,  0.75  ,  0.5   ],
                          [ 0.25  ,  0.75  ,  0.75  ],
                          [ 0.5   ,  0.    ,  0.    ],
                          [ 0.5   ,  0.    ,  0.25  ],
                          [ 0.5   ,  0.    ,  0.5   ],
                          [ 0.5   ,  0.    ,  0.75  ],
                          [ 0.5   ,  0.25  ,  0.    ],
                          [ 0.5   ,  0.25  ,  0.25  ],
                          [ 0.5   ,  0.25  ,  0.5   ],
                          [ 0.5   ,  0.25  ,  0.75  ],
                          [ 0.5   ,  0.5   ,  0.    ],
                          [ 0.5   ,  0.5   ,  0.25  ],
                          [ 0.5   ,  0.5   ,  0.5   ],
                          [ 0.5   ,  0.5   ,  0.75  ],
                          [ 0.5   ,  0.75  ,  0.    ],
                          [ 0.5   ,  0.75  ,  0.25  ],
                          [ 0.5   ,  0.75  ,  0.5   ],
                          [ 0.5   ,  0.75  ,  0.75  ],
                          [ 0.75  ,  0.    ,  0.    ],
                          [ 0.75  ,  0.    ,  0.25  ],
                          [ 0.75  ,  0.    ,  0.5   ],
                          [ 0.75  ,  0.    ,  0.75  ],
                          [ 0.75  ,  0.25  ,  0.    ],
                          [ 0.75  ,  0.25  ,  0.25  ],
                          [ 0.75  ,  0.25  ,  0.5   ],
                          [ 0.75  ,  0.25  ,  0.75  ],
                          [ 0.75  ,  0.5   ,  0.    ],
                          [ 0.75  ,  0.5   ,  0.25  ],
                          [ 0.75  ,  0.5   ,  0.5   ],
                          [ 0.75  ,  0.5   ,  0.75  ],
                          [ 0.75  ,  0.75  ,  0.    ],
                          [ 0.75  ,  0.75  ,  0.25  ],
                          [ 0.75  ,  0.75  ,  0.5   ],
                          [ 0.75  ,  0.75  ,  0.75  ],
                          [ 0.0625,  0.0625,  0.0625],
                          [ 0.0625,  0.0625,  0.3125],
                          [ 0.0625,  0.0625,  0.5625],
                          [ 0.0625,  0.0625,  0.8125],
                          [ 0.0625,  0.3125,  0.0625],
                          [ 0.0625,  0.3125,  0.3125],
                          [ 0.0625,  0.3125,  0.5625],
                          [ 0.0625,  0.3125,  0.8125],
                          [ 0.0625,  0.5625,  0.0625],
                          [ 0.0625,  0.5625,  0.3125],
                          [ 0.0625,  0.5625,  0.5625],
                          [ 0.0625,  0.5625,  0.8125],
                          [ 0.0625,  0.8125,  0.0625],
                          [ 0.0625,  0.8125,  0.3125],
                          [ 0.0625,  0.8125,  0.5625],
                          [ 0.0625,  0.8125,  0.8125],
                          [ 0.3125,  0.0625,  0.0625],
                          [ 0.3125,  0.0625,  0.3125],
                          [ 0.3125,  0.0625,  0.5625],
                          [ 0.3125,  0.0625,  0.8125],
                          [ 0.3125,  0.3125,  0.0625],
                          [ 0.3125,  0.3125,  0.3125],
                          [ 0.3125,  0.3125,  0.5625],
                          [ 0.3125,  0.3125,  0.8125],
                          [ 0.3125,  0.5625,  0.0625],
                          [ 0.3125,  0.5625,  0.3125],
                          [ 0.3125,  0.5625,  0.5625],
                          [ 0.3125,  0.5625,  0.8125],
                          [ 0.3125,  0.8125,  0.0625],
                          [ 0.3125,  0.8125,  0.3125],
                          [ 0.3125,  0.8125,  0.5625],
                          [ 0.3125,  0.8125,  0.8125],
                          [ 0.5625,  0.0625,  0.0625],
                          [ 0.5625,  0.0625,  0.3125],
                          [ 0.5625,  0.0625,  0.5625],
                          [ 0.5625,  0.0625,  0.8125],
                          [ 0.5625,  0.3125,  0.0625],
                          [ 0.5625,  0.3125,  0.3125],
                          [ 0.5625,  0.3125,  0.5625],
                          [ 0.5625,  0.3125,  0.8125],
                          [ 0.5625,  0.5625,  0.0625],
                          [ 0.5625,  0.5625,  0.3125],
                          [ 0.5625,  0.5625,  0.5625],
                          [ 0.5625,  0.5625,  0.8125],
                          [ 0.5625,  0.8125,  0.0625],
                          [ 0.5625,  0.8125,  0.3125],
                          [ 0.5625,  0.8125,  0.5625],
                          [ 0.5625,  0.8125,  0.8125],
                          [ 0.8125,  0.0625,  0.0625],
                          [ 0.8125,  0.0625,  0.3125],
                          [ 0.8125,  0.0625,  0.5625],
                          [ 0.8125,  0.0625,  0.8125],
                          [ 0.8125,  0.3125,  0.0625],
                          [ 0.8125,  0.3125,  0.3125],
                          [ 0.8125,  0.3125,  0.5625],
                          [ 0.8125,  0.3125,  0.8125],
                          [ 0.8125,  0.5625,  0.0625],
                          [ 0.8125,  0.5625,  0.3125],
                          [ 0.8125,  0.5625,  0.5625],
                          [ 0.8125,  0.5625,  0.8125],
                          [ 0.8125,  0.8125,  0.0625],
                          [ 0.8125,  0.8125,  0.3125],
                          [ 0.8125,  0.8125,  0.5625],
                          [ 0.8125,  0.8125,  0.8125],
                          [ 0.125 ,  0.125 ,  0.    ],
                          [ 0.125 ,  0.125 ,  0.25  ],
                          [ 0.125 ,  0.125 ,  0.5   ],
                          [ 0.125 ,  0.125 ,  0.75  ],
                          [ 0.125 ,  0.375 ,  0.    ],
                          [ 0.125 ,  0.375 ,  0.25  ],
                          [ 0.125 ,  0.375 ,  0.5   ],
                          [ 0.125 ,  0.375 ,  0.75  ],
                          [ 0.125 ,  0.625 ,  0.    ],
                          [ 0.125 ,  0.625 ,  0.25  ],
                          [ 0.125 ,  0.625 ,  0.5   ],
                          [ 0.125 ,  0.625 ,  0.75  ],
                          [ 0.125 ,  0.875 ,  0.    ],
                          [ 0.125 ,  0.875 ,  0.25  ],
                          [ 0.125 ,  0.875 ,  0.5   ],
                          [ 0.125 ,  0.875 ,  0.75  ],
                          [ 0.375 ,  0.125 ,  0.    ],
                          [ 0.375 ,  0.125 ,  0.25  ],
                          [ 0.375 ,  0.125 ,  0.5   ],
                          [ 0.375 ,  0.125 ,  0.75  ],
                          [ 0.375 ,  0.375 ,  0.    ],
                          [ 0.375 ,  0.375 ,  0.25  ],
                          [ 0.375 ,  0.375 ,  0.5   ],
                          [ 0.375 ,  0.375 ,  0.75  ],
                          [ 0.375 ,  0.625 ,  0.    ],
                          [ 0.375 ,  0.625 ,  0.25  ],
                          [ 0.375 ,  0.625 ,  0.5   ],
                          [ 0.375 ,  0.625 ,  0.75  ],
                          [ 0.375 ,  0.875 ,  0.    ],
                          [ 0.375 ,  0.875 ,  0.25  ],
                          [ 0.375 ,  0.875 ,  0.5   ],
                          [ 0.375 ,  0.875 ,  0.75  ],
                          [ 0.625 ,  0.125 ,  0.    ],
                          [ 0.625 ,  0.125 ,  0.25  ],
                          [ 0.625 ,  0.125 ,  0.5   ],
                          [ 0.625 ,  0.125 ,  0.75  ],
                          [ 0.625 ,  0.375 ,  0.    ],
                          [ 0.625 ,  0.375 ,  0.25  ],
                          [ 0.625 ,  0.375 ,  0.5   ],
                          [ 0.625 ,  0.375 ,  0.75  ],
                          [ 0.625 ,  0.625 ,  0.    ],
                          [ 0.625 ,  0.625 ,  0.25  ],
                          [ 0.625 ,  0.625 ,  0.5   ],
                          [ 0.625 ,  0.625 ,  0.75  ],
                          [ 0.625 ,  0.875 ,  0.    ],
                          [ 0.625 ,  0.875 ,  0.25  ],
                          [ 0.625 ,  0.875 ,  0.5   ],
                          [ 0.625 ,  0.875 ,  0.75  ],
                          [ 0.875 ,  0.125 ,  0.    ],
                          [ 0.875 ,  0.125 ,  0.25  ],
                          [ 0.875 ,  0.125 ,  0.5   ],
                          [ 0.875 ,  0.125 ,  0.75  ],
                          [ 0.875 ,  0.375 ,  0.    ],
                          [ 0.875 ,  0.375 ,  0.25  ],
                          [ 0.875 ,  0.375 ,  0.5   ],
                          [ 0.875 ,  0.375 ,  0.75  ],
                          [ 0.875 ,  0.625 ,  0.    ],
                          [ 0.875 ,  0.625 ,  0.25  ],
                          [ 0.875 ,  0.625 ,  0.5   ],
                          [ 0.875 ,  0.625 ,  0.75  ],
                          [ 0.875 ,  0.875 ,  0.    ],
                          [ 0.875 ,  0.875 ,  0.25  ],
                          [ 0.875 ,  0.875 ,  0.5   ],
                          [ 0.875 ,  0.875 ,  0.75  ],
                          [ 0.1875,  0.1875,  0.0625],
                          [ 0.1875,  0.1875,  0.3125],
                          [ 0.1875,  0.1875,  0.5625],
                          [ 0.1875,  0.1875,  0.8125],
                          [ 0.1875,  0.4375,  0.0625],
                          [ 0.1875,  0.4375,  0.3125],
                          [ 0.1875,  0.4375,  0.5625],
                          [ 0.1875,  0.4375,  0.8125],
                          [ 0.1875,  0.6875,  0.0625],
                          [ 0.1875,  0.6875,  0.3125],
                          [ 0.1875,  0.6875,  0.5625],
                          [ 0.1875,  0.6875,  0.8125],
                          [ 0.1875,  0.9375,  0.0625],
                          [ 0.1875,  0.9375,  0.3125],
                          [ 0.1875,  0.9375,  0.5625],
                          [ 0.1875,  0.9375,  0.8125],
                          [ 0.4375,  0.1875,  0.0625],
                          [ 0.4375,  0.1875,  0.3125],
                          [ 0.4375,  0.1875,  0.5625],
                          [ 0.4375,  0.1875,  0.8125],
                          [ 0.4375,  0.4375,  0.0625],
                          [ 0.4375,  0.4375,  0.3125],
                          [ 0.4375,  0.4375,  0.5625],
                          [ 0.4375,  0.4375,  0.8125],
                          [ 0.4375,  0.6875,  0.0625],
                          [ 0.4375,  0.6875,  0.3125],
                          [ 0.4375,  0.6875,  0.5625],
                          [ 0.4375,  0.6875,  0.8125],
                          [ 0.4375,  0.9375,  0.0625],
                          [ 0.4375,  0.9375,  0.3125],
                          [ 0.4375,  0.9375,  0.5625],
                          [ 0.4375,  0.9375,  0.8125],
                          [ 0.6875,  0.1875,  0.0625],
                          [ 0.6875,  0.1875,  0.3125],
                          [ 0.6875,  0.1875,  0.5625],
                          [ 0.6875,  0.1875,  0.8125],
                          [ 0.6875,  0.4375,  0.0625],
                          [ 0.6875,  0.4375,  0.3125],
                          [ 0.6875,  0.4375,  0.5625],
                          [ 0.6875,  0.4375,  0.8125],
                          [ 0.6875,  0.6875,  0.0625],
                          [ 0.6875,  0.6875,  0.3125],
                          [ 0.6875,  0.6875,  0.5625],
                          [ 0.6875,  0.6875,  0.8125],
                          [ 0.6875,  0.9375,  0.0625],
                          [ 0.6875,  0.9375,  0.3125],
                          [ 0.6875,  0.9375,  0.5625],
                          [ 0.6875,  0.9375,  0.8125],
                          [ 0.9375,  0.1875,  0.0625],
                          [ 0.9375,  0.1875,  0.3125],
                          [ 0.9375,  0.1875,  0.5625],
                          [ 0.9375,  0.1875,  0.8125],
                          [ 0.9375,  0.4375,  0.0625],
                          [ 0.9375,  0.4375,  0.3125],
                          [ 0.9375,  0.4375,  0.5625],
                          [ 0.9375,  0.4375,  0.8125],
                          [ 0.9375,  0.6875,  0.0625],
                          [ 0.9375,  0.6875,  0.3125],
                          [ 0.9375,  0.6875,  0.5625],
                          [ 0.9375,  0.6875,  0.8125],
                          [ 0.9375,  0.9375,  0.0625],
                          [ 0.9375,  0.9375,  0.3125],
                          [ 0.9375,  0.9375,  0.5625],
                          [ 0.9375,  0.9375,  0.8125],
                          [ 0.125 ,  0.    ,  0.125 ],
                          [ 0.125 ,  0.    ,  0.375 ],
                          [ 0.125 ,  0.    ,  0.625 ],
                          [ 0.125 ,  0.    ,  0.875 ],
                          [ 0.125 ,  0.25  ,  0.125 ],
                          [ 0.125 ,  0.25  ,  0.375 ],
                          [ 0.125 ,  0.25  ,  0.625 ],
                          [ 0.125 ,  0.25  ,  0.875 ],
                          [ 0.125 ,  0.5   ,  0.125 ],
                          [ 0.125 ,  0.5   ,  0.375 ],
                          [ 0.125 ,  0.5   ,  0.625 ],
                          [ 0.125 ,  0.5   ,  0.875 ],
                          [ 0.125 ,  0.75  ,  0.125 ],
                          [ 0.125 ,  0.75  ,  0.375 ],
                          [ 0.125 ,  0.75  ,  0.625 ],
                          [ 0.125 ,  0.75  ,  0.875 ],
                          [ 0.375 ,  0.    ,  0.125 ],
                          [ 0.375 ,  0.    ,  0.375 ],
                          [ 0.375 ,  0.    ,  0.625 ],
                          [ 0.375 ,  0.    ,  0.875 ],
                          [ 0.375 ,  0.25  ,  0.125 ],
                          [ 0.375 ,  0.25  ,  0.375 ],
                          [ 0.375 ,  0.25  ,  0.625 ],
                          [ 0.375 ,  0.25  ,  0.875 ],
                          [ 0.375 ,  0.5   ,  0.125 ],
                          [ 0.375 ,  0.5   ,  0.375 ],
                          [ 0.375 ,  0.5   ,  0.625 ],
                          [ 0.375 ,  0.5   ,  0.875 ],
                          [ 0.375 ,  0.75  ,  0.125 ],
                          [ 0.375 ,  0.75  ,  0.375 ],
                          [ 0.375 ,  0.75  ,  0.625 ],
                          [ 0.375 ,  0.75  ,  0.875 ],
                          [ 0.625 ,  0.    ,  0.125 ],
                          [ 0.625 ,  0.    ,  0.375 ],
                          [ 0.625 ,  0.    ,  0.625 ],
                          [ 0.625 ,  0.    ,  0.875 ],
                          [ 0.625 ,  0.25  ,  0.125 ],
                          [ 0.625 ,  0.25  ,  0.375 ],
                          [ 0.625 ,  0.25  ,  0.625 ],
                          [ 0.625 ,  0.25  ,  0.875 ],
                          [ 0.625 ,  0.5   ,  0.125 ],
                          [ 0.625 ,  0.5   ,  0.375 ],
                          [ 0.625 ,  0.5   ,  0.625 ],
                          [ 0.625 ,  0.5   ,  0.875 ],
                          [ 0.625 ,  0.75  ,  0.125 ],
                          [ 0.625 ,  0.75  ,  0.375 ],
                          [ 0.625 ,  0.75  ,  0.625 ],
                          [ 0.625 ,  0.75  ,  0.875 ],
                          [ 0.875 ,  0.    ,  0.125 ],
                          [ 0.875 ,  0.    ,  0.375 ],
                          [ 0.875 ,  0.    ,  0.625 ],
                          [ 0.875 ,  0.    ,  0.875 ],
                          [ 0.875 ,  0.25  ,  0.125 ],
                          [ 0.875 ,  0.25  ,  0.375 ],
                          [ 0.875 ,  0.25  ,  0.625 ],
                          [ 0.875 ,  0.25  ,  0.875 ],
                          [ 0.875 ,  0.5   ,  0.125 ],
                          [ 0.875 ,  0.5   ,  0.375 ],
                          [ 0.875 ,  0.5   ,  0.625 ],
                          [ 0.875 ,  0.5   ,  0.875 ],
                          [ 0.875 ,  0.75  ,  0.125 ],
                          [ 0.875 ,  0.75  ,  0.375 ],
                          [ 0.875 ,  0.75  ,  0.625 ],
                          [ 0.875 ,  0.75  ,  0.875 ],
                          [ 0.1875,  0.0625,  0.1875],
                          [ 0.1875,  0.0625,  0.4375],
                          [ 0.1875,  0.0625,  0.6875],
                          [ 0.1875,  0.0625,  0.9375],
                          [ 0.1875,  0.3125,  0.1875],
                          [ 0.1875,  0.3125,  0.4375],
                          [ 0.1875,  0.3125,  0.6875],
                          [ 0.1875,  0.3125,  0.9375],
                          [ 0.1875,  0.5625,  0.1875],
                          [ 0.1875,  0.5625,  0.4375],
                          [ 0.1875,  0.5625,  0.6875],
                          [ 0.1875,  0.5625,  0.9375],
                          [ 0.1875,  0.8125,  0.1875],
                          [ 0.1875,  0.8125,  0.4375],
                          [ 0.1875,  0.8125,  0.6875],
                          [ 0.1875,  0.8125,  0.9375],
                          [ 0.4375,  0.0625,  0.1875],
                          [ 0.4375,  0.0625,  0.4375],
                          [ 0.4375,  0.0625,  0.6875],
                          [ 0.4375,  0.0625,  0.9375],
                          [ 0.4375,  0.3125,  0.1875],
                          [ 0.4375,  0.3125,  0.4375],
                          [ 0.4375,  0.3125,  0.6875],
                          [ 0.4375,  0.3125,  0.9375],
                          [ 0.4375,  0.5625,  0.1875],
                          [ 0.4375,  0.5625,  0.4375],
                          [ 0.4375,  0.5625,  0.6875],
                          [ 0.4375,  0.5625,  0.9375],
                          [ 0.4375,  0.8125,  0.1875],
                          [ 0.4375,  0.8125,  0.4375],
                          [ 0.4375,  0.8125,  0.6875],
                          [ 0.4375,  0.8125,  0.9375],
                          [ 0.6875,  0.0625,  0.1875],
                          [ 0.6875,  0.0625,  0.4375],
                          [ 0.6875,  0.0625,  0.6875],
                          [ 0.6875,  0.0625,  0.9375],
                          [ 0.6875,  0.3125,  0.1875],
                          [ 0.6875,  0.3125,  0.4375],
                          [ 0.6875,  0.3125,  0.6875],
                          [ 0.6875,  0.3125,  0.9375],
                          [ 0.6875,  0.5625,  0.1875],
                          [ 0.6875,  0.5625,  0.4375],
                          [ 0.6875,  0.5625,  0.6875],
                          [ 0.6875,  0.5625,  0.9375],
                          [ 0.6875,  0.8125,  0.1875],
                          [ 0.6875,  0.8125,  0.4375],
                          [ 0.6875,  0.8125,  0.6875],
                          [ 0.6875,  0.8125,  0.9375],
                          [ 0.9375,  0.0625,  0.1875],
                          [ 0.9375,  0.0625,  0.4375],
                          [ 0.9375,  0.0625,  0.6875],
                          [ 0.9375,  0.0625,  0.9375],
                          [ 0.9375,  0.3125,  0.1875],
                          [ 0.9375,  0.3125,  0.4375],
                          [ 0.9375,  0.3125,  0.6875],
                          [ 0.9375,  0.3125,  0.9375],
                          [ 0.9375,  0.5625,  0.1875],
                          [ 0.9375,  0.5625,  0.4375],
                          [ 0.9375,  0.5625,  0.6875],
                          [ 0.9375,  0.5625,  0.9375],
                          [ 0.9375,  0.8125,  0.1875],
                          [ 0.9375,  0.8125,  0.4375],
                          [ 0.9375,  0.8125,  0.6875],
                          [ 0.9375,  0.8125,  0.9375],
                          [ 0.    ,  0.125 ,  0.125 ],
                          [ 0.    ,  0.125 ,  0.375 ],
                          [ 0.    ,  0.125 ,  0.625 ],
                          [ 0.    ,  0.125 ,  0.875 ],
                          [ 0.    ,  0.375 ,  0.125 ],
                          [ 0.    ,  0.375 ,  0.375 ],
                          [ 0.    ,  0.375 ,  0.625 ],
                          [ 0.    ,  0.375 ,  0.875 ],
                          [ 0.    ,  0.625 ,  0.125 ],
                          [ 0.    ,  0.625 ,  0.375 ],
                          [ 0.    ,  0.625 ,  0.625 ],
                          [ 0.    ,  0.625 ,  0.875 ],
                          [ 0.    ,  0.875 ,  0.125 ],
                          [ 0.    ,  0.875 ,  0.375 ],
                          [ 0.    ,  0.875 ,  0.625 ],
                          [ 0.    ,  0.875 ,  0.875 ],
                          [ 0.25  ,  0.125 ,  0.125 ],
                          [ 0.25  ,  0.125 ,  0.375 ],
                          [ 0.25  ,  0.125 ,  0.625 ],
                          [ 0.25  ,  0.125 ,  0.875 ],
                          [ 0.25  ,  0.375 ,  0.125 ],
                          [ 0.25  ,  0.375 ,  0.375 ],
                          [ 0.25  ,  0.375 ,  0.625 ],
                          [ 0.25  ,  0.375 ,  0.875 ],
                          [ 0.25  ,  0.625 ,  0.125 ],
                          [ 0.25  ,  0.625 ,  0.375 ],
                          [ 0.25  ,  0.625 ,  0.625 ],
                          [ 0.25  ,  0.625 ,  0.875 ],
                          [ 0.25  ,  0.875 ,  0.125 ],
                          [ 0.25  ,  0.875 ,  0.375 ],
                          [ 0.25  ,  0.875 ,  0.625 ],
                          [ 0.25  ,  0.875 ,  0.875 ],
                          [ 0.5   ,  0.125 ,  0.125 ],
                          [ 0.5   ,  0.125 ,  0.375 ],
                          [ 0.5   ,  0.125 ,  0.625 ],
                          [ 0.5   ,  0.125 ,  0.875 ],
                          [ 0.5   ,  0.375 ,  0.125 ],
                          [ 0.5   ,  0.375 ,  0.375 ],
                          [ 0.5   ,  0.375 ,  0.625 ],
                          [ 0.5   ,  0.375 ,  0.875 ],
                          [ 0.5   ,  0.625 ,  0.125 ],
                          [ 0.5   ,  0.625 ,  0.375 ],
                          [ 0.5   ,  0.625 ,  0.625 ],
                          [ 0.5   ,  0.625 ,  0.875 ],
                          [ 0.5   ,  0.875 ,  0.125 ],
                          [ 0.5   ,  0.875 ,  0.375 ],
                          [ 0.5   ,  0.875 ,  0.625 ],
                          [ 0.5   ,  0.875 ,  0.875 ],
                          [ 0.75  ,  0.125 ,  0.125 ],
                          [ 0.75  ,  0.125 ,  0.375 ],
                          [ 0.75  ,  0.125 ,  0.625 ],
                          [ 0.75  ,  0.125 ,  0.875 ],
                          [ 0.75  ,  0.375 ,  0.125 ],
                          [ 0.75  ,  0.375 ,  0.375 ],
                          [ 0.75  ,  0.375 ,  0.625 ],
                          [ 0.75  ,  0.375 ,  0.875 ],
                          [ 0.75  ,  0.625 ,  0.125 ],
                          [ 0.75  ,  0.625 ,  0.375 ],
                          [ 0.75  ,  0.625 ,  0.625 ],
                          [ 0.75  ,  0.625 ,  0.875 ],
                          [ 0.75  ,  0.875 ,  0.125 ],
                          [ 0.75  ,  0.875 ,  0.375 ],
                          [ 0.75  ,  0.875 ,  0.625 ],
                          [ 0.75  ,  0.875 ,  0.875 ],
                          [ 0.0625,  0.1875,  0.1875],
                          [ 0.0625,  0.1875,  0.4375],
                          [ 0.0625,  0.1875,  0.6875],
                          [ 0.0625,  0.1875,  0.9375],
                          [ 0.0625,  0.4375,  0.1875],
                          [ 0.0625,  0.4375,  0.4375],
                          [ 0.0625,  0.4375,  0.6875],
                          [ 0.0625,  0.4375,  0.9375],
                          [ 0.0625,  0.6875,  0.1875],
                          [ 0.0625,  0.6875,  0.4375],
                          [ 0.0625,  0.6875,  0.6875],
                          [ 0.0625,  0.6875,  0.9375],
                          [ 0.0625,  0.9375,  0.1875],
                          [ 0.0625,  0.9375,  0.4375],
                          [ 0.0625,  0.9375,  0.6875],
                          [ 0.0625,  0.9375,  0.9375],
                          [ 0.3125,  0.1875,  0.1875],
                          [ 0.3125,  0.1875,  0.4375],
                          [ 0.3125,  0.1875,  0.6875],
                          [ 0.3125,  0.1875,  0.9375],
                          [ 0.3125,  0.4375,  0.1875],
                          [ 0.3125,  0.4375,  0.4375],
                          [ 0.3125,  0.4375,  0.6875],
                          [ 0.3125,  0.4375,  0.9375],
                          [ 0.3125,  0.6875,  0.1875],
                          [ 0.3125,  0.6875,  0.4375],
                          [ 0.3125,  0.6875,  0.6875],
                          [ 0.3125,  0.6875,  0.9375],
                          [ 0.3125,  0.9375,  0.1875],
                          [ 0.3125,  0.9375,  0.4375],
                          [ 0.3125,  0.9375,  0.6875],
                          [ 0.3125,  0.9375,  0.9375],
                          [ 0.5625,  0.1875,  0.1875],
                          [ 0.5625,  0.1875,  0.4375],
                          [ 0.5625,  0.1875,  0.6875],
                          [ 0.5625,  0.1875,  0.9375],
                          [ 0.5625,  0.4375,  0.1875],
                          [ 0.5625,  0.4375,  0.4375],
                          [ 0.5625,  0.4375,  0.6875],
                          [ 0.5625,  0.4375,  0.9375],
                          [ 0.5625,  0.6875,  0.1875],
                          [ 0.5625,  0.6875,  0.4375],
                          [ 0.5625,  0.6875,  0.6875],
                          [ 0.5625,  0.6875,  0.9375],
                          [ 0.5625,  0.9375,  0.1875],
                          [ 0.5625,  0.9375,  0.4375],
                          [ 0.5625,  0.9375,  0.6875],
                          [ 0.5625,  0.9375,  0.9375],
                          [ 0.8125,  0.1875,  0.1875],
                          [ 0.8125,  0.1875,  0.4375],
                          [ 0.8125,  0.1875,  0.6875],
                          [ 0.8125,  0.1875,  0.9375],
                          [ 0.8125,  0.4375,  0.1875],
                          [ 0.8125,  0.4375,  0.4375],
                          [ 0.8125,  0.4375,  0.6875],
                          [ 0.8125,  0.4375,  0.9375],
                          [ 0.8125,  0.6875,  0.1875],
                          [ 0.8125,  0.6875,  0.4375],
                          [ 0.8125,  0.6875,  0.6875],
                          [ 0.8125,  0.6875,  0.9375],
                          [ 0.8125,  0.9375,  0.1875],
                          [ 0.8125,  0.9375,  0.4375],
                          [ 0.8125,  0.9375,  0.6875],
                          [ 0.8125,  0.9375,  0.9375]]

# Set up configuration
bulk_configuration = BulkConfiguration(
    bravais_lattice=lattice,
    elements=elements,
    fractional_coordinates=fractional_coordinates
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------

potentialSet = Tersoff_Si_1988()
calculator = TremoloXCalculator(parameters=potentialSet)
calculator.setVerletListsDelta(0.25*Angstrom)

bulk_configuration.setCalculator(calculator)
bulk_configuration.update()

# -------------------------------------------------------------
# Molecular Dynamics
# -------------------------------------------------------------

initial_velocity = MaxwellBoltzmannDistribution(
    temperature=600.0*Kelvin,
    remove_center_of_mass_momentum=True
)

method = NVEVelocityVerlet(
    time_step=1*femtoSecond,
    initial_velocity=initial_velocity,
)

md_trajectory = MolecularDynamics(
    bulk_configuration,
    constraints=[],
    trajectory_filename='traj-md-nve-1fs.hdf5',
    steps=1000,
    log_interval=5,
    method=method
)

bulk_configuration = md_trajectory.lastImage()
