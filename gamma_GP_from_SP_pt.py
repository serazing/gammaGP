import numpy as np

def gamma_GP_from_SP_pt(SP, pt, p, lon, lat):
    """
    Global Polynomial of Neutral Density with respect to Practical
    Salinity and potential temperature

    Calculates the Global Polynomial of Neutral Density gammma_GP using
    an approximate form gamma_poly of Neutral Density on each oceanic basin:
    North Atlantic, South Atlantic, Pacific, Indian and Southern Ocean.
    Each function is a polynomial of 28 terms, function of Practical Salinity,
    potential temperature. The function on the Southern Ocean contains another
    function which is a polynomial of 15 terms of Practical Salinity and potential
    temperature times by pressure and a potential temperature terms which is
    effective close to Antarctica in shallow waters. The polynomials on each ocean
    basins are combined to form the Global Polynomial using tests in lattitude
    and longitude to recognize the basin and weighting functions to make the functions
    zip together where the oceanic basins communicate.

    Parameters
    ----------
    SP : array_like
        Practical Salinity (psu)
    pt : array_like
        Potential temperature (deg C)
    p : array_like
        Sea pressure (dbar)
    lon : array_like
        Longitude
    lat : array_like
        Latitude

    Returns
    -------
    gamma_GP : array_like
     The Global Polynomial of Neutral Density (kg/m^3)
    """

    # Normalization of variables
    # --------------------------
    basins = ['NATL', 'SATL', 'PAC', 'IND', 'SO']



def gamma_G(SP, pt, basin='NATL'):
    """
    Parameters
    ----------
    SP : array_like
        Practical Salinity (psu)
    pt : array_like
        Potential temperature (deg C)
    basin : {'NATL', 'SATL', 'PAC', 'IND', 'SO'}, optional
        Name of the basin. Default is 'NATL'.
    """
    # Normalization of variables
    # --------------------------
    SP = SP / 42.
    pt = pt / 40.
    coeffs = gammaG_coeffs[basin]
    gamma_G = coeffs(1, 3) * np.ones(np.shape(SP))
    for k in np.range(2, len(coeffs)):
        i = coeffs(k, 1);
        j = coeffs(k, 2);
        gamma_G = gamma_G + coeffs(k, 3) * SP ** i * pt ** j
    return gamma_G


def basin_mask():
    io_long = [100, 100, 55, 22, 22, 146, 146, 133.9, 126.94,
               123.62, 120.92, 117.42, 114.11, 107.79, 102.57,
               102.57, 98.79, 100]
    io_lat = [20, 40, 40, 20, -90, -90, -41, -12.48, -8.58, -8.39,
              -8.7, -8.82, - 8.02, -7.04, -3.784, 2.9, 10, 20]
    po_long = [100, 140, 240, 260, 272.59, 276.5, 278.65, 280.73, 295.217,
               290, 300, 294, 290, 146, 146, 133.9, 126.94, 123.62, 120.92,
               117.42, 114.11, 107.79, 102.57, 102.57, 98.79, 100]
    po_lat = [20, 66, 66, 19.55, 13.97, 9.6, 8.1, 9.33, 0, -52, -64.5,
              - 67.5, -90, -90, -41, -12.48, -8.58, -8.39, -8.7, -8.82, -8.02,
              - 7.04, -3.784, 2.9, 10, 20];

    i_inter_indian_pacific = (inpolygon(lon, lat, io_long, io_lat) *
                              inpolygon(lon, lat, po_long, po_lat))
    i_indian = inpolygon(long, lat, io_long, io_lat) - i_inter_indian_pacific;
    i_pacific = inpolygon(long, lat, po_long, po_lat);
    i_atlantic = (1 - i_pacific). * (1 - i_indian);