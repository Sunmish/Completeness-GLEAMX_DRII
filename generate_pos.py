#!/usr/bin/env python

# Generates a list of random RA and Dec positions in specified region of sky

import numpy as np
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS

### gcd() was written by Paul Hancock
# The following functions are explained at http://www.movable-type.co.uk/scripts/latlong.html
# phi ~ lat ~ Dec
# lambda ~ lon ~ RA
def gcd(ra1, dec1, ra2, dec2):
    """
  Great circle distance as calculated by the haversine formula
  ra/dec in degrees
  returns:
  sep in degrees
  """
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(np.radians(dlat) / 2) ** 2
    a += (
        np.cos(np.radians(dec1))
        * np.cos(np.radians(dec2))
        * np.sin(np.radians(dlon) / 2) ** 2
    )
    sep = np.degrees(2 * np.arcsin(min(1, np.sqrt(a))))
    return sep

def get_args():
    # Read input parameters
    description = "Generates a set os source positions with a minimum separation"
    parser = ArgumentParser(description=description)
    parser.add_argument(
        "outfile", type=str, help="Output path to write the source positions to"
    )
    parser.add_argument(
        "--nsrc",
        type=int,
        default=1000,
        help="Number of random source positions to simulate [default=%default]",
    )
    parser.add_argument(
        "--region",
        type=str,
        default="0,360,-90,90",
        help="Region of sky. Enter ra_min, ra_max, dec_min, "
        "dec_max, in deg. For example, for 60 < RA < 300 enter: 60,300,-90,90. For RA < 60 or RA > 300, and -40 < Dec < -10, enter: "
        "300,60,-40,-10. [default=%default]",
    )
    parser.add_argument("--template_image",
        type=str,
        default=None,
        help="Template image to define region, instead of `region`."
    )
    parser.add_argument(
        "--sep-min",
        type=float,
        default=0.0,
        help="Minimum separation between simulated sources, in arcmin "
        "[default=%default]",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=100,
        help="Maximum number of passes to make over generating valid source positions before failing. ",
    )
    parser.add_argument(
        "--return_if_not_converged",
        action="store_true",
        help="Return results even if the placement loop does not converge."
    )
    return parser.parse_args()

def cli(args):
    if args.template_image is None:
        get_pos(
            output_file=args.outfile,
            nsrc=args.nsrc,
            region=args.region,
            sep_min=args.sep_min,
            max_attempts=args.max_attempts,
        )

    else:
        get_pos_from_image(
            output_file=args.outfile,
            nsrc=args.nsrc,
            min_sep=args.sep_min,
            max_attempts=args.max_attempts,
            template_image=args.template_image,
            return_if_not_converged=args.return_if_not_converged
        )



def get_pos_from_image(output_file, nsrc, min_sep, max_attempts, template_image,
    return_if_not_converged=False):
    """
    """

    hdu = fits.open(template_image)
    wcs = WCS(hdu[0].header).celestial

    indices = np.indices((hdu[0].header["NAXIS1"], hdu[0].header["NAXIS2"]))
    
    # im_ra, im_dec = wcs.all_pix2world(indices_y, indices_x, 0)

    source_x = np.zeros(nsrc)
    source_y = np.zeros(nsrc)
    mask = source_x == 0

    success = False
    for c in range(max_attempts):
        source_x[mask] = np.random.uniform(low=1, high=hdu[0].header["NAXIS1"]-1, size=np.sum(mask))
        source_y[mask] = np.random.uniform(low=1, high=hdu[0].header["NAXIS2"]-1, size=np.sum(mask))

        print(source_x)

        source_ra, source_dec = wcs.all_pix2world(source_x, source_y, 0)
        pos = SkyCoord(ra=source_ra, dec=source_dec, unit=(u.deg, u.deg))
        _, seps, _ = match_coordinates_sky(pos, pos, nthneighbor=2)

        mask = seps.value*60. < min_sep
        print(f"Pass {c+1} : Accepted {np.sum(~mask)}")
        if np.all(~mask):
            success = True
            break

    if not success and not return_if_not_converged:
        print("Injection did not converge. Was --sep-min too large?")
        exit(1)


    # jitter positions to avoid source lying directly on a pixel centre:
    # not required because np.random.uniform does not return integers anyway
    if "CDELT2" in hdu[0].header:
        cdelt = hdu[0].header["CDELT2"]
    else:
        cdelt = hdu[0].header["CD2_2"]

    source_ra, source_dec = wcs.all_pix2world(source_x, source_y, 0)
    jitter_ra = np.random.uniform(low=0., high=1., size=len(source_ra))
    jitter_dec = np.random.uniform(low=0., high=1., size=len(source_dec))

    ras = source_ra + jitter_ra*cdelt
    decs = source_dec + jitter_dec*cdelt

    # Print simulated source positions to file
    with open(output_file, "w") as f:
        print("# ra dec", file=f)
        for i in range(0, nsrc):
            print("%.6f" % ras[i], "%.6f" % decs[i], file=f)


def get_pos(output_file, nsrc, region, sep_min, max_attempts, 
):

    # Read region parameter
    ra_min = float(region.split(",")[0])
    ra_max = float(region.split(",")[1])
    dec_min = float(region.split(",")[2])
    dec_max = float(region.split(",")[3])

    # Check RA and Dec ranges
    if (
        dec_min < -90.0
        or dec_min > 90.0
        or dec_max < -90.0
        or dec_max > 90.0
        or dec_min >= dec_max
    ):
        print("Error: Dec range entered incorrectly. Aborting.")
        exit()
    # if ra_min < 0.0 or ra_min > 360.0 or ra_max < 0.0 or ra_max > 360.0 or ra_min == ra_max:
    #     print("Error: RA range entered incorrectly. Aborting.")
    #     exit()
    if ra_max < ra_min:
        ra_max = ra_max + 360.0

    # Convert dec_min and dec_max from deg to rad
    dec_min = np.radians(dec_min)
    dec_max = np.radians(dec_max)

    # Convert sep_min from arcmin to deg
    sep_min = sep_min / 60.0

    ras = np.zeros(nsrc)
    decs = np.zeros(nsrc)
    mask = ras == 0
    min_sep = sep_min * u.arcmin

    c = 0
    success = False
    while (success is False) and (c < max_attempts):
        ras[mask] = np.random.uniform(low=ra_min, high=ra_max, size=np.sum(mask))
        ras = np.where(ras >= 360, ras - 360, ras)

        r = np.random.uniform(low=0.0, high=1.0, size=np.sum(mask))
        decs[mask] = np.rad2deg(
            np.arcsin((np.sin(dec_max) - np.sin(dec_min)) * r + np.sin(dec_min))
        )

        pos = SkyCoord(ras * u.deg, decs * u.deg)
        seps = match_coordinates_sky(pos, pos, nthneighbor=2)

        mask = seps[1] < min_sep
        print(f"Pass {c+1} : Accepted {np.sum(~mask)}")
        if np.all(~mask):
            success = True

        c += 1


    if not success:
        print("Injection did not converge. Was --sep-min to large?")
        exit(1)

    # Print simulated source positions to file
    with open(output_file, "w") as f:
        print("# ra dec", file=f)
        for i in range(0, nsrc):
            print("%.6f" % ras[i], "%.6f" % decs[i], file=f)


if __name__ == "__main__":
    cli(get_args())
