PHANGS SAMPLE TABLE README

VERSION: v1p5 (first astropy revision)

This table provides survey coordinates, orientations, and some basic
physical parameters for current and possible targets of PHANGS
surveys, including PHANGS-ALMA, PHANGS-HST, and PHANGS-MUSE.

-------------------
FORMAT
-------------------

As of v1p5 (and v2p0) full release, the sample table is constructed by
a network of astropy programs and built as an astropy table. Then it
is written as an ecsv file and a FITS table. The metadata includes
units, though a few unit types may not work with FITS. Most common
astronomy languages include ways to read either FITS binary tables or
ecsv, but if other formats are desired please request them.

-------------------
VERSIONING
-------------------

The sample table is versioned independent of data releases, but we
should be clear about which version of the sample table a data release
uses.

E.g., The version 4 PHANGS-ALMA release will probably the version
"2p0" sample table.

Similar to the ALMA scheme, provisional releases (for comment and with
incremental changes) will be labeled 1p1, 1p2, etc. and major
releases, maintained for backwards compatibility, will be 1p0, 2p0,
etc. It will be expected that v3 or v4, e.g., has more fields than v2.

The first full release on the python-only version will be v2p0.

-------------------
FEEDBACK
-------------------

Right now the overall sample table is maintained by Adam Leroy
(leroy.42@osu.edu). The overall data status page that feeds the
"survey" fields is maintained by I-Ting Ho (iting@mpia.de).

Individual fields and surveys are maintained by individual
contributors. See below or on the Google spreadsheet for the data
status page.

Feedback by email is okay, but opening an issue on the phangs_sample
repository is the ideal way to give feedback. There is also a
sample_table slack channel.

--------------------------------------
DESCRIPTION OF INDIVIDUAL QUANTITIEs
--------------------------------------

The rest of this README includes notes describing individual
quantities or groups of quantities in the sample table.

IDENTIFIERS: 

The table includes columns for a 'name', 'pgc', and a list of
semicolon-delimited aliases ('alias'). The accompanying module can
look up a galaxy in the table given any of these aliases (the
'row_for_gal') routine.

DISTANCE:

The table include distances, which currently descend from the
Extragalactic Distance Database (EDD) through several rounds of
curation and literature input. Uncertainties as recorded in dex (be
careful of this) and a quality (TRGB, QUALITY, etc.) is recorded for
each distance. 

ORIENTATION: 

Position (RA center, Dec center, systemic velocity) and orientation
(position angle and inclination) on the sky of each galaxy. Current
preference is given to rotation curve based orientations (esp. Lang,
Meidt et al. 2020) and near-infrared photometric centers.

ROTATION CURVES

Analytic fit parameters for the rotation curves. These allow a
galactocentric radius to be translated to a predicted circular
rotation. Currently these are based on simple fits to rotation curve
models from Lang, Meidt et al. 2020. The sample_tab_utils module has
functions for expanding these and calculating their derivatives.


INTEGRATED PROPERTIES

Star formation rate, stellar mass, integrated gas masses, and offset
from the main sequence. The integrated HI masses are from the
literature, the H2 masses will be from PHANGS with an aperture
correction but require a little extra work.

The SFR and Mstar values are from "Z0MGS" (Leroy, Sandstrom et
al. resubmitted) based on WISE and GALEX measurements. The advantage
is that these on an understood can be readily compared to large
samples (the full z0mgs or the even larger GSWLC/SDSS sample from
Salim et al.). The disadvantage is that the accuracy of the
measurements may be lower than what we would get from bespoke
analysis.

There is a lot to be written here on details of physical parameter
estimation that goes beyond the scope of a README file.

METALLICITY

Metallicity at the effective radius, r_eff, currently predicted from a
scaling relation (Sanchez et al. 2019) and on the PP04 oxygen scale.

