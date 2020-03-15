# Data for the Galaxies project

## Full data

Particle-level data for 3 galaxies at 801 timepoints ('snaps') is supplied by Prof Besla, at various resolutions. The text files are about 133 GB for HighRes. 

This was copied to my local PostgreSQL server, table `galaxy.simdata`, taking up about 200 GB on disk but providing somewhat faster and significantly more flexible data access. Schema-only CREATE commands are supplied in this directory, not the raw data.

## CoM data

Center of mass is calculated for disk (type 2) particles within an appropriate distance of the galactic center. Calculations are in `source/centerofmass.py`. For HighRes data these can be slow, so full results (3 galaxies, all timepoints) were saved to the `com_{galname}.txt` files and in the `galaxy.centerofmass` postgres table.

## Angular momentum data

