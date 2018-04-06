# Astrometric Analysis 

This code calculates the astrometric motion of star over a specified time frame using JPL Horizons ephermerides and SIMBAD proper motion and parallax measurements. A directly imaged candidate companion can be determined to be gravitationally bound or an unbound background source in the same line of sight by assessing whether or not the candidate shares common proper motion with the primary star. A background star will remain effectively staionary while the primary star moves across the sky with its known proper motion and parallax. Therefore, it is possible to project the relative, time-varying, relative separation that we expect between the primary and background star, called the "background track" (i.e. the evolution of the companion’s separation as a function of time if it was a background object).        

The primary star's parallactic motion is calculated using its celestial coordinates from SIMBAD and the Earth ephemerides from JPL Horizons. The proper motion is also calculated from SIMBAD. The background track uncertanties are computed using Monte Carlo error propagation that takes into account uncertanties in the primary star’s celestial coordinates, proper motion, and parallax in addition to our measurement uncertainties in astrometric separation.

## Getting Started

This code is written in Python 2, and therefore first requires a copy of that in order to run.

### Package dependencies 

The following Python packages must be installed: 
1. numpy 
2. matplotlib
3. csv
4. astropy
5. jplephem
6. de421

### File dependencies 
The following files must also be located in the main running directory: 
1. de430.bsp 
Note: The JPL planetary ephemerides are saved as files of Chebyshev polynomials fit to the Cartesian positions and velocities of the planets, Sun, and Moon, typically in 32-day intervals. DE430 covers 1549-12-21 to 2650-01-25. Referred to the ICRF version 2.0.
2. starlist.txt - example named CompendiumStars.txt
Note: This is a file of all target stars used in the analysis. It is used as a look up table for SIMBAD values. The file is formatted in the following way (with a few examples): 
```
Name,RA,DEC,PRLX(mas),dPRLX(mas),PM_RA(mas/yr),dPM_RA(mas/yr),PM_DEC(mas/yr),dPM_DEC(mas/yr)
HAT-P-7,19h28m59.3616s,+47d58m10.26s,+03.12, 0.44,-14.8,1.5,8.7,1.4
HAT-P-33,07h32m44.218s,+33d50m06.12s,+02.58, 0.06,-2.2,-4.9, 1.3, 1.0
HD129814,14h44m11.69s,+18d27m43.56s,24.32, 0.60,-70.88, 0.54,-151.63, 0.6 
HD142229,15h53m20.01s,+04d15m11.50s,24.52, 1.25,-23.19, 1.08,+8.48, 1.24 
```
3. A folder with the same name as the star that contains the astrometric measurements in the format:
```
Julian Date, NS(mas),dNS(mas),EW(mas),dEW(mas)
2455803.8,-453.055,0.687,-110.458,0.714
2456522.9,-454.782,0.075,-111.798,0.201
2457213.3,-456.027,0.168,-77.298,0.811
```

### Example run

Set target star and starlist: 

```
target = 'HD224983AB'
params = 'CompendiumStars.txt'
```

Then set the number of points to compute the background track at and the number of tracks for Monte Carlo error propagation:
```
npoints = 1000            
ntracks = 1000
```
The program will automatically search for your target in the paramater file and plot the background track against the measured relative astrometry. 

## Authors

**Eric Bechter** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* This code was originally developed for the TRENDS program, PI Justin Crepp , University of Notre Dame.  
* I would like to acknowledge Henry Ngo for his very helpful explanations on getting started and specifically calculating the projected parallactic motion and integrating the ephemerides into Python. 
