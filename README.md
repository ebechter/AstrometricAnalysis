# Astrometric Analysis code

This code calculates the astrometric motion of star over a specified time frame using JPL Horizons ephermerides and SIMBAD proper motion and parallax measurements. A directly imaged candidate companion can be determined to be gravitationally bound or an unbound background source in the same line of sight by assessing whether or not the candidate shares common proper motion with the primary star. A background star will remain effectively staionary while the primary star moves across the sky with its known proper motion and parallax. Therefore, it is possible to project the relative, time-varying, relative separation that we expect between the primary and background star, called the "background track" (i.e. the evolution of the companion’s separation as a function of time if it was a background object).        

The primary star's parallactic motion is calculated using its celestial coordinates from SIMBAD and the Earth ephemerides from JPL Horizons. The proper motion is also calculated from SIMBAD. The background track uncertanties are computed using Monte Carlo error propagation that takes into account uncertanties in the primary star’s celestial coordinates, proper motion, and parallax in addition to our measurement uncertainties in astrometric separation.

## Getting Started


### Dependencies 

There are a lot of extraneous files needed to run this code.




### Example Simulation

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
