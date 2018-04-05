# Astrometric Analysis code

This code calculates the astrometric motion of star over a specified time frame using JPL Horizons ephermerides and SIMBAD proper motion and parallax measurements. A directly imaged candidate companion can be determined to be gravitationally bound or an unbound foreground star using this analysis.   

First, multiple relative astrometric values are computed from direct images of the target and companion(s). Then, we propagate the expected motion of the primary star, in time, according to its known proper motion velocity and parallax measurments. This astrometric path, split into $\delta$ East and \Delta North motion, is plotted with 1 and 2 sigma uncertainty regions shaded against the measured relative position of the companion. If the  


Monte Carlo error propagation:
Sample proper motion, parallax, and initial position (in North and East) from a gaussian 
to generate the simulated paths. Use percentile to take 68 and 95 positions on each timestep. 



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
