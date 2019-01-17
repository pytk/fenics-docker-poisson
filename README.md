# fenics-docker-poisson

 This repository's goal is computing partical differencies equation especially poisson equation. It can work in Docker for Mac like below:

 ```bash
 # mount a volume in local machine to share some file with a docker container
 $ docker run -it -v $(pwd):/home/fenics/shared [fenicspath] --name [name]

 # satart container
 $ docker satart

 # execute container with command line
 $ docker exec -it [name] /bin/bash -l
 ```

### Boundary function for Dirichlet condition

Using FEniCS which is a science compuing library to calculate with fenite element method for solving partial difference equation. The goal is solving self-consistent problem about poisson equation and schrodinger equation in an inversion layer of Germanium channel ultran thin mosfet.


![electrostatic_potential](https://user-images.githubusercontent.com/27273842/51325994-84240f00-1ab1-11e9-8bf8-306328d8e211.png)

we have to consider boudary conditon based on semiconductor physics not only considering electron distribution.

And we assume the device structure is like below

![2019-01-17 1 17 24](https://user-images.githubusercontent.com/27273842/51262424-af492880-19f5-11e9-9aa2-9bbdb9a36334.png)

this project is still going on. Eventually, I'm going to apply to electron transport problem with Monte Carlo Simlation.

