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

![2019-01-04 18 40 29](https://user-images.githubusercontent.com/27273842/50681867-499d8980-1050-11e9-9008-99c54d6142d6.png)

 ```python
f = Expression(’x[0]>=0 && x[1]>=0 ? pow(x[0], 2) : 2’, degree=2)
 ```

![poisson](https://user-images.githubusercontent.com/27273842/50679244-0c33fe80-1046-11e9-9e0e-53160b3e77c4.png)

