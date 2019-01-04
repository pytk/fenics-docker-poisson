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

![poisson](https://user-images.githubusercontent.com/27273842/50679244-0c33fe80-1046-11e9-9e0e-53160b3e77c4.png)

