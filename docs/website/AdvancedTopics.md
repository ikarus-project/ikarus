<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Advanced topics
The following topics are probably not relevant for your daily work, but you might need this knowledge at some point.

## Add additional dependencies or update existing dependencies to a newer version
All the dependencies of Ikarus are shipped in a docker container. The corresponding docker image is available at 
[https://github.com/ikarus-project/ikarus-docker-container](https://github.com/ikarus-project/ikarus-docker-container).
In order to modify the dependencies, you need to create a modified version of this dockerfile, create a new docker container
and execute your code in this docker container (i.e. change the docker image to be used in the CLion settings).
