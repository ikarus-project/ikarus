<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Installation of Ikarus

!!! warning "Change links on this website when final accounts for repository and docker container are fixed and remove this warning."

Graphical output is currently not supported on Windows 10 (but will probably be available in the future). Therefore,
working on Windows 11 is recommended.

The installations on Windows relies on WSL 2, i.e. although working with Windows, the code is compiled and executed in Linux.

## Installation on Windows using Docker Container
1. Install WSL: Open the PowerShell as an admin and execute the following commands. Reboot afterwards, if requested.
```sh
wsl --install
wsl --set-default-version 2 #(Is not needed for Windows 11)
```
2. [Download and install Docker for Windows](https://docs.docker.com/desktop/windows/install/).
  During the installation, select the option "Install required Windows components for WSL 2"
3. Install debian from [WindowsAppStore](https://www.microsoft.com/en-us/p/debian/9msvkqc78pk6#activetab=pivot:overviewtab)
    1. Open the debian app
    2. Give yourself a username and password
    3. Close the debian app
4. Open the PowerShell and execute:
```sh
wsl --list --all
```
`Debian` should appear as one of the available Linux distributions.
5. In the PowerShell execute:
```sh
wsl --setdefault Debian
```
6. Try to start Docker. If it works, continue with the next step. If a message occurs that you are not allowed to use docker because
    you are not in the docker user group, follow [these instructions](https://icij.gitbook.io/datashare/faq-errors/you-are-not-allowed-to-use-docker-you-must-be-in-the-docker-users-group-.-what-should-i-do).
    In short:
    1. Open **computer management** as admin
    2. Go to **Local users and groups** and find **docker-users**
    3. Add your Account (or a group of which you are a member) to the group.
    4. Restart your computer
7. In Docker, go to Settings --> General and select autostart for docker
    (otherwise you have to start it manually each time you want to work with Ikarus).
8. In the Docker settings, select that Docker uses your WSL2 distribution Debian as shown in the picture.
    ![DockerWslSettings.png](auxiliaryImages/Installation/DockerWslSettings.png)

    In cases docker says that you don't have a WSL 2 distribution, go to the PowerShell and execute
```sh
wsl --set-default-version 2 #(just to be sure that you didn't forgot this at the beginning)
wsl --set-version Debian 2 #(Converts debian to version 2)
```
You should now be able to change the docker settings according to the picture above.

9. Open the PowerShell and execute:
```sh
docker pull rath3t/ikarus-dev:latest #if you want to develop in Ikarus
docker pull rath3t/ikarus:latest #if you want to use Ikarus to run your own main file as in https://github.com/ikarus-project/ikarus-examples
```
6. Download and install [CLion](https://www.jetbrains.com/clion). You need a version >=2022.1.
7. In CLion, go to **File** and **Settings** and apply the following settings for the toolchain:
    ![img.png](auxiliaryImages/Installation/CLionToolchainSettings.png)
    Edit the Container settings and paste the following command into `Run options`:
```
-e DISPLAY=:0 -v \\wsl$\debian\mnt\wslg\.X11-unix:/tmp/.X11-unix -v \\wsl$\debian\mnt\wslg:/mnt/wslg --cap-add=SYS_PTRACE
```
8. [Clone Ikarus](#clone-ikarus)

## Clone Ikarus

- Clone the Ikarus repository as you do it with any other repository (e.g. using GitKraken)
- ToDo: Describe here how to access it from Github.com
- Open the CMake tab `CMake` in the CLion footer:
  ![ClionFooter.png](auxiliaryImages/Installation/ClionFooter.png)
- Click on `Reload CMake project` (refresh symbol)  
  ![ReloadCmakeProject.png](auxiliaryImages/Installation/ReloadCmakeProject.png)
- CMake now detects all required sources automatically. The output should look similar to
  the screenshot below
  ![CMakeOutput.png](auxiliaryImages/Installation/CMakeOutput.png)

## Installation on Windows using WSL
!!! warning "This installation procedure is not recommended"
    The installation using Docker described above has several advantages and should be the standard.
    This section will be removed in the future
1. Install WSL: Open the PowerShell as an admin and execute the following commands. Reboot afterwards, if requested.
  ```sh 
  wsl --install
  wsl --set-default-version 2 #(Is not needed for Windows 11)
  ```
2. Install debian from `.tar` file OR alternatively (3)
  ```sh 
  wsl --import debian <install location> debian_bookworm.tar
  ```
3. Install debian manually
    1. Install from [WindowsAppStore](https://www.microsoft.com/en-us/p/debian/9msvkqc78pk6#activetab=pivot:overviewtab)
    2. Open the debian app
    3. Give yourself a username and password
    4. Execute in debian `sudo sed -i 's/bullseye/bookworm/g' /etc/apt/sources.list`
    5. Execute the following list of commands in debian
  ```sh 
  sudo apt update && \
  sudo apt full-upgrade -y && \
  sudo apt -y install lsb-release && \
  sudo apt -y install build-essential \
  libssl-dev \
  git \
  wget \
  apt-utils \
  software-properties-common \
  gfortran \
  gcc-11 \
  g++-11 \
  gcovr \
  clang \
  libmetis-dev \
  clang-tidy \
  libclang-13-dev \
  clang-format-13 \
  libc++-13-dev \
  libc++abi-13-dev \
  llvm-13-dev \
  liblapack-dev \
  libopenblas-dev \
  libsuitesparse-dev \
  libdune-common-dev \
  libdune-geometry-dev \
  libdune-grid-dev \
  libdune-functions-dev \
  libdune-typetree-dev \
  libdune-localfunctions-dev \
  libdune-uggrid-dev \
  libdune-grid-glue-dev \
  libdune-istl-dev \
  libspdlog-dev \
  libbenchmark-dev \
  libgtest-dev \
  gnuplot \
  python3 \
  pip \
  clang-format-12 \
  gnuplot-x11 \
  curl \
  cppcheck && \
  sudo apt-get install libayatana-appindicator3-1 -y && \
  sudo apt-get -y -f install && \
  sudo apt install libasound2 xvfb -y && \
  wget https://github.com/jgraph/drawio-desktop/releases/download/v16.5.1/drawio-amd64-16.5.1.deb && \
  sudo dpkg -i drawio-amd64-16.5.1.deb && \
  pip install cmakelang==0.6.13 pyyaml && \
  pip install mkdocs && \
  pip install mkdocs-material && \
  pip install mkdocs-macros-plugin && \
  pip install mkdocs-drawio-exporter && \ && \
  sudo cp /usr/bin/clang-format-12 /usr/bin/clang-format && \
  cd /usr/local/bin && \
  sudo ln -s $HOME/.local/bin/cmake-format cmake-format && \
  sudo ln -s $HOME/.local/bin/mkdocs mkdocs && \
  cd ~ && \
  mkdir -p iwyu && \
  cd iwyu && \
  git clone https://github.com/include-what-you-use/include-what-you-use.git && \
  cd include-what-you-use && \
  git checkout clang_13 && \
  cd .. && \
  mkdir -p build && cd build && \
  cmake -G "Unix Makefiles" -DIWYU_LLVM_ROOT_PATH=/usr/lib/llvm-13 ../include-what-you-use && \
  make && \
  sudo make install && \
  cd /usr/src/googletest && \
  cmake . && \
  sudo cmake --build . --target install && \
  cd ~ && \
  git clone https://gitlab.com/libeigen/eigen.git && \
  cd eigen && \
  git checkout 3.4 && \
  mkdir build && \
  cd build && \
  cmake ../ && \
  sudo make install && \
  cd ~ && \
  rm -rf eigen && \
  git clone https://github.com/alandefreitas/matplotplusplus.git && \
  cd matplotplusplus && \
  mkdir -p build && \
  cd build && \
  cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF && \
  cmake --build . --parallel 4 --config Release && \
  sudo cmake --install . && \
  cd ~ && \
  rm -rf matplotplusplus && \
  git clone https://github.com/autodiff/autodiff && \
  cd autodiff/ && \
  mkdir .build && \
  cd .build/ && \
  cmake .. -DAUTODIFF_BUILD_PYTHON=0 -DAUTODIFF_BUILD_EXAMPLES=0 -DAUTODIFF_BUILD_DOCS=0 -DAUTODIFF_BUILD_TESTS=0 && \
  sudo cmake --build . --target install && \
  cd ../.. && \
  mkdir -p dune && \
  cd dune && \
  git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git && \
  git clone https://gitlab.dune-project.org/extensions/dune-foamgrid.git && \
  dunecontrol git checkout releases/2.8 && \
  git clone https://github.com/rath3t/dune-iga.git && \
  dunecontrol cmake "-DCMAKE_BUILD_TYPE=Release" && \
  dunecontrol make && \
  sudo dunecontrol make install && \
  cd .. && \
  rm -rf dune && \
  sudo apt-get auto-remove -y && \
  sudo apt-get clean
  ```
  ```sh
  wget https://raw.githubusercontent.com/JetBrains/clion-wsl/master/ubuntu_setup_env.sh && bash ubuntu_setup_env.sh
  ```
6. In Clion, go to **File** and **Settings** and apply the following settings for the toolchain:
  - Build, Execution, Deployment --> Toolchains: Add with the `+`-sign a WSL configuration
  - Make sure it is used as defaultm i.e. it has to be the first item in the list. Move it up with the arrow buttons otherwise.
