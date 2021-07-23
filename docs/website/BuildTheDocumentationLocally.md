# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer, 
  see [the installation page](../installation/#clone-ikarus).
- Install [Python](https://www.python.org/downloads/)  
- Install [Drawio](https://github.com/jgraph/drawio-desktop/releases)  
- Install Drawio and other prerequisites on linux:
```sh 
wget https://github.com/jgraph/drawio-desktop/releases/download/v14.6.13/drawio-amd64-14.6.13.deb
sudo apt-get install libappindicator3-1
sudo dpkg -i drawio-amd64-14.6.13.deb
sudo apt-get -y -f install
sudo apt install libasound2 xvfb
```

## Preview the documentation locally
- Execute `docs/BuildLocally/InstalldepsForDocBuild.sh`
- Add `-DBUILD_DOCS=TRUE` to your cmake options
- Build target `localSite`
- Open the link [localhost](http://127.0.0.1:8000/)
- Now you should see a live preview of the docs
- Cancel the build process to stop the live preview
