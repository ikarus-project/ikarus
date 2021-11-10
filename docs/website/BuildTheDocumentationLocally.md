# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer, 
  see [the installation page](../installation/#clone-ikarus).
- **On Windows:** Make sure [Python](https://www.python.org/downloads/),
  [Drawio](https://github.com/jgraph/drawio-desktop/releases) and [Git](https://git-scm.com/downloads) are installed on 
  your computer. Admin rights are required for the installation, contact your
  admin if it's not installed.
- **On Ubuntu:** execute (one-by-one):
```sh 
sudo apt-get install git -y
sudo apt-get install python3 -y
wget https://github.com/jgraph/drawio-desktop/releases/download/v14.6.13/drawio-amd64-14.6.13.deb
sudo apt-get install libappindicator3-1
sudo dpkg -i drawio-amd64-14.6.13.deb
sudo apt-get -y -f install
sudo apt install libasound2 xvfb
```

## Preview the documentation locally
- Execute (double-click) `Ikarus/docs/BuildLocally/InstalldepsForDocBuild.sh`
- Another window should open now which installs python packages. Wait until it closes itself.
- Changing cmake option: E.g. In Clion: Open `File --> Settings --> Build,Execution,Deployment --> Cmake` 
  Add `-DBUILD_DOCS=TRUE` to your cmake options 
  ![CmakeOptions.png](images/Build documentation locally/CmakeOptions.png)
- Choose target `localSite` and build it (click on the hammer)
  
  ![localSite.png](images/Build documentation locally/localSite.png)
  
- [Click on this link](http://127.0.0.1:8000/)
- Now you should see a live preview of the documentation in your browser
- You can edit the documentation in CLion. `STRG` + `s` saves the documentation and updates it in
  your browser window.
- Cancel the build process to stop the live preview
