# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer,
  see [the download page](../../download).

## Preview the documentation locally
- Changing cmake option: E.g. In Clion: Open `File --> Settings --> Build,Execution,Deployment --> Cmake` 
  Add `-DBUILD_DOCS=TRUE` to your cmake options 
  ![CmakeOptions.png](../images/Build%20documentation%20locally/CmakeOptions.png)
- Choose target `localSite` and build it (click on the hammer)
  
  ![localSite.png](../images/Build%20documentation%20locally/localSite.png)
  
- [Click on this link](http://127.0.0.1:8000/)
- Now you should see a live preview of the documentation in your browser
- You can edit the documentation in CLion. `STRG` + `s` saves the documentation and updates it in
  your browser window.
- Cancel the build process to stop the live preview
