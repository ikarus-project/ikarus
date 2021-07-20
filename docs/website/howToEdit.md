# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer, 
  see [the installation page](../installation/#clone-ikarus).

## Edit a page
- Open Ikarus in Clion
- go the folder `docs\website`
- find the markdown file which corresponds to the page you want to edit
- Apply your changes (can be done directly in CLion or in external tools e.g. Sublime Text)
- Create a pull request
- Once the pull request is accepted, the website is automatically updated

## Add a new page
- Open Ikarus in Clion
- go the folder `docs\website` and create a new markdown-File, 
  e.g. `MyAdditionalPage.md`
- Open the file `docs\mkdocs.yml`
- Find the navigation section which starts with `# Navigation`
- The navigation section describes the navigation on the left side of the
website. Add `MyAdditionalPage.md` where you want it to appear
- Create a pull request
- Once the pull request is accepted, the website is automatically updated
