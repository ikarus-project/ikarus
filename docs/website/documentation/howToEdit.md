<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer, 
  see [the download page](../../download).

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

## Insert a latex formula
`$$ \mathbf{X} \left( \xi,\eta \right) = \begin{bmatrix} \xi^2 \\ 5\xi\eta \end{bmatrix} $$` 
is compiled to

$$ \mathbf{X} \left( \xi,\eta \right) = \begin{bmatrix} \xi^2 \\ 5\xi\eta \end{bmatrix} $$

## Insert C++ code
The C++ code:
```cpp
double complicatedCalculation(double number, double anotherNumber) 
{
  return number*anotherNumber;
};
```
How it needs to be written in markdown:
```
    ```cpp
    double complicatedCalculation(double number, double anotherNumber) 
    {
      return number*anotherNumber;
    };
    ```
```

## Insert a table
Look at the markdown file of this page to see how a table can be inserted.

| Grid Entity Interface        ||
| :------------ | :-----------: |
| `#!cpp GridViewType leafGridView()`     |
| `#!cpp GridViewType levelGridView(int level)`     |


!!! warning "Insert a warning"
    Note that the **four** spaces at the beginning of this line are essential for the warning to be displayed
    correctly.

!!! note "References"
    For available features in the documentation see [Mkdocs-Material](https://squidfunk.github.io/mkdocs-material/) and [Mkdocs](https://www.mkdocs.org/user-guide/).