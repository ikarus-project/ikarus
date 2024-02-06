# How to edit this documentation

## Prerequisites

- Ikarus cloned on your computer,
  see [the download page](../download.md).

## Edit a page

- Open Ikarus
- Go to the folder `docs\website`
- Go to the Markdown file that corresponds to the page to be edited
- Apply your changes using any desired tool
- Create a pull request
- Once the pull request is accepted, the website is automatically updated

## Add a new page

- Open Ikarus
- Go to the folder `docs\website` and create a new Markdown file,
  e.g. `MyAdditionalPage.md`. The new Markdown file could be added in any relevant existing folder or added to a new folder starting
  with a consecutive folder number, e.g., `XX_myFolder`
- Open the file `docs\mkdocs.yml`
- Find the navigation section which starts with `# Navigation`
- The navigation section describes the navigation on the left side of the
website. Add `XX_myFolder/MyAdditionalPage.md` where you want it to appear
- Create a pull request
- Once the pull request is accepted, the website is automatically updated

## Insert a LaTeX formula

The Markdown format:

`$$ \mathbf{X} \left( \xi,\eta \right) = \begin{bmatrix} \xi^2 \\ 5\xi\eta \end{bmatrix} $$`

The compiled output:

$$ \mathbf{X} \left( \xi,\eta \right) = \begin{bmatrix} \xi^2 \\ 5\xi\eta \end{bmatrix} $$

## Insert a C++ code

The Markdown format:

```md
    ```cpp
    double complicatedCalculation(double number, double anotherNumber)
    {
      return number*anotherNumber;
    };
    ```
```

The compiled output:

```cpp
double complicatedCalculation(double number, double anotherNumber)
{
  return number*anotherNumber;
};
```

## Insert tables, warnings and notes

Look at the Markdown file (`03_contribution/howToEdit.md`) to see how tables, warnings and notes can be inserted.

|          Grid Entity Interface                   |
|:------------------------------------------------:|
|       `#!cpp GridViewType leafGridView()`        |
|  `#!cpp GridViewType levelGridView(int level)`   |

!!! warning "Insert a warning"
    Note that the **four** spaces at the beginning of this line are essential for the warning to be displayed
    correctly.

!!! note "References"
    For available features in the documentation see [Mkdocs-Material](https://squidfunk.github.io/mkdocs-material/) and [Mkdocs](https://www.mkdocs.org/user-guide/).
