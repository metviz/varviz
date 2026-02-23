# Contributing to VarViz

Thank you for your interest in contributing to VarViz! We welcome contributions from the community — whether it's reporting bugs, suggesting features, improving documentation, or submitting code.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Features](#suggesting-features)
  - [Submitting Code Changes](#submitting-code-changes)
- [Development Setup](#development-setup)
- [Coding Guidelines](#coding-guidelines)
- [Pull Request Process](#pull-request-process)
- [Getting Help](#getting-help)

---

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the maintainers.

---

## How to Contribute

### Reporting Bugs

If you encounter a bug, please [open an issue](https://github.com/metviz/varviz/issues/new) with the following information:

- **Description** — a clear summary of the problem
- **Steps to reproduce** — gene name, variants entered, tracks selected, browser used
- **Expected behavior** — what you expected to happen
- **Actual behavior** — what actually happened (screenshots are very helpful)
- **Environment** — R version, OS, browser, whether running locally or on varviz.org
- **Console output** — any error messages or warnings from the R console (if running locally)

### Suggesting Features

We welcome feature requests! Please [open an issue](https://github.com/metviz/varviz/issues/new) with:

- **Use case** — describe the problem you're trying to solve
- **Proposed solution** — how you envision the feature working
- **Alternatives considered** — other approaches you've thought about
- **Priority** — how important this is for your workflow

### Submitting Code Changes

1. Fork the repository
2. Create a feature branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. Make your changes (see [Coding Guidelines](#coding-guidelines))
4. Test your changes locally with `shiny::runApp()`
5. Commit with clear, descriptive messages:
   ```bash
   git commit -m "Add conservation score tooltip to density plot"
   ```
6. Push to your fork and open a Pull Request

---

## Development Setup

### Prerequisites

- R ≥ 4.2
- RStudio (recommended but not required)

### Install Dependencies

```r
install.packages(c(
  "shiny", "shinyjs", "shinycssloaders", "DT", "plotly",
  "ggplot2", "dplyr", "stringr", "purrr", "data.table",
  "jsonlite", "curl", "httr", "httr2", "readr",
  "cowplot", "grid", "stringi", "wordcloud", "tm"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bedr")
```

### Data Files

Ensure these files are in the `data/` directory:
- `VarViz.RData` — gene coordinate mapping
- `ccrs_b38_ucon.bed.gz` + `.tbi` — constrained coding regions

### Run Locally

```r
shiny::runApp()
```

### Test with a Known Gene

After launching, try:
- **Gene:** `SLC13A5`
- **Variants:** `p.G219R, p.S427L, p.T227M`
- **Tracks:** all enabled

Verify that all eight tracks render correctly and the variant table populates.

---

## Coding Guidelines

### General Principles

- **Readability** — write clear, well-commented R code
- **Consistency** — follow existing code style and naming conventions
- **Modularity** — keep functions focused on single responsibilities
- **Logging** — add `cat()` debug messages with `[Module]` prefixes (e.g., `[ClinVar]`, `[gnomAD]`)

### R Style

- Use `snake_case` for variable and function names
- Use `<-` for assignment (not `=`)
- Indent with 2 spaces
- Keep lines under 120 characters where practical
- Use explicit `package::function()` calls for non-base packages where ambiguity may arise (e.g., `ggplot2::annotate()` vs `NLP::annotate()`)

### API Functions

- All API fetch functions should check the in-memory cache first
- Log cache hits: `cat("[Source] Cache hit for:", key, "\n")`
- Handle errors gracefully — return `NULL` or empty data frames, never crash the app
- Failed API calls should NOT be cached

### Plotting

- All tracks share the same x-axis (amino acid position)
- Use consistent color schemes as defined in existing code
- New tracks should integrate with the plotly subplot system
- Provide both interactive (plotly) and static (ggplot2) versions for download

### UI Changes

- Follow the existing card-based design (rounded corners, `#003f5c` headers)
- Use the established color palette: `#003f5c`, `#58508d`, `#bc5090`, `#ff6361`, `#ffa600`
- Keep the sidebar panel clean — avoid cluttering with too many options

---

## Pull Request Process

1. **Ensure your code works** — test with at least 3 different genes before submitting
2. **Update documentation** — if your change affects user-facing features, update `help.html`
3. **Describe your changes** — write a clear PR description explaining what changed and why
4. **Link related issues** — reference any issues your PR addresses (e.g., "Fixes #42")
5. **One feature per PR** — keep pull requests focused; split large changes into multiple PRs
6. **Be responsive** — address review feedback promptly

### PR Review Criteria

Maintainers will evaluate:
- Does the code work correctly?
- Is it consistent with the existing codebase?
- Does it handle edge cases (missing data, API failures, unusual gene structures)?
- Is the UI impact appropriate?
- Is documentation updated?

---

## Areas Where Help Is Especially Welcome

- **Testing** — automated tests for API functions and plot generation
- **Documentation** — tutorials, use-case examples, video walkthroughs
- **Accessibility** — improving color contrast, screen reader support
- **Performance** — optimizing rendering for genes with many variants
- **New data sources** — integrating additional annotation databases
- **Docker** — containerization and deployment improvements
- **Internationalization** — supporting non-English interfaces

---

## Getting Help

- **Questions about using VarViz:** [open a discussion](https://github.com/metviz/varviz/discussions) or email the corresponding author
- **Bug reports:** [open an issue](https://github.com/metviz/varviz/issues/new)
- **General inquiries:** contact Agasthya Metpally OR Sarathbabu Krishnamurthy (corresponding author)

---

Thank you for helping make VarViz better for the research community!
