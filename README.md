# Leanna's Software

An open-source home for software, analysis tools, scripts, and small applications created and maintained by Leanna Pfeffer.

Each project lives in its own folder, with its own setup instructions and dependencies. This repository is a collection of independent tools; there is no single installation command for the whole collection.

## Projects

Each project below has its own setup instructions.

| Project | What it does | Requirements | Status |
| --- | --- | --- | --- |
| [Raman raw-WDF analysis and report figures](projects/raman-report/) | PP/PE Raman mapping, spectral analysis, sampling comparisons, and report figures | Python 3.11+ | Validated on the report dataset |

## Add a project

1. Create a folder under `projects/`, using a short name such as `my-analysis-tool`.
2. Put the project's source files in that folder, preserving their existing structure.
3. Copy the [project README template](templates/project/README.md) into the folder as `README.md` and fill in its sections.
4. Include the project's dependency files and a small example when useful.
5. Replace the placeholder row above, or add a new row, with a link to the project's folder.

### Upload through GitHub

Open the `projects` folder, select **Add file > Upload files**, and drag in your project folder. Enter a brief description of the change, then commit it or propose it on a new branch.

GitHub's browser uploader accepts up to 100 files at a time, with a 25 MiB limit per file. See the [official upload instructions](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository) for larger uploads and other options.

### Update a project

Edit or replace files inside the same project folder and commit the changes. GitHub keeps the history, so there is no need to create folders named `final`, `final-v2`, or `latest-copy`.

## Find your way around

| Path | Purpose |
| --- | --- |
| [`projects/`](projects/) | Software projects and their individual instructions |
| [`templates/project/README.md`](templates/project/README.md) | Copyable starting point for each project's documentation |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | How to report problems and contribute changes |
| [`LICENSE`](LICENSE) | MIT license for original code in this repository |

## What to include

Include source code, dependency manifests, instructions, and small shareable examples. Keep credentials, local environment files, and private datasets out of public uploads. The `.gitignore` covers common local files when using Git; browser uploads still need to be checked before submitting.

## License

Original code in this repository is available under the [MIT License](LICENSE), unless a project states otherwise. Keep existing license notices with third-party code and assets. Dependencies retain their own licenses.

